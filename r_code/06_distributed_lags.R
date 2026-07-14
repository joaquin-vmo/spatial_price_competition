
# objective: distributed-lag ols estimate of retail-to-wholesale (mepco)
# cost pass-through, in logs, and how the cumulative pass-through differs by
# brand and by n of competitors within 2km.
#
# this is the single-regression, distributed-lag counterpart of the
# horizon-by-horizon local projections in 05_local_projection.R: instead of
# one regression per horizon, one regression per group with K shock lags,
# where the cumulative pass-through at horizon K is PT(K) = sum_{k<=K} beta_k.
# the grid construction below is identical to 05_local_projection.R
#
# the data pools every mepco update since 2014-08-07 as a shock

library(data.table)
library(fixest)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))
mepco <- fread(file.path("data/processed", "whosale_prices_mepco.csv"))

K_DL <- 26 # max lag for the distributed-lag pass-through specification below
fuels <- c("93", "97", "di") # 95 has no wholesale reference

dbf <- db[fuel %in% fuels & !is.na(wholesale_w)]

#####

# resample each station-fuel's price series onto the regular mepco weekly
# calendar (last observation carried forward), so horizon-h differences are
# well defined even though raw reporting is irregular

mepco_dates <- sort(unique(as.Date(mepco$date)))

# dates of each fuel first and last observations by station
id_fuel_bounds <- dbf[,
                      .(min_date = min(date), max_date = max(date)),
                      by = .(id, fuel)]
# grid for the date of each fuel and station's observation
grid <- id_fuel_bounds[,
                       .(date = mepco_dates[mepco_dates >= min_date & mepco_dates <= max_date]),
                        by = .(id, fuel)]

price_obs <- dbf[, .(id, fuel, date, price)]
setkey(price_obs, id, fuel, date)
setkey(grid, id, fuel, date)
grid <- price_obs[grid, roll = TRUE]

# now the same process applied to mepco

# wholesale price and its log-change, straight from mepco (one row per fuel x date)
mepco_long <- rbindlist(lapply(fuels, function(f) {
  data.table(fuel = f, date = mepco$date, wholesale_w = mepco[[paste0(f, "_w/")]])
}))
setorder(mepco_long, fuel, date)
mepco_long[, dlog_wholesale := log(wholesale_w) - shift(log(wholesale_w)), by = fuel]

# the distributed-lag spec needs the shock's own lags k = 1..K_DL as
# regressors (the k = 0 term is dlog_wholesale itself)
for (l in 1:K_DL) {
  mepco_long[, paste0("dlog_wholesale_l", l) := shift(dlog_wholesale, l), by = fuel]
}

grid <- merge(grid, mepco_long, by = c("fuel", "date"), all.x = TRUE)

# n° competitors within 2km, resampled the same way (locf of the station's
# own already-time-varying competitors_2km, computed in 04_build_database.R)
comp_hist <- unique(dbf[, .(id, date, competitors_2km)])
setkey(comp_hist, id, date)
setkey(grid, id, date)
grid <- comp_hist[grid, roll = TRUE]

# brand: modal distributor over the station's history, same big-3-vs-rest
# split used in 05_passthrough_analysis.R
brand <- dbf[, .N, by = .(id, distributor)]
setorder(brand, id, -N)
brand <- brand[, .SD[1], by = id]
brand[, brand_group := fifelse(distributor %in% c("copec", "shell", "aramco_petrobras"),
                                distributor, "independiente")]
grid <- merge(grid, brand[, .(id, brand_group)], by = "id", all.x = TRUE)

comp_breaks <- quantile(grid$competitors_2km, probs = 0:3 / 3, na.rm = TRUE)
grid[, comp_tercile := cut(competitors_2km, breaks = unique(comp_breaks),
                            include.lowest = TRUE, labels = c("baja", "media", "alta"))]

# highway/rest-stop proxy: showers or a truck-specific pump are amenities a
# corner urban station wouldn't bother with, but a long-haul route stop
# would. station-level attributes, so not fe-collinearity relevant, but
# used as an extra heterogeneity split (station fe already absorbs any
# purely time-invariant characteristic, so this isn't a "control" against
# confounding -- it's asking whether the *elasticity itself* differs)
carretera <- unique(dbf[, .(id, duchas, surtidor_camiones)])
carretera <- carretera[, .SD[1], by = id]
carretera[, es_carretera := fifelse(duchas | surtidor_camiones, "carretera", "no_carretera")]
grid <- merge(grid, carretera[, .(id, es_carretera)], by = "id", all.x = TRUE)

# macrozone: 16 regions collapsed to 3, since several (arica, aysen,
# magallanes) have too few stations for a standalone estimate
macrozona_map <- c(
  arica = "norte", tarapaca = "norte", antofagasta = "norte",
  atacama = "norte", coquimbo = "norte",
  valparaiso = "centro", metropolitana = "centro", ohiggins = "centro",
  maule = "centro", nuble = "centro", biobio = "centro",
  araucania = "sur", los_rios = "sur", los_lagos = "sur",
  aysen = "sur", magallanes_antartida = "sur"
)
region_static <- dbf[, .N, by = .(id, region)]
setorder(region_static, id, -N)
region_static <- region_static[, .SD[1], by = id]
region_static[, macrozona := macrozona_map[region]]
grid <- merge(grid, region_static[, .(id, macrozona)], by = "id", all.x = TRUE)

#####

# weekly log price change, the dependent variable for the distributed-lag
# spec (unlike the cumulative y_h used in 05_local_projection.R, this is not
# cumulative -- each beta_k already carries its own horizon)

setorder(grid, id, fuel, date)
grid[, log_price := log(price)]
grid[, dlog_price := log_price - shift(log_price, 1), by = .(id, fuel)]

#####

# one distributed-lag ols per group,
# dlog_price_it ~ sum_{k=0}^{K} beta_k dlog_wholesale_{t-k} | id^fuel,
# with the cumulative pass-through at horizon K given by PT(K) = sum_{k<=K} beta_k.
# se(PT(K)) uses the model's clustered vcov (not a naive sum of se_k), since the
# beta_k are correlated with each other.
#
# tried adding lags of dlog_price itself (ardl) as controls for the price's
# own momentum, but dlog_price_{t-j} is largely generated by
# dlog_wholesale_{t-j-k} that's already in the model -- the resulting
# collinearity blew up sum(beta) to >2x pass-through, offset by a similarly
# large negative sum(gamma) on the price lags (sum(beta)/(1-sum(gamma)) ~
# 0.76 pooled, matching this simpler spec almost exactly). not worth the
# added complexity for the same answer -- dropped
run_dl <- function(data, group_col = NULL, K = K_DL) {
  if (!is.null(group_col)) data <- data[!is.na(get(group_col))]
  groups <- if (is.null(group_col)) list(pooled = data) else split(data, data[[group_col]])
  lag_names <- c("dlog_wholesale", paste0("dlog_wholesale_l", 1:K))
  rbindlist(lapply(names(groups), function(g) {
    d <- groups[[g]]
    fml <- paste0("dlog_price ~ ", paste(lag_names, collapse = " + "), " | id^fuel")
    m <- feols(as.formula(fml), data = d, cluster = ~id)
    b <- coef(m)[lag_names]
    V <- vcov(m)[lag_names, lag_names]
    # cumulative pass-through PT(k) and its se, for k = 0..K
    rbindlist(lapply(0:K, function(k) {
      idx <- lag_names[1:(k + 1)]
      w <- rep(1, k + 1)
      data.table(
        horizon = k, group = g,
        estimate = sum(b[idx]),
        se = sqrt(as.numeric(t(w) %*% V[idx, idx] %*% w))
      )
    }))
  }))
}

irf_dl_pooled <- run_dl(grid)[, dimension := "pooled"]
irf_dl_fuel <- run_dl(grid, "fuel")[, dimension := "combustible"]
irf_dl_brand <- run_dl(grid, "brand_group")[, dimension := "marca"]
irf_dl_comp <- run_dl(grid, "comp_tercile")[, dimension := "competidores_2km"]
irf_dl_carretera <- run_dl(grid, "es_carretera")[, dimension := "carretera"]
irf_dl_macrozona <- run_dl(grid, "macrozona")[, dimension := "macrozona"]

irf_dl <- rbind(irf_dl_pooled, irf_dl_fuel, irf_dl_brand, irf_dl_comp, irf_dl_carretera, irf_dl_macrozona)
irf_dl[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

fwrite(irf_dl, file.path("results/tables", "local_projections_dl_irf.csv"))

#####

# graphs

dl_y_label <- expression(PT(K) * ": " * sum(beta[k], k == 0, K) * ", " *
                            Delta * log(precio)[t] * " por " * Delta * log(mayorista)[t-k])

p_dl_brand <- ggplot(irf_dl[dimension == "marca"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "K: semanas de rezagos del shock incluidas",
    y = dl_y_label,
    color = "marca", fill = "marca",
    title = "distributed lag: pass-through de costo a precio, por marca",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "dl_irf_brand.pdf"), p_dl_brand, width = 7.5, height = 4.5)

p_dl_fuel <- ggplot(irf_dl[dimension == "combustible"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "K: semanas de rezagos del shock incluidas",
    y = dl_y_label,
    color = "combustible", fill = "combustible",
    title = "distributed lag: pass-through de costo a precio, completamente pooled por combustible",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "dl_irf_fuel.pdf"), p_dl_fuel, width = 7.5, height = 4.5)

p_dl_comp <- ggplot(irf_dl[dimension == "competidores_2km"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "K: semanas de rezagos del shock incluidas",
    y = dl_y_label,
    color = "competidores\na 2km (tercil)", fill = "competidores\na 2km (tercil)",
    title = "distributed lag: pass-through de costo a precio, por competencia",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "dl_irf_competitors.pdf"), p_dl_comp, width = 7.5, height = 4.5)

p_dl_carretera <- ggplot(irf_dl[dimension == "carretera"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "K: semanas de rezagos del shock incluidas",
    y = dl_y_label,
    color = "tipo de estacion", fill = "tipo de estacion",
    title = "distributed lag: pass-through de costo a precio, carretera vs no",
    subtitle = "carretera = tiene duchas o surtidor de camiones; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "dl_irf_carretera.pdf"), p_dl_carretera, width = 7.5, height = 4.5)

p_dl_macrozona <- ggplot(irf_dl[dimension == "macrozona"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "K: semanas de rezagos del shock incluidas",
    y = dl_y_label,
    color = "macrozona", fill = "macrozona",
    title = "distributed lag: pass-through de costo a precio, por macrozona",
    subtitle = "norte: arica-coquimbo; centro: valparaiso-biobio; sur: araucania-magallanes",
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "dl_irf_macrozona.pdf"), p_dl_macrozona, width = 7.5, height = 4.5)
