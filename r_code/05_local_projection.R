
# objective: local projections (jorda 2005) impulse-response of retail
# prices to wholesale (mepco) cost changes, in logs, and how the response
# differs by brand and by n of competitors within 2km.
#
# the data pools every mepco update since 2014-08-07 as a shock
#
# this file builds grid so to estimate ols at every period with a function.
# the distributed-lag counterpart of this analysis lives in
# 06_distributed_lags.R and reuses the same grid construction

library(data.table)
library(fixest)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))
mepco <- fread(file.path("data/processed", "whosale_prices_mepco.csv"))

H <- 12 # horizons, in mepco weeks (~3 months)
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

# dlog_wholesale is likely serially correlated (mepco tends to move several
# weeks in a row in the same direction), so a horizon-h regression on the
# contemporaneous shock alone partly picks up the response to later,
# correlated shocks too -- not a clean single-shock irf. following jorda
# (2005), control for n_lags of the shock in every lp regression below, so
# beta_h is purged of the part of the shock predictable from its own recent
# history. lag order picked once by aic on an ar(p) fit per fuel, so it's
# data-driven rather than arbitrary
ar_orders <- sapply(fuels, function(f) {
  s <- na.omit(mepco_long[fuel == f, dlog_wholesale])
  ar(s, aic = TRUE, order.max = 20)$order
})
n_lags <- max(ar_orders, 1)

for (l in 1:n_lags) {
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
# magallanes) have too few stations for a standalone local projection
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

# cumulative log price response at each horizon: y_h = log(p_{t+h}) - log(p_{t-1})

setorder(grid, id, fuel, date)
grid[, log_price := log(price)]
for (h in 0:H) {
  grid[, paste0("y_h", h) := shift(log_price, -h) - shift(log_price, 1), by = .(id, fuel)]
}

#####

# run one local projection per horizon, for the pooled sample and for each
# brand / competition subgroup separately (id^fuel fe only)

# for the data (grid) i and group (brand, competition, etc) i

run_lp <- function(data, group_col = NULL, n_lags = 0) {
  # create groups
  if (!is.null(group_col)) data <- data[!is.na(get(group_col))]
  groups <- if (is.null(group_col)) list(pooled = data) else split(data, data[[group_col]])
  # shock lag controls, so beta_h isn't contaminated by the predictable
  # (autocorrelated) part of dlog_wholesale -- see note above n_lags
  lag_controls <- if (n_lags > 0) {
    paste0(" + ", paste(paste0("dlog_wholesale_l", 1:n_lags), collapse = " + "))
  } else ""
  # loops for group selection (g)
  rbindlist(lapply(names(groups), function(g) {
    d <- groups[[g]]
    # group for horizon selection (h)
    rbindlist(lapply(0:H, function(h) {
      fml <- paste0("y_h", h, " ~ dlog_wholesale", lag_controls, " | id^fuel")
      m <- feols(as.formula(fml), data = d, cluster = ~id)
      ct <- coeftable(m)
      data.table(
        horizon = h, group = g,
        estimate = ct["dlog_wholesale", "Estimate"],
        se = ct["dlog_wholesale", "Std. Error"]
      )
    }))
  }))
}

irf_pooled <- run_lp(grid, n_lags = n_lags)[, dimension := "pooled"]
irf_fuel <- run_lp(grid, "fuel", n_lags = n_lags)[, dimension := "combustible"]
irf_brand <- run_lp(grid, "brand_group", n_lags = n_lags)[, dimension := "marca"]
irf_comp <- run_lp(grid, "comp_tercile", n_lags = n_lags)[, dimension := "competidores_2km"]
irf_carretera <- run_lp(grid, "es_carretera", n_lags = n_lags)[, dimension := "carretera"]
irf_macrozona <- run_lp(grid, "macrozona", n_lags = n_lags)[, dimension := "macrozona"]

irf <- rbind(irf_pooled, irf_fuel, irf_brand, irf_comp, irf_carretera, irf_macrozona)
irf[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

fwrite(irf, file.path("results/tables", "local_projections_irf.csv"))

#####

# graphs

p_brand <- ggplot(irf[dimension == "marca"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "marca", fill = "marca",
    title = "local projections: pass-through de costo a precio, por marca",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "lp_irf_brand.pdf"), p_brand, width = 7.5, height = 4.5)

p_fuel <- ggplot(irf[dimension == "combustible"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "combustible", fill = "combustible",
    title = "local projections: pass-through de costo a precio, completamente pooled por combustible",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "lp_irf_fuel.pdf"), p_fuel, width = 7.5, height = 4.5)

p_comp <- ggplot(irf[dimension == "competidores_2km"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "competidores\na 2km (tercil)", fill = "competidores\na 2km (tercil)",
    title = "local projections: pass-through de costo a precio, por competencia",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "lp_irf_competitors.pdf"), p_comp, width = 7.5, height = 4.5)

p_carretera <- ggplot(irf[dimension == "carretera"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "tipo de estacion", fill = "tipo de estacion",
    title = "local projections: pass-through de costo a precio, carretera vs no",
    subtitle = "carretera = tiene duchas o surtidor de camiones; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "lp_irf_carretera.pdf"), p_carretera, width = 7.5, height = 4.5)

p_macrozona <- ggplot(irf[dimension == "macrozona"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "macrozona", fill = "macrozona",
    title = "local projections: pass-through de costo a precio, por macrozona",
    subtitle = "norte: arica-coquimbo; centro: valparaiso-biobio; sur: araucania-magallanes",
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "lp_irf_macrozona.pdf"), p_macrozona, width = 7.5, height = 4.5)
