
# objective: local projections for the single march 2026 mepco shock only
# (not pooled over history like 06_local_projections.R). the point of doing
# this with dlog_wholesale as the regressor, instead of a plain pre/post
# dummy like 05_event_study.R, is that the shock wasn't the same size for
# every fuel that week:
#   93: 1083.3 -> 1455.5  (dlog = 0.296)
#   97: 1133.1 -> 1524.6  (dlog = 0.297)
#   di:  831.2 -> 1411.5  (dlog = 0.530, almost double the gasolines)
# comparing raw price-level responses across fuels (as 05 does) conflates
# "responded more" with "got a bigger cost shock"; using dlog_wholesale as
# the dose fixes that.
#
# caveat: with a single shock date, there's no time-series variation left,
# so this is a cross-sectional regression (one row per station-fuel,
# clustered se by station), not a panel with station fe -- identification
# comes from the 3 dose levels across fuels plus cross-sectional
# heterogeneity (brand, competitors) within each.

library(data.table)
library(fixest)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))
mepco <- fread(file.path("data/processed", "whosale_prices_mepco.csv"))

H <- 12
fuels <- c("93", "97", "di")
event_date <- as.Date("2026-03-26")
prior_date <- as.Date("2026-03-19") # the mepco date right before the shock

dbf <- db[fuel %in% fuels & !is.na(wholesale_w)]

# keep stations alive through the whole window (true entry/exit dates, same
# logic as 05_event_study.R): avoids composition effects contaminating the
# horizon-12 estimate
alive <- unique(dbf[, .(id, first_date, last_date)])
keep_ids <- alive[first_date <= prior_date & last_date >= event_date + H * 7, id]
dbf <- dbf[id %in% keep_ids]

#####

# dose: dlog(wholesale) for each fuel between prior_date and event_date
dose <- mepco[date %in% c(prior_date, event_date), .SD, .SDcols = c("date", paste0(fuels, "_w/"))]
setnames(dose, paste0(fuels, "_w/"), fuels)
dose <- melt(dose, id.vars = "date", variable.name = "fuel", value.name = "wholesale_w")
dose <- dcast(dose, fuel ~ date, value.var = "wholesale_w")
setnames(dose, c("fuel", "wholesale_prior", "wholesale_event"))
dose[, dlog_wholesale := log(wholesale_event) - log(wholesale_prior)]
print(dose)

#####

# price just before the shock, and at h = 0..H weeks after, per station-fuel
# (a few id-fuel-date combos have more than one report that day; keep the
# last one so the rolling join below has a unique key)
price_obs <- dbf[, .(id, fuel, date, price)]
price_obs <- price_obs[, .(price = price[.N]), by = .(id, fuel, date)]
setkey(price_obs, id, fuel, date)

get_price_at <- function(target_date) {
  q <- unique(price_obs[, .(id, fuel)])
  q[, date := target_date]
  setkey(q, id, fuel, date)
  price_obs[q, roll = TRUE][, .(id, fuel, price)]
}

price_prior <- get_price_at(prior_date)
setnames(price_prior, "price", "price_prior")

lp_data <- price_prior
for (h in 0:H) {
  ph <- get_price_at(event_date + h * 7)
  setnames(ph, "price", paste0("price_h", h))
  lp_data <- merge(lp_data, ph, by = c("id", "fuel"))
}

lp_data <- merge(lp_data, dose[, .(fuel, dlog_wholesale)], by = "fuel")

for (h in 0:H) {
  lp_data[, paste0("y_h", h) := log(get(paste0("price_h", h))) - log(price_prior)]
}

#####

# heterogeneity groups, fixed at the last pre-shock observation (same
# construction as 05_event_study.R / 06_local_projections.R)
pre_obs <- dbf[date < event_date][order(id, fuel, -date)][, .SD[1], by = .(id, fuel)]

comp_breaks <- quantile(pre_obs$competitors_2km, probs = 0:3 / 3, na.rm = TRUE)
pre_obs[, comp_tercile := cut(competitors_2km, breaks = unique(comp_breaks),
                               include.lowest = TRUE, labels = c("baja", "media", "alta"))]
pre_obs[, brand_group := fifelse(distributor %in% c("copec", "shell", "aramco_petrobras"),
                                  distributor, "independiente")]

lp_data <- merge(lp_data, pre_obs[, .(id, fuel, comp_tercile, brand_group)], by = c("id", "fuel"))

#####

run_lp_cross_section <- function(data, group_col = NULL) {
  if (!is.null(group_col)) data <- data[!is.na(get(group_col))]
  groups <- if (is.null(group_col)) list(pooled = data) else split(data, data[[group_col]])
  rbindlist(lapply(names(groups), function(g) {
    d <- groups[[g]]
    rbindlist(lapply(0:H, function(h) {
      m <- feols(as.formula(paste0("y_h", h, " ~ dlog_wholesale")), data = d, cluster = ~id)
      ct <- coeftable(m)
      data.table(
        horizon = h, group = g,
        estimate = ct["dlog_wholesale", "Estimate"],
        se = ct["dlog_wholesale", "Std. Error"],
        n = nrow(d)
      )
    }))
  }))
}

irf_pooled <- run_lp_cross_section(lp_data)[, dimension := "pooled"]
irf_brand <- run_lp_cross_section(lp_data, "brand_group")[, dimension := "marca"]
irf_comp <- run_lp_cross_section(lp_data, "comp_tercile")[, dimension := "competidores_2km"]

irf <- rbind(irf_pooled, irf_brand, irf_comp)
irf[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

fwrite(irf, file.path("results/tables", "local_projections_evento2026_irf.csv"))

#####

p_brand <- ggplot(irf[dimension == "marca"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas desde el shock (2026-03-26)",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "marca", fill = "marca",
    title = "local projections del shock mepco de marzo 2026, por marca",
    subtitle = "dosis = dlog(mayorista) especifico de cada combustible ese dia; solo esta fecha, no la historia"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "lp_evento2026_brand.pdf"), p_brand, width = 7.5, height = 4.5)

p_comp <- ggplot(irf[dimension == "competidores_2km"], aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas desde el shock (2026-03-26)",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "competidores\na 2km (tercil)", fill = "competidores\na 2km (tercil)",
    title = "local projections del shock mepco de marzo 2026, por competencia",
    subtitle = "dosis = dlog(mayorista) especifico de cada combustible ese dia; solo esta fecha, no la historia"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "lp_evento2026_competitors.pdf"), p_comp, width = 7.5, height = 4.5)
