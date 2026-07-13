
# objective: "rockets and feathers" -- is cost pass-through faster/more
# complete when the mepco wholesale price rises than when it falls? the
# central question in this literature (bacon 1991; borenstein-cameron-
# gilbert 1997; bachmeier-griffin 2003; peltzman 2000). same local
# projections machinery as 06_local_projections.R (pooled over every mepco
# update 2014-2026, no date fe -- see that script's header note on why),
# just split the sample by the sign of dlog_wholesale instead of by brand
# or competition.

library(data.table)
library(fixest)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))
mepco <- fread(file.path("data/processed", "whosale_prices_mepco.csv"))

H <- 12
fuels <- c("93", "97", "di")

dbf <- db[fuel %in% fuels & !is.na(wholesale_w)]

#####

# resample onto the mepco weekly calendar (locf), same construction as 06
mepco_dates <- sort(unique(as.Date(mepco$date)))
id_fuel_bounds <- dbf[, .(min_date = min(date), max_date = max(date)), by = .(id, fuel)]
grid <- id_fuel_bounds[, .(date = mepco_dates[mepco_dates >= min_date & mepco_dates <= max_date]),
                        by = .(id, fuel)]

price_obs <- dbf[, .(id, fuel, date, price)]
setkey(price_obs, id, fuel, date)
setkey(grid, id, fuel, date)
grid <- price_obs[grid, roll = TRUE]

mepco_long <- rbindlist(lapply(fuels, function(f) {
  data.table(fuel = f, date = mepco$date, wholesale_w = mepco[[paste0(f, "_w/")]])
}))
setorder(mepco_long, fuel, date)
mepco_long[, log_w := log(wholesale_w)]
mepco_long[, dlog_wholesale := log_w - shift(log_w), by = fuel]

# cumulative wholesale change over [t-1, t+h], same window as y_h below.
# using only the origin week's dlog_wholesale as the regressor at every
# horizon is the wrong specification here: it's correlated ~0.40-0.44 with
# how much wholesale keeps moving over the following weeks (mepco doesn't
# mean-revert), so that single-week shock acts as a proxy for a larger
# cumulative move and inflates beta_h at longer horizons (checked: gave
# beta_h > 2 for "baja" at h=8, not economically sensible). regressing on
# the cumulative change over the same window as the price response is the
# standard fix (as in exchange-rate pass-through regressions).
for (h in 0:H) {
  mepco_long[, paste0("cumw_h", h) := shift(log_w, -h) - shift(log_w, 1), by = fuel]
}

grid <- merge(grid, mepco_long, by = c("fuel", "date"), all.x = TRUE)

setorder(grid, id, fuel, date)
grid[, log_price := log(price)]
for (h in 0:H) {
  grid[, paste0("y_h", h) := shift(log_price, -h) - shift(log_price, 1), by = .(id, fuel)]
}

# sign of that week's shock; drop exact-zero weeks (uninformative either way)
grid[, sentido := fifelse(dlog_wholesale > 0, "sube", fifelse(dlog_wholesale < 0, "baja", NA_character_))]

#####

run_lp <- function(data, group_col) {
  data <- data[!is.na(get(group_col))]
  groups <- split(data, data[[group_col]])
  rbindlist(lapply(names(groups), function(g) {
    d <- groups[[g]]
    rbindlist(lapply(0:H, function(h) {
      x_h <- paste0("cumw_h", h)
      m <- feols(as.formula(paste0("y_h", h, " ~ ", x_h, " | id^fuel")), data = d, cluster = ~id)
      ct <- coeftable(m)
      data.table(horizon = h, group = g,
                 estimate = ct[x_h, "Estimate"],
                 se = ct[x_h, "Std. Error"],
                 n = nrow(d))
    }))
  }))
}

irf <- run_lp(grid, "sentido")
irf[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

fwrite(irf, file.path("results/tables", "rockets_feathers_irf.csv"))
print(irf[horizon %in% c(0, 1, 4, 8, 12)])

#####

p <- ggplot(irf, aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c(sube = "firebrick", baja = "steelblue")) +
  scale_fill_manual(values = c(sube = "firebrick", baja = "steelblue")) +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "sentido del\ncambio mayorista", fill = "sentido del\ncambio mayorista",
    title = "rockets and feathers: pass-through segun el signo del shock de costo",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "rockets_feathers_irf.pdf"), p, width = 7.5, height = 4.5)
