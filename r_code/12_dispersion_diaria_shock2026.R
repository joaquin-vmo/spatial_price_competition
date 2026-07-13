
# objective: same question as 11_dispersion_shock2026.R (does price
# dispersion spike right when the march 2026 mepco shock hits, then close
# as laggards catch up?) but at daily resolution instead of weekly. the
# weekly version came back flat -- no visible spike -- which is plausible
# given 08_leader_follower_descriptivo.R found the chains-vs-independents
# adjustment gap resolves within ~40 hours (under 2 days): a 7-day bin
# would average right over a transient that short. resampling to a daily
# calendar (locf, same technique as the mepco-week grids in
# 06/07/10) should be fine-grained enough to catch it if it's there.

library(data.table)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))

event_date <- as.Date("2026-03-26")
window_days <- 20
fuels <- c("93", "97", "di")

dbf <- db[fuel %in% fuels]

#####

# resample each station-fuel's price series onto a daily calendar around
# the shock (locf), analogous to the mepco-week grids elsewhere, but here
# the grid is plain calendar days since we want daily, not mepco-weekly,
# resolution
day_grid <- seq(event_date - window_days, event_date + window_days, by = "day")

id_fuel_bounds <- dbf[, .(min_date = min(date), max_date = max(date)), by = .(id, fuel)]
id_fuel_bounds <- id_fuel_bounds[min_date <= min(day_grid) & max_date >= max(day_grid)]

grid <- id_fuel_bounds[, .(date = day_grid), by = .(id, fuel)]

price_obs <- dbf[, .(id, fuel, date, price, municipality)]
price_obs <- price_obs[, .(price = price[.N], municipality = municipality[.N]), by = .(id, fuel, date)]
setkey(price_obs, id, fuel, date)
setkey(grid, id, fuel, date)
grid <- price_obs[grid, roll = TRUE]

grid[, days_to_event := as.integer(date - event_date)]
grid[, log_price := log(price)]

#####

# dispersion within each comuna-fuel-day (only where >=3 stations, same
# minimum as the weekly version)
comuna_day <- grid[, .(sd_log_price = sd(log_price), n_stations = .N),
                    by = .(municipality, fuel, days_to_event)]
comuna_day <- comuna_day[n_stations >= 3]

cat("=== comuna-fuel-dia con dispersion calculable ===\n")
cat(nrow(comuna_day), "celdas\n")

dispersion <- comuna_day[, .(
  median_sd = median(sd_log_price),
  ci_low = quantile(sd_log_price, 0.25),
  ci_high = quantile(sd_log_price, 0.75),
  n_comunas = .N
), by = days_to_event]

setorder(dispersion, days_to_event)
fwrite(dispersion, file.path("results/tables", "dispersion_diaria_shock2026.csv"))
print(dispersion)

#####

p <- ggplot(dispersion, aes(days_to_event, median_sd)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15) +
  geom_line() +
  geom_point(size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "dias relativos al shock (2026-03-26)",
    y = "dispersion: sd(log precio) dentro de comuna, mediana entre comunas",
    title = "dispersion de precios alrededor del shock mepco de marzo 2026 (resolucion diaria)",
    subtitle = "comunas con >=3 estaciones activas ese dia; 93/97/di agrupados; banda = P25-P75"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "dispersion_diaria_shock2026.pdf"), p, width = 8, height = 4.5)
