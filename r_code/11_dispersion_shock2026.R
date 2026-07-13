
# objective: does price dispersion spike right when the march 2026 mepco
# shock hits, then close back down as laggards catch up to the fast
# movers? motivated by Solis Soto (2025, U. Chile) "Dispersion de precios
# en el mercado minorista de gasolina en Chile", which documents that
# price dispersion (rank reversals between nearby stations) persists in
# steady state -- but studies 2013-2022 with no single common shock. we
# have one, and already know from 08_leader_follower_descriptivo.R that
# chains reprice within ~1-2h while independents take ~10h-days, so a
# temporary widening of dispersion right at the shock (as some have
# repriced and others haven't yet) is a direct, testable prediction.
#
# unlike that paper's "etapa 1" (residualizing out cadena x comuna fixed
# effects), we use *raw* dispersion within comuna: residualizing on chain
# would remove exactly the variation this hypothesis is about (chains vs.
# independents resolving at different speeds).

library(data.table)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))

event_date <- as.Date("2026-03-26")
window_weeks <- 12
fuels <- c("93", "97", "di") # the shocked fuels

db <- db[fuel %in% fuels]
db <- db[abs(days_to_event) <= window_weeks * 7]
db[, event_week := floor(days_to_event / 7)]

# one price per station-comuna-fuel-week (last report that week), so a
# station reporting more than once doesn't get extra weight in the sd
db[, log_price := log(price)]
setorder(db, id, fuel, event_week, date)
station_week <- db[, .SD[.N], by = .(id, fuel, event_week), .SDcols = c("municipality", "log_price")]

#####

# dispersion within each comuna-fuel-week (only where >=3 stations report,
# the minimum for an sd that means anything)
comuna_week <- station_week[, .(sd_log_price = sd(log_price), n_stations = .N),
                             by = .(municipality, fuel, event_week)]
comuna_week <- comuna_week[n_stations >= 3]

cat("=== comuna-fuel-semana con dispersion calculable ===\n")
cat(nrow(comuna_week), "celdas, de", uniqueN(station_week$municipality), "comunas totales\n")

# median across comunas per event week, not mean: a handful of comunas hit
# by the same kind of data-loading error documented in
# 05_passthrough_analysis.R (e.g. diesel ~0.25-0.30 log-point comuna sd's
# concentrated right at event_week -1, 2026-03-25) would otherwise blow up
# a simple average. iqr band instead of a normal-approximation ci for the
# same reason.
dispersion <- comuna_week[, .(
  mean_sd = median(sd_log_price),
  ci_low = quantile(sd_log_price, 0.25),
  ci_high = quantile(sd_log_price, 0.75),
  n_comunas = .N
), by = event_week]

setorder(dispersion, event_week)
fwrite(dispersion, file.path("results/tables", "dispersion_shock2026.csv"))
print(dispersion)

#####

p <- ggplot(dispersion, aes(event_week, mean_sd)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "dispersion: sd(log precio) dentro de comuna, mediana entre comunas",
    title = "dispersion de precios alrededor del shock mepco de marzo 2026",
    subtitle = "comunas con >=3 estaciones reportando esa semana; 93/97/di agrupados; banda = P25-P75"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "dispersion_shock2026.pdf"), p, width = 7.5, height = 4.5)
