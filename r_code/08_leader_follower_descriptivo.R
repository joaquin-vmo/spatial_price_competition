
# objective: descriptive look at *when* (date + time-of-day) each
# distributor's stations adopted the march 2026 mepco shock in their retail
# price, as a first pass before any formal leader/follower model.
#
# unit of analysis is the distributor, not the individual station: see
# 07/08 data-quality check (2026-07-13) -- a large share of timestamps are
# simultaneous batch uploads (e.g. 189 shell stations at the identical
# second 10:06:33-34 on 2026-03-26; 1,551 shell stations within a 20-second
# window at ~01:06), so intra-day timing looks like "when the distributor's
# system pushed the update," not a per-station decision. station-level
# leader/follower would just be measuring batch-upload infrastructure.

library(data.table)
library(fixest)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))

event_date <- as.Date("2026-03-26")
fuels <- c("93", "97", "di")

# each station-fuel's *first* price change on or after the shock (not every
# change in the window -- a station may reprice again later for unrelated
# reasons, which would just dilute "time to adopt this shock")
post <- db[fuel %in% fuels & date >= event_date & date <= event_date + 10 &
             price_changed == TRUE & time != ""]

post[, timestamp := as.POSIXct(paste(date, time), tz = "America/Santiago")]
setorder(post, id, fuel, timestamp)
window <- post[, .SD[1], by = .(id, fuel)]

shock_timestamp <- as.POSIXct(paste(event_date, "00:00:00"), tz = "America/Santiago")
window[, hours_since_shock := as.numeric(difftime(timestamp, shock_timestamp, units = "hours"))]

#####

# keep distributors with enough change events in the window to say
# anything; batch-flag distributors whose typical simultaneous-timestamp
# cluster size is large (see 04's distributor_n_stations / the check above)
dist_summary <- window[, .(
  n_changes = .N,
  n_stations = uniqueN(id),
  first_change = min(timestamp),
  median_hours = median(hours_since_shock),
  p25_hours = quantile(hours_since_shock, 0.25),
  p75_hours = quantile(hours_since_shock, 0.75)
), by = distributor]

dist_summary <- dist_summary[n_stations >= 15]
setorder(dist_summary, median_hours)

fwrite(dist_summary, file.path("results/tables", "leader_follower_distributor_summary.csv"))

print(dist_summary)

#####

plot_data <- window[distributor %in% dist_summary$distributor]
plot_data[, distributor := factor(distributor, levels = dist_summary$distributor)]

p <- ggplot(plot_data, aes(x = distributor, y = hours_since_shock)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    x = NULL,
    y = "horas desde el shock (2026-03-26 00:00) hasta el cambio de precio",
    title = "cuando cada distribuidor adopto el shock mepco de marzo 2026",
    subtitle = "ordenado por mediana; caution: varios distribuidores suben precios en lotes masivos\nvia sistema central, no estacion por estacion (ver comentario al inicio del script)"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "leader_follower_boxplot.pdf"), p, width = 8, height = 6)

#####

# grande-vs-chico gap, using the *full* sample of stations (not just the 5
# named distributors above -- most of "independiente" is individually too
# small to name, but pooled together they're a real comparison group).
# is_franchise / distributor_n_stations come straight from
# 04_build_database.R (sin bandera and single-station operators = chico).

pre_obs <- db[fuel %in% fuels & date < event_date][order(id, fuel, -date)][, .SD[1], by = .(id, fuel)]
pre_vars <- pre_obs[, .(id, fuel, is_franchise_pre = is_franchise,
                         distributor_n_stations_pre = distributor_n_stations,
                         competitors_2km_pre = competitors_2km,
                         competition_intensity_pre = competition_intensity)]

window <- merge(window, pre_vars, by = c("id", "fuel"), all.x = TRUE)
window[, grupo := fifelse(is_franchise_pre, "grande (cadena)", "chico (independiente)")]

cat("\n=== n estaciones-combustible con respuesta valida, por grupo ===\n")
print(window[, .N, by = grupo])

# adoption curve: % of stations that have already repriced by hour h
p_ecdf <- ggplot(window, aes(hours_since_shock, color = grupo)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 200)) +
  labs(
    x = "horas desde el shock (2026-03-26 00:00)",
    y = "% de estaciones que ya ajustaron su precio",
    color = NULL,
    title = "velocidad de adopcion del shock mepco de marzo 2026: grandes vs chicos",
    subtitle = "grande = cadena (is_franchise); chico = sin bandera u operador de 1 sola estacion"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "leader_follower_adoption_curve.pdf"), p_ecdf, width = 7.5, height = 4.5)

# quantify the gap: log(hours+1) so same-day changes (hours near 0) don't
# blow up the scale; clustered by distributor because a distributor's
# stations don't reprice independently of each other (batch uploads)
m_binary <- feols(log(hours_since_shock + 1) ~ is_franchise_pre, data = window, cluster = ~distributor)

# continuous version: does adoption speed scale with how big the chain is,
# not just big-3-or-not? also checks whether local competition independently
# predicts speed (a race-to-reprice story, distinct from chain size)
m_continuous <- feols(
  log(hours_since_shock + 1) ~ log(distributor_n_stations_pre) + competitors_2km_pre,
  data = window, cluster = ~distributor
)

etable(
  m_binary, m_continuous,
  tex = TRUE,
  file = file.path("results/tables", "leader_follower_regressions.tex"),
  replace = TRUE
)

print(etable(m_binary, m_continuous))

#####

# does competition speed up (or slow down) cost pass-through? uses
# competition_intensity (the inverse-distance-weighted continuous measure
# from 04_build_database.R, not the discrete competitors_2km count) as the
# predictor of speed = hours_since_shock for that station's first repricing.
# model 1 is unconditional; model 2 nets out chain size, since 08 already
# showed chain size is a strong, independent driver of speed and it could
# also be correlated with local competition (e.g. majors clustering in
# denser urban markets), which would confound a naive reading of model 1.

m_speed_uncontrolled <- feols(
  log(hours_since_shock + 1) ~ competition_intensity_pre,
  data = window, cluster = ~distributor
)

m_speed_controlled <- feols(
  log(hours_since_shock + 1) ~ competition_intensity_pre + log(distributor_n_stations_pre),
  data = window, cluster = ~distributor
)

etable(
  m_speed_uncontrolled, m_speed_controlled,
  tex = TRUE,
  file = file.path("results/tables", "competition_speed_regressions.tex"),
  replace = TRUE
)

print(etable(m_speed_uncontrolled, m_speed_controlled))

p_speed <- ggplot(window, aes(competition_intensity_pre, hours_since_shock)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, color = "steelblue") +
  scale_y_continuous(limits = c(0, 200)) +
  labs(
    x = "competition_intensity (pre-shock)",
    y = "horas desde el shock hasta el primer cambio de precio",
    title = "competencia vs. velocidad de traspaso del shock mepco de marzo 2026",
    subtitle = "cada punto = una estacion-combustible; recta = ajuste lineal simple (ver regresion clusterizada en la tabla)"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "competition_speed_scatter.pdf"), p_speed, width = 7, height = 4.5)
