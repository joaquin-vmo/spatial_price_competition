
# objective: descriptive analysis of the march 2026 mepco shock -- the
# retail margin (price minus the vigente wholesale price with mepco) for
# copec / shell / aramco-petrobras / independientes, pre and post shock

library(data.table)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))

event_date <- as.Date("2026-03-26")
window_weeks <- 14 # data currently run through 2026-07-12: ~15 post weeks max

# mepco only stabilizes 93/97/di; margin needs a wholesale reference, so 95
# (no mepco quote) is dropped
db <- db[fuel %in% c("93", "97", "di") & !is.na(wholesale_w)]

# keep stations alive through the whole window (true entry/exit dates, not
# just whether they happened to report exactly at the window's edge): avoids
# composition effects (entry/exit) contaminating the margin comparison
alive <- unique(db[, .(id, first_date, last_date)])
keep_ids <- alive[first_date <= event_date - window_weeks * 7 &
                    last_date >= event_date + window_weeks * 7, id]
db <- db[id %in% keep_ids]

db <- db[abs(days_to_event) <= window_weeks * 7]
db[, event_week := floor(days_to_event / 7)]

# brand fixed at the last pre-shock observation (constant per station in
# practice, but kept consistent with the other event-time scripts). the big
# 3 get their own group; every other distributor (small chains and "sin
# bandera" alike) is pooled into "independiente"
pre_obs <- db[date < event_date][order(id, fuel, -date)][, .SD[1], by = .(id, fuel)]
pre_obs[, brand_group_pre := fifelse(
  distributor %in% c("copec", "shell", "aramco_petrobras"),
  distributor,
  "independiente"
)]
pre_vars <- pre_obs[, .(id, fuel, brand_group_pre)]
db <- merge(db, pre_vars, by = c("id", "fuel"), all.x = TRUE)
db <- db[!is.na(brand_group_pre)]

db[, margin := price - wholesale_w]

# note on outliers (investigated 2026-07-13): ~665 rows nationwide have
# margin > 400, way above the typical ~100-180 range. these are not stale
# prices (days_since_last_price_change == 0 for all of them, i.e. they're
# freshly reported). two distinct sources: (1) a real regional effect --
# aysen/magallanes run genuinely higher margins across their *entire*
# sample (~150-180 avg vs ~80-95 in valparaiso/metropolitana), consistent
# with higher transport/logistics costs, not an error; and (2) isolated
# data-loading errors, e.g. the independiente spike right before the march
# 2026 shock was driven by 8 "custom service" stations (quilpue/villa
# alemana) all reporting the identical diesel price of 1530 on 2026-03-25,
# one day before mepco's wholesale price actually moved -- looks like a
# bulk reporting glitch from that one distributor, not real behavior. only
# ~21% of the 665 rows show a clean spike-then-revert error signature, and
# the rest are concentrated in 2025-2026 without a single clear cause.
# using the median (vs. mean) below is a simple, robust-enough fix for this
# descriptive plot without deciding on a global cleaning rule yet.

#####

# two-step aggregation: median within station-week first (so a station
# reporting more often doesn't dominate), then median across stations
# within each brand group. iqr-based band instead of a normal-approximation
# ci, since that pairs more naturally with a median than a mean/se does

station_week <- db[, .(margin = median(margin)), by = .(id, event_week, brand_group_pre)]

group_week <- station_week[, .(
  median_margin = median(margin),
  ci_low = quantile(margin, 0.25),
  ci_high = quantile(margin, 0.75),
  n_stations = .N
), by = .(event_week, brand_group_pre)]

fwrite(group_week, file.path("results/tables", "margin_by_brand_week.csv"))

#####

p_margin <- ggplot(group_week, aes(event_week, median_margin, color = brand_group_pre, fill = brand_group_pre)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "margen: precio - costo mayorista con mepco ($/L)",
    color = "marca",
    fill = "marca",
    title = "margen retail alrededor del shock mepco de marzo 2026",
    subtitle = "copec / shell / aramco-petrobras / independientes, mediana entre 93/97/di (banda = P25-P75)"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "margin_by_brand_week.pdf"), p_margin, width = 7.5, height = 4.5)

#####

# same shock, but the aggregate margin by fuel instead of by brand: pool all
# brands together and split only on 93 / 97 / di. same two-step median
# aggregation and iqr band as above, but keeping fuel (not collapsing over it)
#
# note: di at event_week -1 shows an inflated upper band (P75 ~ 675) driven by
# the 2026-03-25 reporting glitch documented in the outlier note above. the
# median (~162) is unaffected, so the point estimate is fine -- only that one
# cell's band is contaminated. nothing to fix here without committing to a
# global outlier-cleaning rule (deliberately deferred, see note above); left
# as-is on purpose

station_week_fuel <- db[, .(margin = median(margin)), by = .(id, event_week, fuel)]

group_week_fuel <- station_week_fuel[, .(
  median_margin = median(margin),
  ci_low = quantile(margin, 0.25),
  ci_high = quantile(margin, 0.75),
  n_stations = .N
), by = .(event_week, fuel)]

fwrite(group_week_fuel, file.path("results/tables", "margin_by_fuel_week.csv"))

p_margin_fuel <- ggplot(group_week_fuel, aes(event_week, median_margin, color = fuel, fill = fuel)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "margen: precio - costo mayorista con mepco ($/L)",
    color = "combustible",
    fill = "combustible",
    title = "margen retail agregado por combustible alrededor del shock mepco de marzo 2026",
    subtitle = "todas las marcas agrupadas, mediana entre estaciones (banda = P25-P75)"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "margin_by_fuel_week.pdf"), p_margin_fuel, width = 7.5, height = 4.5)

#####

# margin by competition: does more local competition compress the margin
# and/or slow its post-shock re-expansion? competitors_2km is time-varying and
# station-level (fuel-invariant), so -- like brand above -- fix it at the last
# pre-shock observation per station, then split into low/med/high terciles

comp_pre <- db[date < event_date][order(id, -date)][, .SD[1], by = id]
comp_breaks <- quantile(comp_pre$competitors_2km, probs = 0:3 / 3, na.rm = TRUE)
comp_pre[, comp_tercile_pre := cut(competitors_2km, breaks = unique(comp_breaks),
                                    include.lowest = TRUE, labels = c("baja", "media", "alta"))]
db <- merge(db, comp_pre[, .(id, comp_tercile_pre)], by = "id", all.x = TRUE)

# same two-step median aggregation and iqr band as the brand plot (grouping by
# a station-level key collapses over fuels within each station-week)
station_week_comp <- db[!is.na(comp_tercile_pre),
                        .(margin = median(margin)), by = .(id, event_week, comp_tercile_pre)]

group_week_comp <- station_week_comp[, .(
  median_margin = median(margin),
  ci_low = quantile(margin, 0.25),
  ci_high = quantile(margin, 0.75),
  n_stations = .N
), by = .(event_week, comp_tercile_pre)]

fwrite(group_week_comp, file.path("results/tables", "margin_by_competition_week.csv"))

p_margin_comp <- ggplot(group_week_comp, aes(event_week, median_margin, color = comp_tercile_pre, fill = comp_tercile_pre)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "margen: precio - costo mayorista con mepco ($/L)",
    color = "competidores\na 2km (tercil)",
    fill = "competidores\na 2km (tercil)",
    title = "margen retail por competencia local alrededor del shock mepco de marzo 2026",
    subtitle = "tercil de competidores a 2km fijado pre-shock, mediana entre 93/97/di (banda = P25-P75)"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "margin_by_competition_week.pdf"), p_margin_comp, width = 7.5, height = 4.5)

#####

# the "scissor": retail price and wholesale (mepco) cost on the same axis, by
# fuel, over event time. the margin is literally the vertical gap between the
# two lines -- the cost jumps at week 0 and the price chases it with a lag.
# wholesale_w is common per fuel-date, so its median across stations is just
# the vigente mepco cost; one facet per fuel since levels differ

scissor_sw <- db[, .(price = median(price), cost = median(wholesale_w)),
                 by = .(id, event_week, fuel)]
scissor <- scissor_sw[, .(price = median(price), cost = median(cost)),
                      by = .(event_week, fuel)]
fwrite(scissor, file.path("results/tables", "price_cost_by_fuel_week.csv"))

scissor_long <- melt(scissor, id.vars = c("event_week", "fuel"),
                     measure.vars = c("price", "cost"),
                     variable.name = "serie", value.name = "valor")
scissor_long[, serie := fifelse(serie == "price", "precio minorista", "costo mayorista (mepco)")]

p_scissor <- ggplot(scissor_long, aes(event_week, valor, color = serie)) +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  geom_line() +
  geom_point(size = 0.8) +
  facet_wrap(~fuel, scales = "free_y") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "$/L",
    color = NULL,
    title = "precio minorista vs costo mayorista mepco, por combustible, alrededor del shock de marzo 2026",
    subtitle = "el margen es la brecha vertical entre ambas lineas; mediana entre estaciones"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path("results/figures", "price_cost_by_fuel_week.pdf"), p_scissor, width = 9, height = 4.5)

#####

# margin normalized to the pre-shock baseline: re-center each brand group to
# its own mean margin over weeks -4..-1, so the plot shows the *change* from
# baseline (delta margin) rather than the level. this strips out the level
# gaps between groups (shell always runs higher, aysen/magallanes higher on
# transport cost) and isolates the speed and size of the adjustment -- the
# thing actually being compared across brands. band re-centered by the same
# baseline so it moves with the line

baseline_brand <- group_week[event_week %in% -4:-1, .(baseline = mean(median_margin)), by = brand_group_pre]
group_week_norm <- merge(group_week, baseline_brand, by = "brand_group_pre")
group_week_norm[, `:=`(
  d_median = median_margin - baseline,
  d_low = ci_low - baseline,
  d_high = ci_high - baseline
)]

fwrite(group_week_norm, file.path("results/tables", "margin_by_brand_week_normalized.csv"))

p_margin_norm <- ggplot(group_week_norm, aes(event_week, d_median, color = brand_group_pre, fill = brand_group_pre)) +
  geom_ribbon(aes(ymin = d_low, ymax = d_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "cambio en el margen vs baseline pre-shock ($/L)",
    color = "marca",
    fill = "marca",
    title = "ajuste del margen retail al shock mepco de marzo 2026, por marca (normalizado)",
    subtitle = "cada grupo re-centrado a su margen medio de las semanas -4 a -1; mediana entre 93/97/di (banda = P25-P75)"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "margin_by_brand_week_normalized.pdf"), p_margin_norm, width = 7.5, height = 4.5)
