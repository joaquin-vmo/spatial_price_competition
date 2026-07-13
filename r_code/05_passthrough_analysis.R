
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
