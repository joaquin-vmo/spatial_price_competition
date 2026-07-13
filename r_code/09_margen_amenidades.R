
# objective: retail margin (price minus the vigente wholesale price with
# mepco) around the march 2026 mepco shock, by amenity level -- simple
# stations vs. those with a convenience store vs. full highway rest-stops
# (store + showers/truck pump). same median-based approach as
# 05_passthrough_analysis.R (see the outlier note there: raw margin has a
# long right tail from data-loading errors, so mean would be misleading).

library(data.table)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))

event_date <- as.Date("2026-03-26")
window_weeks <- 14 # data currently run through 2026-07-12: ~15 post weeks max

# mepco only stabilizes 93/97/di; margin needs a wholesale reference, so 95
# (no mepco quote) is dropped
db <- db[fuel %in% c("93", "97", "di") & !is.na(wholesale_w)]

# keep stations alive through the whole window (true entry/exit dates)
alive <- unique(db[, .(id, first_date, last_date)])
keep_ids <- alive[first_date <= event_date - window_weeks * 7 &
                     last_date >= event_date + window_weeks * 7, id]
db <- db[id %in% keep_ids]

db <- db[abs(days_to_event) <= window_weeks * 7]
db[, event_week := floor(days_to_event / 7)]

# amenity tier, fixed per station (physical attribute, not time-varying):
# simple = no convenience store; con_tienda = store but no showers/truck
# pump; completa = store + at least one of showers/truck pump (the
# highway-rest-stop proxy from 06_local_projections.R)
amenidad <- unique(db[, .(id, tienda_conveniencia, duchas, surtidor_camiones)])
amenidad <- amenidad[!is.na(tienda_conveniencia)]
amenidad[, nivel_amenidades := fifelse(
  !tienda_conveniencia, "simple",
  fifelse(tienda_conveniencia & !duchas & !surtidor_camiones, "con_tienda", "completa")
)]
db <- merge(db, amenidad[, .(id, nivel_amenidades)], by = "id", all.x = TRUE)
db <- db[!is.na(nivel_amenidades)]
db[, nivel_amenidades := factor(nivel_amenidades, levels = c("simple", "con_tienda", "completa"))]

db[, margin := price - wholesale_w]

#####

# two-step aggregation: median within station-week first (so a station
# reporting more often doesn't dominate), then median across stations
# within each amenity tier. iqr band instead of a normal-approximation ci

station_week <- db[, .(margin = median(margin)), by = .(id, event_week, nivel_amenidades)]

group_week <- station_week[, .(
  median_margin = median(margin),
  ci_low = quantile(margin, 0.25),
  ci_high = quantile(margin, 0.75),
  n_stations = .N
), by = .(event_week, nivel_amenidades)]

fwrite(group_week, file.path("results/tables", "margin_by_amenidades_week.csv"))

cat("=== n estaciones por nivel de amenidades ===\n")
print(unique(db[, .(id, nivel_amenidades)])[, .N, by = nivel_amenidades])

#####

p_margin <- ggplot(group_week, aes(event_week, median_margin, color = nivel_amenidades, fill = nivel_amenidades)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "margen: precio - costo mayorista con mepco ($/L)",
    color = "nivel de amenidades",
    fill = "nivel de amenidades",
    title = "margen retail alrededor del shock mepco de marzo 2026, por amenidades",
    subtitle = "simple: sin tienda; con_tienda: tienda pero sin duchas/surtidor camiones; completa: tienda + duchas o surtidor camiones\nmediana entre 93/97/di (banda = P25-P75)"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "margin_by_amenidades_week.pdf"), p_margin, width = 7.5, height = 5)
