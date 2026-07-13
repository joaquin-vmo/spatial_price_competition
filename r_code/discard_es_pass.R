
# objective: event-study of retail price response to the march 2026 mepco
# wholesale shock (+372 on 93 w/, 2026-03-26)
#
# identification: the shock hits every station on the same calendar date, so
# there is no untreated control group -- station and (fuel x date) fixed
# effects absorb the common shock entirely. what's identified is the
# *differential* response: event-time dummies interacted with a station's
# pre-shock competition intensity / franchise status. see rockets-and-
# feathers literature (borenstein-cameron-gilbert 1997, lewis 2011) for the
# same logic applied to distributed-lag pass-through regressions.

library(data.table)
library(fixest)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))

event_date <- as.Date("2026-03-26")
window_weeks <- 14 # data currently run through 2026-07-12: ~15 post weeks max

# mepco only stabilizes 93/97/di; 95 has no wholesale reference series
db <- db[fuel %in% c("93", "97", "di")]

# keep stations alive through the whole window (using each station's true
# entry/exit dates, not just whether it happened to report exactly at the
# window's edge): avoids composition effects (entry/exit) contaminating the
# estimated price path
alive <- unique(db[, .(id, first_date, last_date)])
keep_ids <- alive[first_date <= event_date - window_weeks * 7 &
                     last_date >= event_date + window_weeks * 7, id]
db <- db[id %in% keep_ids]

db <- db[abs(days_to_event) <= window_weeks * 7]
db[, event_week := floor(days_to_event / 7)]

# heterogeneity groups fixed at the last pre-shock observation, so a
# station's own group can't drift across the event window
pre_obs <- db[date < event_date][order(id, fuel, -date)][, .SD[1], by = .(id, fuel)]

comp_breaks <- quantile(pre_obs$competition_intensity, probs = 0:3 / 3, na.rm = TRUE)
pre_obs[, comp_tercile := cut(competition_intensity, breaks = comp_breaks,
                               include.lowest = TRUE, labels = c("baja", "media", "alta"))]

pre_vars <- pre_obs[, .(id, fuel, comp_tercile, is_franchise_pre = is_franchise)]
db <- merge(db, pre_vars, by = c("id", "fuel"), all.x = TRUE)
db <- db[!is.na(comp_tercile)]

#####

# model a: raw average price path (station fe only, no date fe -- purely
# descriptive, since with date fe the common path isn't separately
# identified from the fixed effect; picks up general inflation/trend too)

m_raw <- feols(price ~ i(event_week, ref = -1) | id^fuel,
                data = db, cluster = ~id)

# model b: differential path by pre-shock competition tercile (vs "baja",
# the reference group) -- identified within date because a station's
# tercile varies across stations reporting on the same date

m_competition <- feols(
  price ~ i(event_week, comp_tercile, ref = -1, ref2 = "baja") | id^fuel + date^fuel,
  data = db, cluster = ~id
)

# model c: differential path, franchise vs independent

m_franchise <- feols(
  price ~ i(event_week, is_franchise_pre, ref = -1) | id^fuel + date^fuel,
  data = db, cluster = ~id
)

# model d: both interactions together

m_combined <- feols(
  price ~ i(event_week, comp_tercile, ref = -1, ref2 = "baja") +
    i(event_week, is_franchise_pre, ref = -1) | id^fuel + date^fuel,
  data = db, cluster = ~id
)

etable(
  m_raw, m_competition, m_franchise, m_combined,
  tex = TRUE,
  file = file.path("results/tables", "event_study.tex"),
  replace = TRUE
)

#####

# event-time plots

extract_path <- function(model, pattern, label) {
  ct <- as.data.table(coeftable(model), keep.rownames = "term")
  ct <- ct[grepl(pattern, term)]
  ct[, event_week := as.integer(gsub(".*event_week::(-?[0-9]+).*", "\\1", term))]
  ct[, `:=`(
    estimate = Estimate,
    ci_low = Estimate - 1.96 * `Std. Error`,
    ci_high = Estimate + 1.96 * `Std. Error`,
    series = label
  )]
  rbind(
    ct[, .(event_week, estimate, ci_low, ci_high, series)],
    data.table(event_week = -1L, estimate = 0, ci_low = 0, ci_high = 0, series = label)
  )
}

path_raw <- extract_path(m_raw, "event_week", "precio promedio")

p_raw <- ggplot(path_raw, aes(event_week, estimate)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "precio, relativo a la semana previa al shock ($/L)",
    title = "trayectoria de precios alrededor del shock mepco de marzo 2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "event_study_raw.pdf"), p_raw, width = 7, height = 4.5)

extract_tercile_path <- function(model) {
  ct <- as.data.table(coeftable(model), keep.rownames = "term")
  ct <- ct[grepl("^event_week::.*:comp_tercile::", term)]
  ct[, event_week := as.integer(gsub("event_week::(-?[0-9]+):comp_tercile::.*", "\\1", term))]
  ct[, tercile := gsub("event_week::-?[0-9]+:comp_tercile::(.*)", "\\1", term)]
  ct[, `:=`(
    estimate = Estimate,
    ci_low = Estimate - 1.96 * `Std. Error`,
    ci_high = Estimate + 1.96 * `Std. Error`
  )]
  rbind(
    ct[, .(event_week, tercile, estimate, ci_low, ci_high)],
    data.table(event_week = -1L, tercile = c("media", "alta"), estimate = 0, ci_low = 0, ci_high = 0)
  )
}

path_comp <- extract_tercile_path(m_competition)

p_comp <- ggplot(path_comp, aes(event_week, estimate, color = tercile, fill = tercile)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "diferencia vs. tercil de competencia baja ($/L)",
    color = "tercil de competencia\npre-shock (competition_intensity)",
    fill = "tercil de competencia\npre-shock (competition_intensity)",
    title = "pass-through diferencial por tercil de competencia pre-shock"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "event_study_competition.pdf"), p_comp, width = 7, height = 4.5)

path_franchise <- extract_path(m_franchise, "event_week", "franquicia vs independiente")

p_franchise <- ggplot(path_franchise, aes(event_week, estimate)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = -0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = "semanas relativas al shock (2026-03-26)",
    y = "diferencia franquicia - independiente ($/L)",
    title = "pass-through diferencial: franquicia vs independiente"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "event_study_franchise.pdf"), p_franchise, width = 7, height = 4.5)
