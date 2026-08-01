
# objective: the two identification checks that defend the maintained
# assumption -- that the TIMING of a rival's entry is exogenous to the
# incumbent, conditional on fixed effects.
#
#   A. event study with pre-trends (fischer, martin & schmidt-dengler 2025).
#      if entry were timed to local conditions that also move prices, treated
#      and control stations would diverge BEFORE the entry. flat pre-trends do
#      not prove exogeneity, but a joint wald test on the pre-period
#      coefficients is the sharpest available refutation, and it is reported
#      here rather than left to the eye.
#
#   B. fe ladder (arcidiacono, ellickson, mela & singleton 2020, table 6).
#      the selection of WHERE entrants locate is real and large. running the
#      same regression with progressively richer fixed effects makes it
#      visible: the cross-sectional specification loads on "entrants go to
#      cheap, dense markets", and only the station fixed effect isolates the
#      within-station change that the design actually uses.
#
# TAKES (data/processed, from 10_build_event_panel.R):
#   stack_entrada.csv / stack_entrada_far.csv, panel_mensual.csv
#
# PRODUCES:
#   results/tables/es_entrada_coefs.csv      event-study paths
#   results/tables/es_entrada_att.csv        static att per design/outcome
#   results/tables/es_entrada_pretrend.csv   joint wald test on pre-period
#   results/tables/fe_ladder_att.csv         the fe ladder
#   results/figures/es_entrada_precio.pdf, es_entrada_margen.pdf,
#     es_entrada_twfe.pdf, fe_ladder.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

PRE <- 12L
POST <- 12L
FOCAL_FROM <- 2014L

precios  <- c(p93 = "gasolina 93", p95 = "gasolina 95",
              p97 = "gasolina 97", pdi = "diesel")
margenes <- c(m93 = "gasolina 93", m97 = "gasolina 97", mdi = "diesel")

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------

# fixest coeftable -> data.table
tidy_coefs <- function(m) {
  ct <- coeftable(m)
  data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2], p = ct[, 4])
}

# event-study path from an i(<pattern>, ...) term, with tau = -1 added back at
# zero. scale = 100 turns a log-price coefficient into a percentage
tidy_es <- function(m, patron, scale = 100) {
  re <- paste0(patron, "::(-?\\d+)")
  d <- tidy_coefs(m)[grepl(re, term)]
  d[, event_time := as.integer(sub(paste0(".*", re, ".*"), "\\1", term))]
  d <- d[, .(event_time, estimate = estimate * scale, se = se * scale)]
  d <- rbind(d, data.table(event_time = -1L, estimate = 0, se = 0))
  d[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]
  setorder(d, event_time)[]
}

# joint test that every pre-period coefficient is zero. with ref = -1 omitted,
# all remaining negative event times ARE the pre-period, so the regex is enough
test_pretend <- function(m, patron) {
  w <- wald(m, paste0(patron, "::-"), print = FALSE)
  data.table(stat = w$stat, df = w$df1, p = w$p)
}

# static att: first coefficient of the formula, se clustered at comuna level
att_estatico <- function(fml, data) {
  ct <- coeftable(feols(fml, data = data, cluster = ~comuna))
  data.table(est = ct[1, 1], se = ct[1, 2], p = ct[1, 4])
}

fct_fuel <- function(x, lbl) factor(x, levels = names(lbl), labels = lbl)

plot_es <- function(d, ylab, subtitulo, color_var = NULL) {
  p <- if (is.null(color_var)) {
    ggplot(d, aes(event_time, estimate)) +
      geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15,
                  fill = AZUL, color = NA) +
      geom_line(color = AZUL) +
      geom_point(color = AZUL, size = 0.8)
  } else {
    ggplot(d, aes(event_time, estimate,
                  color = .data[[color_var]], fill = .data[[color_var]])) +
      geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.12, color = NA) +
      geom_line() +
      geom_point(size = 0.8)
  }
  p +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey55") +
    facet_wrap(~combustible, scales = "free_y") +
    labs(x = "meses desde la entrada del competidor",
         y = ylab, color = NULL, fill = NULL,
         title = "event study de la entrada de un competidor",
         subtitle = subtitulo) +
    tema()
}

# ==============================================================================
# A. EVENT STUDY WITH PRE-TRENDS
# ==============================================================================

stacks <- lapply(list(anillo = "stack_entrada", lejano = "stack_entrada_far"),
                 function(f) {
                   d <- fread(file.path("data/processed", paste0(f, ".csv")))
                   d[, year := year(as.IDate(ym))]
                   d[]
                 })
ctrl_lbl <- c(anillo = "control anillo (2-5 km)",
              lejano = "control lejano (>5 km, misma region)")

# stacked did: fe are event x station and event x month, so identification is
# strictly within an experiment. se clustered at comuna level
# brand x year absorbs chain-level pricing shocks. event_id^ym is a time effect
# common to the whole experiment, so a national move by one franchise in a given
# month survives it; in a market where three chains hold ~80% of stations that
# is exactly the confounder to worry about. distribuidor is a station-level
# constant by construction of the build, so this is not a bad control
FE_STACK <- "event_id^station_key + event_id^ym + distribuidor^year"

fml_es <- function(y, en_log) as.formula(sprintf(
  "%s ~ i(event_time, treated, ref = -1) | %s",
  if (en_log) sprintf("log(%s)", y) else y, FE_STACK
))
fml_att <- function(y, en_log) as.formula(sprintf(
  "%s ~ treated:post | %s",
  if (en_log) sprintf("log(%s)", y) else y, FE_STACK
))

run_stacked <- function(y, en_log, tipo) {
  rbindlist(lapply(names(stacks), function(cl) {
    d <- stacks[[cl]][!is.na(get(y))]
    m <- feols(fml_es(y, en_log), data = d, cluster = ~comuna)
    es <- tidy_es(m, "event_time", scale = if (en_log) 100 else 1)
    a  <- att_estatico(fml_att(y, en_log), d)
    pt <- test_pretend(m, "event_time")
    es[, `:=`(outcome = y, tipo = tipo, control = cl,
              att = a$est * if (en_log) 100 else 1,
              att_se = a$se * if (en_log) 100 else 1, att_p = a$p,
              pre_stat = pt$stat, pre_df = pt$df, pre_p = pt$p,
              n_obs = nobs(m), n_eventos = uniqueN(d$event_id))]
    es[]
  }))
}

res <- rbind(
  rbindlist(lapply(names(precios),  run_stacked, en_log = TRUE,  tipo = "precio")),
  rbindlist(lapply(names(margenes), run_stacked, en_log = FALSE, tipo = "margen"))
)

# ------------------------------------------------------------------------------
# A.2 twfe on the full monthly panel (fischer et al. baseline)
#
# treated = incumbent with an entry within 2 km; observations from a second
# nearby entry onward are dropped, so the estimate is the effect of a FIRST
# entry. the control pool excludes the entrants themselves: their own post-
# opening path is not a valid counterfactual
# ------------------------------------------------------------------------------

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g2_entry = as.IDate(g2_entry))]
panel[, year := year(ym)]

mi <- function(d) year(d) * 12L + (month(d) - 1L)

d_tw <- panel[is.na(g2_entry) | ym < g2_entry]
d_tw <- d_tw[!(role_entry == "treated" & year(g_entry) < FOCAL_FROM)]
d_tw <- d_tw[role_entry == "treated" | base_2012 == TRUE]
d_tw[, treated := as.integer(role_entry == "treated")]
d_tw[, et := miym - mi(g_entry)]
# endpoints binned: everything beyond the window is folded into +-12
d_tw[, rel := fifelse(treated == 1L, pmax(-PRE, pmin(POST, et)), -1L)]
d_tw[, post := fifelse(treated == 1L & !is.na(g_entry) & ym >= g_entry, 1L, 0L)]

FE_TW <- "station_key + ym + region^year + distribuidor^year"

res_tw <- rbindlist(lapply(names(precios), function(y) {
  d <- d_tw[!is.na(get(y))]
  m <- feols(as.formula(sprintf("log(%s) ~ i(rel, treated, ref = -1) | %s",
                                y, FE_TW)), data = d, cluster = ~comuna)
  a <- att_estatico(as.formula(sprintf("log(%s) ~ post | %s", y, FE_TW)), d)
  pt <- test_pretend(m, "rel")
  es <- tidy_es(m, "rel")
  es[, `:=`(outcome = y, tipo = "precio", control = "twfe panel",
            att = a$est * 100, att_se = a$se * 100, att_p = a$p,
            pre_stat = pt$stat, pre_df = pt$df, pre_p = pt$p,
            n_obs = nobs(m),
            # the panel design has no stacked experiments: what plays the role
            # of an "event" here is a treated station with its own cohort
            n_eventos = uniqueN(d[treated == 1L, station_key]))]
  es[]
}))

res <- rbind(res, res_tw)

# ------------------------------------------------------------------------------
# A.3 export
# ------------------------------------------------------------------------------

fwrite(res[, .(tipo, outcome, control, event_time, estimate, se,
               ci_low, ci_high)],
       file.path(TAB, "es_entrada_coefs.csv"))

att <- unique(res[, .(tipo, outcome, control, att, att_se, att_p,
                      n_obs, n_eventos)])
setorder(att, tipo, outcome, control)
fwrite(att, file.path(TAB, "es_entrada_att.csv"))

pretend <- unique(res[, .(tipo, outcome, control, pre_stat, pre_df, pre_p)])
setorder(pretend, tipo, outcome, control)
fwrite(pretend, file.path(TAB, "es_entrada_pretrend.csv"))

cat("\n=== ATT estatico de la entrada ===\n")
print(att[, .(tipo, outcome, control, att = round(att, 3),
              se = round(att_se, 3), p = round(att_p, 4), n_eventos)])

cat("\n=== TEST CONJUNTO DE PRE-TENDENCIAS (H0: todos los coef. tau<-1 = 0) ===\n")
cat("un p alto es lo que sostiene el supuesto: no hay divergencia previa\n")
print(pretend[, .(tipo, outcome, control, F = round(pre_stat, 2), df = pre_df,
                  p = round(pre_p, 4))])

# ------------------------------------------------------------------------------
# A.4 figures
# ------------------------------------------------------------------------------

d_fig <- res[tipo == "precio" & control %in% c("anillo", "lejano")]
d_fig[, `:=`(combustible = fct_fuel(outcome, precios),
             control = factor(ctrl_lbl[control], levels = ctrl_lbl))]
guardar(       plot_es(d_fig, "efecto sobre el precio (%)",
               "stacked did; fe evento x estacion y evento x mes; ee cluster comuna",
               "control"),
       file.path(FIG, "es_entrada_precio.pdf"), 9, 5.4)

d_mar <- res[tipo == "margen"]
d_mar[, `:=`(combustible = fct_fuel(outcome, margenes),
             control = factor(ctrl_lbl[control], levels = ctrl_lbl))]
guardar(       plot_es(d_mar, "efecto sobre el margen ($/L)",
               "margen = precio - costo mayorista mepco (nacional, semanal)",
               "control"),
       file.path(FIG, "es_entrada_margen.pdf"), 9, 4.2)

d_twf <- res[control == "twfe panel"]
d_twf[, combustible := fct_fuel(outcome, precios)]
guardar(       plot_es(d_twf, "efecto sobre el precio (%)",
               "twfe panel completo; fe estacion + mes + region x anio; extremos agrupados en +-12"),
       file.path(FIG, "es_entrada_twfe.pdf"), 9, 5.4)

# ==============================================================================
# B. FE LADDER (arcidiacono et al. 2020, table 6)
#
# same sample and same regressor throughout; only the fixed effects change.
# the point is not that one rung is right and the others wrong, but that the
# gap between them MEASURES the selection of entry locations. if the estimate
# collapses once station fixed effects go in, the cross-sectional association
# was picking up where entrants choose to be, not what their entry does.
# ==============================================================================

escalera <- list(
  `1. solo mes`                        = "ym",
  `2. mes + grupo tratado`             = "ym",           # + treated regressor
  `3. mes + comuna + distribuidor`     = "ym + comuna + distribuidor",
  `4. mes + estacion`                  = "ym + station_key",
  `5. mes + estacion + region x anio`  = "ym + station_key + region^year",
  `6. + marca x anio (preferida)`      = FE_TW,
  `7. estacion + comuna x mes + marca` = "station_key + comuna^ym + distribuidor^year"
)

fe_ladder <- rbindlist(lapply(names(precios), function(y) {
  d <- d_tw[!is.na(get(y))]
  rbindlist(lapply(names(escalera), function(nm) {
    rhs <- if (nm == "2. mes + grupo tratado") "post + treated" else "post"
    m <- feols(as.formula(sprintf("log(%s) ~ %s | %s", y, rhs, escalera[[nm]])),
               data = d, cluster = ~comuna)
    ct <- coeftable(m)
    data.table(outcome = y, spec = nm,
               est = ct["post", 1] * 100, se = ct["post", 2] * 100,
               p = ct["post", 4], n_obs = nobs(m))
  }))
}))
fe_ladder[, `:=`(ci_low = est - 1.96 * se, ci_high = est + 1.96 * se)]

fwrite(fe_ladder, file.path(TAB, "fe_ladder_att.csv"))

cat("\n=== ESCALERA DE EFECTOS FIJOS (efecto de la entrada, % del precio) ===\n")
cat("filas 1-3 no controlan por estacion: mezclan DONDE entran con QUE hace la entrada\n")
print(dcast(fe_ladder, spec ~ outcome, value.var = "est")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

fe_ladder[, `:=`(combustible = fct_fuel(outcome, precios),
                 spec = factor(spec, levels = rev(names(escalera))))]
p_lad <- ggplot(fe_ladder, aes(est, spec)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.15, color = AZUL) +
  geom_point(color = AZUL, size = 1.8) +
  facet_wrap(~combustible, scales = "free_x") +
  labs(x = "efecto de la entrada sobre el precio (%)", y = NULL,
       title = "escalera de efectos fijos: seleccion de ubicacion vs. efecto de la entrada",
       subtitle = "misma muestra y mismo regresor; solo cambian los efectos fijos") +
  tema()
guardar(p_lad, file.path(FIG, "fe_ladder.pdf"), 8.5, 4.5)

message("11_event_study_entrada.R: figuras en ", FIG, " y tablas en ", TAB, ".")
