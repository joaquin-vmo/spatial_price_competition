
# objective: the raw data behind the event study, and a diagnosis of what the
# fixed effects are doing to it.
#
# this is NOT the parallel-trends test -- the leads of the event study are, and
# the joint wald plus the slope test in 11/13 are sharper than any figure.
#
# and the raw series should NOT be presented as parallel-trends evidence, for a
# reason this script documents: the crude within-event difference carries a
# significant pre-entry slope (about -0.11 $/L per month, p < 0.005) that the
# estimated event studies do not have. section 3 below shows why. the drift is
# killed by the experiment-specific time effect (event_id^ym), and the
# specification without it returns att estimates of inconsistent sign, so it is
# not a credible unconditional benchmark -- it is just contaminated by the
# national mepco cost path. showing the raw figure as reassurance would hand a
# reader a picture that appears to contradict the paper's own tests.
#
# why event time and not calendar time: treatment is staggered across 2014-2026,
# so there is no common "before" and "after" on a calendar axis -- at any given
# date the treated group mixes stations that have already received entry with
# stations that have not. a calendar plot would show the national mepco cost
# path and the level gap between groups (which the station fixed effects absorb
# anyway), and nothing about parallel trends.
#
# the sample is the stacked design, where every control is matched to a specific
# entry event and therefore has a well-defined event time. pairs are kept only
# if observed at EVERY tau in the window, so that a change across the x axis is
# a change in prices and not a change in which stations are being averaged.
#
# TAKES:  data/processed/stack_entrada_far.csv  (from 10_build_event_panel.R)
# PRODUCES:
#   results/figures/descriptivo_niveles.pdf     raw means by event time
#   results/figures/descriptivo_normalizado.pdf same, each series set to 0 at -1
#   results/tables/descriptivo_medias.csv

library(data.table)
library(ggplot2)

FIG <- "results/figures"
TAB <- "results/tables"

VENTANA <- 12L

outcomes <- c(p93 = "precio gasolina 93 ($/L)",
              pdi = "precio diesel ($/L)",
              m93 = "margen gasolina 93 ($/L)",
              mdi = "margen diesel ($/L)")

stk <- fread("data/processed/stack_entrada_far.csv")
stk <- stk[event_time >= -VENTANA & event_time <= VENTANA]

# balanced event-station pairs: without this the series would move because the
# composition of stations changes with event time, not because prices do
stk[, par := paste(event_id, station_key)]
completos <- stk[, .N, by = par][N == 2L * VENTANA + 1L, par]
bal <- stk[par %in% completos]

cat("=== MUESTRA BALANCEADA ===\n")
cat("pares evento-estacion: ", uniqueN(stk$par), " -> ",
    uniqueN(bal$par), " tras exigir los ", 2 * VENTANA + 1,
    " meses completos\n", sep = "")
print(bal[, .(pares = uniqueN(par), eventos = uniqueN(event_id)),
          by = .(grupo = fifelse(treated == 1L, "tratadas", "control"))])

largo <- melt(bal, id.vars = c("par", "event_time", "treated"),
              measure.vars = names(outcomes),
              variable.name = "outcome", value.name = "valor",
              variable.factor = FALSE)
largo <- largo[!is.na(valor)]

medias <- largo[, .(media = mean(valor), n = .N),
                by = .(outcome, event_time,
                       grupo = fifelse(treated == 1L, "tratadas (<2 km)",
                                       "control (>5 km, misma region)"))]
# each series recentred on its own tau = -1, which is the reference the event
# study normalises to; this is what makes the pre-entry comovement readable
medias[, base := media[event_time == -1L], by = .(outcome, grupo)]
medias[, norm := media - base]

fwrite(medias, file.path(TAB, "descriptivo_medias.csv"))

medias[, `:=`(
  serie = factor(outcome, levels = names(outcomes), labels = outcomes),
  grupo = factor(grupo, levels = c("tratadas (<2 km)",
                                   "control (>5 km, misma region)")))]

graficar <- function(yvar, ylab, subt) {
  ggplot(medias, aes(event_time, get(yvar), color = grupo)) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey55") +
    geom_line(linewidth = 0.55) +
    geom_point(size = 0.8) +
    facet_wrap(~serie, scales = "free_y") +
    scale_color_manual(values = c("#b2182b", "#5b8db8")) +
    labs(x = "meses desde la entrada del competidor", y = ylab, color = NULL,
         title = "series crudas alrededor de la entrada, sin efectos fijos",
         subtitle = subt) +
    theme_minimal() + theme(legend.position = "bottom")
}

ggsave(file.path(FIG, "descriptivo_niveles.pdf"),
       graficar("media", "media del grupo ($/L)",
                paste("medias simples sobre pares evento-estacion balanceados;",
                      "la brecha de nivel es lo que absorbe el ef de estacion")),
       width = 8.5, height = 5.4)

ggsave(file.path(FIG, "descriptivo_normalizado.pdf"),
       graficar("norm", "cambio respecto de tau = -1 ($/L)",
                paste("cada serie recentrada en su propio valor en tau = -1,",
                      "la misma referencia del event study")),
       width = 8.5, height = 5.4)

cat("\n=== NIVELES EN tau = -1 (la brecha que absorbe el ef de estacion) ===\n")
print(dcast(medias[event_time == -1L], grupo ~ serie, value.var = "media")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 1) else z)])

cat("\n=== CAMBIO ACUMULADO RESPECTO DE tau = -1 ===\n")
print(dcast(medias[event_time %in% c(-12L, -6L, 6L, 12L)],
            serie + grupo ~ event_time, value.var = "norm")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 2) else z)])

# ==============================================================================
# within-event difference: the descriptive counterpart of the event study
#
# the level plot above is not a like-for-like comparison. the 225 events that
# contribute treated pairs and the 394 that contribute controls fall on
# different calendar mixes, and prices track the national mepco cost, so the two
# series drift apart for reasons that have nothing to do with entry. taking the
# treated-minus-control difference INSIDE each event and only then averaging
# removes that, because within an event both arms face the same dates. this is
# exactly the job that the event_id^ym fixed effect does in the regression
# ==============================================================================

por_evento <- largo[, .(valor = mean(valor)),
                    by = .(event_id = sub(" .*", "", par), outcome, event_time,
                           treated)]
por_evento <- dcast(por_evento, event_id + outcome + event_time ~ treated,
                    value.var = "valor")
setnames(por_evento, c("0", "1"), c("control", "tratada"))
por_evento <- por_evento[!is.na(control) & !is.na(tratada)]
# keep events observed on both arms at every tau, so the average is over a
# fixed set of experiments
ok <- por_evento[, .N, by = .(event_id, outcome)][N == 2L * VENTANA + 1L]
por_evento <- por_evento[ok, on = .(event_id, outcome)]
por_evento[, dif := tratada - control]

cat("\n=== EVENTOS CON AMBOS BRAZOS COMPLETOS ===\n")
print(por_evento[, .(eventos = uniqueN(event_id)), by = outcome])

did <- por_evento[, .(dif = mean(dif), n = .N), by = .(outcome, event_time)]
did[, norm := dif - dif[event_time == -1L], by = outcome]
fwrite(did, file.path(TAB, "descriptivo_did.csv"))

did[, serie := factor(outcome, levels = names(outcomes), labels = outcomes)]

p_did <- ggplot(did, aes(event_time, norm)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey55") +
  geom_line(color = "#2166ac", linewidth = 0.55) +
  geom_point(color = "#2166ac", size = 0.8) +
  facet_wrap(~serie, scales = "free_y") +
  labs(x = "meses desde la entrada del competidor",
       y = "tratadas menos control, respecto de tau = -1 ($/L)",
       title = "diagnostico: diferencia cruda, sin efecto de tiempo por experimento",
       subtitle = paste("no es evidencia de tendencias paralelas; la deriva",
                        "previa desaparece al incluir event_id x mes (seccion 3)")) +
  theme_minimal()
ggsave(file.path(FIG, "descriptivo_did.pdf"), p_did, width = 8.5, height = 5.4)

cat("\n=== DIFERENCIA CRUDA RESPECTO DE tau = -1 ===\n")
print(dcast(did[event_time %in% c(-12L, -6L, 6L, 12L)], serie ~ event_time,
            value.var = "norm")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 2) else z)])

# ==============================================================================
# 3. what closes the gap between the raw series and the estimates
#
# the same balanced sample, adding fixed effects one at a time, reporting the
# slope of the pre-period path. the crude difference in section 2 is the first
# rung with no time effect at all
# ==============================================================================

library(fixest)

stk_fe <- copy(bal)
stk_fe[, year := year(as.IDate(ym))]

pend_es <- function(m, sc) {
  b <- coef(m); V <- vcov(m)
  nm <- grep("^event_time::-", names(b), value = TRUE)
  tt <- as.integer(sub("^event_time::(-?[0-9]+).*", "\\1", nm))
  o <- order(tt); nm <- nm[o]; tt <- tt[o]
  w <- c(tt, -1L) - mean(c(tt, -1L))
  w <- (w / sum(w^2))[-(length(tt) + 1L)]
  e <- sum(w * b[nm]) * sc
  se <- sqrt(as.numeric(t(w) %*% V[nm, nm] %*% w)) * sc
  data.table(pend = e, p = 2 * pnorm(-abs(e / se)))
}

fes <- c(`1. solo evento x estacion` = "event_id^station_key",
         `2. + evento x mes` = "event_id^station_key + event_id^ym",
         `3. + marca x anio (principal)` =
           "event_id^station_key + event_id^ym + distribuidor^year")

escalera <- rbindlist(lapply(names(fes), function(nm) {
  rbindlist(lapply(names(outcomes), function(fv) {
    en_log <- substr(fv, 1, 1) == "p"
    sc  <- if (en_log) 100 else 1
    lhs <- if (en_log) sprintf("log(%s)", fv) else fv
    dd <- stk_fe[!is.na(get(fv))]
    m <- feols(as.formula(sprintf(
      "%s ~ i(event_time, treated, ref = -1) | %s", lhs, fes[[nm]])),
      data = dd, cluster = ~comuna)
    a <- coeftable(feols(as.formula(sprintf("%s ~ treated:post | %s",
                                            lhs, fes[[nm]])),
                         data = dd, cluster = ~comuna))
    cbind(data.table(fe = nm, outcome = fv, att = a[1, 1] * sc),
          pend_es(m, sc))
  }))
}))
fwrite(escalera, file.path(TAB, "descriptivo_escalera_fe.csv"))

cat("\n=== QUE APLANA LA PRE-TENDENCIA (misma muestra balanceada) ===\n")
cat("pendiente pre por mes:\n")
print(dcast(escalera, fe ~ outcome, value.var = "pend")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 4) else z)])
cat("\np de la pendiente:\n")
print(dcast(escalera, fe ~ outcome, value.var = "p")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 4) else z)])
cat("\natt (ojo con el signo en el peldano 1: sin efecto de tiempo el\n")
cat("estimador recoge la trayectoria del costo mepco, no la entrada):\n")
print(dcast(escalera, fe ~ outcome, value.var = "att")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

message("15_tendencias_descriptivas.R: figuras en ", FIG, ".")
