
# objective: two robustness exercises for the annex.
#
#   A. staggered difference-in-differences. with 403 entry cohorts spread over
#      twelve years, the two-way fixed effects estimator is a weighted average
#      of cohort-time comparisons in which already-treated units act as controls
#      for later-treated ones, and the weights can be negative when treatment
#      effects vary across cohorts \parencite{SUNABRAHAM2021}. the
#      interaction-weighted estimator of sun & abraham is estimated here on the
#      same sample and reported next to the twfe one. it is the natural choice
#      in this setting because the design already has a large never-treated
#      pool, which is what the estimator uses as its comparison group.
#
#   B. persistence of the treatment. the design assumes that an entry is a
#      discrete and lasting increase in the number of competitors of the local
#      market. two event studies test it, following figure 2 of fischer, martin
#      & schmidt-dengler (2025): the count of active stations within the
#      treatment radius should jump by one and stay there, and the count that
#      EXCLUDES the entrant should stay flat. the second is the substantive one:
#      if entry pushed incumbents out, the estimated price effect would be
#      partly a composition effect rather than a competitive response.
#
# TAKES:  data/processed/panel_mensual.csv
# PRODUCES:
#   results/tables/sunab_att.csv, sunab_coefs.csv, tab_sunab.tex
#   results/tables/persistencia_coefs.csv
#   results/figures/sunab_vs_twfe.pdf, persistencia.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

RTREAT <- 2
FOCAL_FROM <- 2014L
NSEM  <- 4L          # semesters kept each side; endpoints binned
NUNCA <- 1000000L    # cohort code for the never treated
R_EARTH <- 6371

precios <- c(p93 = "gasolina 93", p95 = "gasolina 95",
             p97 = "gasolina 97", pdi = "diesel")

mi <- function(d) year(d) * 12L + (month(d) - 1L)
sem_i <- function(m) m %/% 6L          # continuous six-month index

celda_tex <- function(est, se, p, dig = 3) {
  st <- fifelse(p < 0.01, "^{***}",
                fifelse(p < 0.05, "^{**}", fifelse(p < 0.10, "^{*}", "")))
  sprintf("$%.*f%s$ ($%.*f$)", dig, est, st, dig, se)
}

guardar_tabla_tex <- function(d, archivo) {
  cuerpo <- apply(d, 1, function(r) paste(paste(r, collapse = " & "), "\\\\"))
  writeLines(
    c(sprintf("\\begin{tabular}{l%s}", strrep("c", ncol(d) - 1L)), "\\toprule",
      paste(paste(names(d), collapse = " & "), "\\\\"), "\\midrule",
      cuerpo, "\\bottomrule", "\\end{tabular}"),
    file.path(TAB, archivo)
  )
}

etiquetas <- c(sprintf("≤−%d", NSEM),
               as.character(-(NSEM - 1):(NSEM - 1)),
               sprintf("≥%d", NSEM))

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g2_entry = as.IDate(g2_entry))]
panel[, year := year(ym)]

# ==============================================================================
# A. SUN & ABRAHAM
#
# the estimator saturates the regression in cohort x relative-period cells and
# then aggregates with the cohort shares, so it never uses an already-treated
# unit as a control. it is estimated on a six-month panel: with 403 monthly
# cohorts the saturated regression is unwieldy, and six months is also the bin
# width of the reported event studies
# ==============================================================================

d <- panel[role_entry %in% c("treated", "ring", "far")]
d <- d[is.na(g2_entry) | ym < g2_entry]
d <- d[!(role_entry == "treated" & year(g_entry) < FOCAL_FROM)]
d <- d[role_entry == "treated" | base_2012 == TRUE]
d[, periodo := sem_i(miym)]
d[, cohorte := fifelse(role_entry == "treated", sem_i(mi(g_entry)), NUNCA)]

colapsar <- function(fv, en_log) {
  x <- d[!is.na(get(fv))]
  x[, y := if (en_log) log(get(fv)) else get(fv)]
  x <- x[is.finite(y)]
  x[, .(y = mean(y), cohorte = first(cohorte), region = first(region),
        comuna = first(comuna), distribuidor = first(distribuidor),
        year = first(year)),
    by = .(station_key, periodo)]
}

FE_SEM <- "station_key + periodo + region^year + distribuidor^year"

corre_sa <- function(fv) {
  en_log <- substr(fv, 1, 1) == "p"
  sc <- if (en_log) 100 else 1
  dq <- colapsar(fv, en_log)
  dq[, treated := as.integer(cohorte < NUNCA)]
  dq[, rel := fifelse(treated == 1L,
                      pmax(-NSEM, pmin(NSEM, periodo - cohorte)), -1L)]
  # sunab needs the never treated flagged with an infinite cohort
  dq[, coh_sa := fifelse(cohorte == NUNCA, 10000L, as.integer(cohorte))]

  m_sa <- feols(as.formula(sprintf("y ~ sunab(coh_sa, periodo) | %s", FE_SEM)),
                data = dq, cluster = ~comuna)
  m_tw <- feols(as.formula(sprintf("y ~ i(rel, treated, ref = -1) | %s", FE_SEM)),
                data = dq, cluster = ~comuna)

  saca <- function(m, patron, etiqueta) {
    ct <- coeftable(m)
    r <- data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2])
    r <- r[grepl(patron, term)]
    r[, rel := as.integer(sub(paste0(".*", patron, "(-?[0-9]+).*"), "\\1", term))]
    r <- rbind(r[, .(rel, estimate, se)], data.table(rel = -1L, estimate = 0, se = 0))
    r[, `:=`(estimate = estimate * sc, se = se * sc, estimador = etiqueta,
             outcome = fv)]
    r[rel >= -NSEM & rel <= NSEM][order(rel)]
  }
  es <- rbind(saca(m_sa, "periodo::", "Sun-Abraham (IW)"),
              saca(m_tw, "rel::", "TWFE"))
  es[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

  at <- summary(m_sa, agg = "att")$coeftable
  aw <- coeftable(feols(as.formula(sprintf(
    "y ~ i(periodo >= cohorte, treated, ref = FALSE) | %s", FE_SEM)),
    data = dq, cluster = ~comuna))
  att <- data.table(outcome = fv,
                    att_sa = at[1, 1] * sc, se_sa = at[1, 2] * sc, p_sa = at[1, 4],
                    att_tw = aw[1, 1] * sc, se_tw = aw[1, 2] * sc, p_tw = aw[1, 4],
                    n_cohortes = uniqueN(dq[treated == 1L, cohorte]),
                    n_nunca = uniqueN(dq[treated == 0L, station_key]))
  list(es = es, att = att)
}

res_sa <- lapply(names(precios), corre_sa)
es_sa  <- rbindlist(lapply(res_sa, `[[`, "es"))
att_sa <- rbindlist(lapply(res_sa, `[[`, "att"))

fwrite(es_sa,  file.path(TAB, "sunab_coefs.csv"))
fwrite(att_sa, file.path(TAB, "sunab_att.csv"))

cat("\n=== ATT: SUN-ABRAHAM vs TWFE (% del precio) ===\n")
print(att_sa[, .(outcome, sa = round(att_sa, 3), se_sa = round(se_sa, 3),
                 twfe = round(att_tw, 3), se_tw = round(se_tw, 3),
                 n_cohortes, n_nunca)])

tab <- data.table(
  Combustible = unname(precios[att_sa$outcome]),
  `Sun--Abraham` = celda_tex(att_sa$att_sa, att_sa$se_sa, att_sa$p_sa),
  `TWFE` = celda_tex(att_sa$att_tw, att_sa$se_tw, att_sa$p_tw))
guardar_tabla_tex(tab, "tab_sunab.tex")

es_sa[, combustible := factor(outcome, levels = names(precios), labels = precios)]
p_sa <- ggplot(es_sa, aes(rel, estimate, colour = estimador, fill = estimador)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.10, colour = NA) +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~combustible, scales = "free_y") +
  scale_x_continuous(breaks = -NSEM:NSEM, labels = etiquetas) +
  scale_colour_manual(values = c("Sun-Abraham (IW)" = AZUL, "TWFE" = NARANJO_OSC)) +
  scale_fill_manual(values = c("Sun-Abraham (IW)" = AZUL, "TWFE" = NARANJO_OSC)) +
  labs(x = "semestres desde la entrada (extremos agrupados)",
       y = "efecto sobre el precio (%)", colour = NULL, fill = NULL,
       title = "estimador de sun-abraham frente a twfe",
       subtitle = "misma muestra y mismos efectos fijos; solo cambia el estimador") +
  tema()
guardar(p_sa, file.path(FIG, "sunab_vs_twfe.pdf"), 9.5, 5.8)

# ==============================================================================
# B. PERSISTENCE OF THE TREATMENT
#
# two counts as outcomes: every active station within RTREAT of the focal one,
# and the same count restricted to stations that were ALREADY operating before
# the focal's own entry event. the first measures whether the treatment happens
# and lasts; the second whether it displaces anybody
# ==============================================================================

sloc <- unique(panel[!is.na(lat), .(station_key, lat, lon)])[, .SD[1], by = station_key]
a1 <- sloc$lat * pi / 180; o1 <- sloc$lon * pi / 180
h <- sin(outer(a1, a1, "-") / 2)^2 +
  outer(cos(a1), cos(a1)) * sin(outer(o1, o1, "-") / 2)^2
h[h > 1] <- 1
D <- 2 * R_EARTH * asin(sqrt(h))
diag(D) <- Inf
ix <- which(D <= RTREAT, arr.ind = TRUE)
edges <- data.table(focal = sloc$station_key[ix[, "row"]],
                    vecina = sloc$station_key[ix[, "col"]])

ini <- panel[, .(vecina_ini = min(miym)), by = .(vecina = station_key)]
edges <- merge(edges, ini, by = "vecina")
presencia <- unique(panel[, .(vecina = station_key, miym)])
ev <- merge(edges, presencia, by = "vecina", allow.cartesian = TRUE)

g_focal <- unique(panel[role_entry == "treated" & !is.na(g_entry),
                        .(focal = station_key, g_mi = mi(g_entry))])
ev <- merge(ev, g_focal, by = "focal", all.x = TRUE)

ncomp <- ev[, .(ncomp = .N), by = .(station_key = focal, miym)]
# neighbours already operating before the focal's treating entry: excludes the
# entrant itself, so this series answers whether anybody was displaced
ninc <- ev[is.na(g_mi) | vecina_ini < g_mi,
           .(ninc = .N), by = .(station_key = focal, miym)]

dp <- panel[role_entry %in% c("treated", "ring", "far")]
dp <- dp[is.na(g2_entry) | ym < g2_entry]
dp <- dp[!(role_entry == "treated" & year(g_entry) < FOCAL_FROM)]
dp <- dp[role_entry == "treated" | base_2012 == TRUE]
dp[, treated := as.integer(role_entry == "treated")]
dp[, rel := fifelse(treated == 1L,
                    pmax(-NSEM, pmin(NSEM,
                      as.integer(floor((miym - mi(g_entry)) / 6L)))), -1L)]
dp <- merge(dp, ncomp, by = c("station_key", "miym"), all.x = TRUE)
dp <- merge(dp, ninc,  by = c("station_key", "miym"), all.x = TRUE)
dp[is.na(ncomp), ncomp := 0L][is.na(ninc), ninc := 0L]

salidas_lbl <- c(ncomp = "Competidores a 2 km (incluye a la entrante)",
                 ninc  = "Competidores a 2 km (excluye a la entrante)")

persist <- rbindlist(lapply(names(salidas_lbl), function(v) {
  m <- feols(as.formula(sprintf(
    "%s ~ i(rel, treated, ref = -1) | station_key + ym + region^year", v)),
    data = dp, cluster = ~comuna)
  ct <- coeftable(m)
  r <- data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2])
  r <- r[grepl("^rel::", term)]
  r[, rel := as.integer(sub("^rel::(-?[0-9]+).*", "\\1", term))]
  r <- rbind(r[, .(rel, estimate, se)], data.table(rel = -1L, estimate = 0, se = 0))
  r[, `:=`(serie = salidas_lbl[[v]], ci_low = estimate - 1.96 * se,
           ci_high = estimate + 1.96 * se)]
  setorder(r, rel)[]
}))
fwrite(persist, file.path(TAB, "persistencia_coefs.csv"))

cat("\n=== PERSISTENCIA DEL TRATAMIENTO (cambio en el n de competidores) ===\n")
print(dcast(persist[rel %in% c(-3, -1, 0, 2, 4)], rel ~ serie,
            value.var = "estimate")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

persist[, serie := factor(serie, levels = salidas_lbl)]
p_per <- ggplot(persist, aes(rel, estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey65") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, fill = AZUL) +
  geom_line(colour = AZUL) + geom_point(colour = AZUL, size = 1.2) +
  facet_wrap(~serie) +
  scale_x_continuous(breaks = -NSEM:NSEM, labels = etiquetas) +
  labs(x = "semestres desde la entrada (extremos agrupados)",
       y = "cambio en el número de estaciones a 2 km",
       title = "persistencia del tratamiento",
       subtitle = "la línea gris marca el salto de una estación") +
  tema()
guardar(p_per, file.path(FIG, "persistencia.pdf"), 9.5, 4.6)

message("19_staggered_persistencia.R: figuras en ", FIG, " y tablas en ", TAB, ".")
