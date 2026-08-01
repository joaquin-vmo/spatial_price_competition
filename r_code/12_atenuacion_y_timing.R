
# objective: the two checks that turn the timing assumption from "defensible"
# into "defended". both use data already built by 10_build_event_panel.R.
#
#   A. distance attenuation. the main design leans on comparing incumbents
#      within 2 km against far controls (>5 km), and that estimate is much
#      larger than the one using ring controls. two readings compete: the ring
#      is partially treated (spillover), or the far controls are systematically
#      more rural and the effect is an urban/rural trend. a smooth decay from
#      0-1 km outwards is only consistent with the first. a flat profile in
#      distance would mean composition, not entry.
#      this is a port of section D of 05-estimate-eventstudy.R (old project),
#      kept faithful: quarterly panel, +-4 quarters, same fixed effects and
#      same estimator; only the panel build underneath is the new one.
#
#   B. entry-timing hazard. the sharpest available test of the assumption
#      itself. among markets that EVENTUALLY receive entry, do lagged local
#      conditions predict WHEN it arrives? the contrast is the point:
#        - incidence sample: every incumbent at risk. local profitability
#          should predict entry -- that is location endogeneity, and the
#          station fixed effect absorbs it.
#        - timing sample: only stations that do get entry. here the same
#          covariates must be uninformative, because that is exactly what the
#          maintained assumption claims.
#
# TAKES:  data/processed/panel_mensual.csv  (from 10_build_event_panel.R)
# PRODUCES:
#   results/tables/atenuacion_distancia.csv, atenuacion_distancia_fina.csv
#   results/tables/hazard_timing.csv
#   results/figures/atenuacion_distancia.pdf, atenuacion_distancia_fina.pdf,
#     hazard_timing.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

RTREAT <- 2
RCTRL  <- 5
FOCAL_FROM <- 2014L
COHORTE_NUNCA <- 1000000L
R_EARTH <- 6371

precios <- c(p93 = "Gasolina 93", p95 = "Gasolina 95",
             p97 = "Gasolina 97", pdi = "Diésel")

trim_i <- function(d) year(d) * 4L + (month(d) - 1L) %/% 3L

tidy_coefs <- function(m) {
  ct <- coeftable(m)
  data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2], p = ct[, 4])
}

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g5_entry = as.IDate(g5_entry), g2_5_entry = as.IDate(g2_5_entry))]
panel[, `:=`(year = year(ym), period_q = trim_i(ym))]

# ==============================================================================
# A. DISTANCE ATTENUATION
#
# every station is assigned to the ring of its FIRST entry within RCTRL km;
# control = never within RCTRL of any entry. the coefficient on each ring is
# the static did effect for that band, so the profile across bands is the
# spatial decay. quarterly for tractability, as in the old build
# ==============================================================================

anillo_lbl <- c("0-1 km", "1-2 km", "2-3 km", "3-4 km", "4-5 km")

armar_anillos <- function(fv, dmax = NULL, brk = NULL, lbl = NULL) {
  d <- panel[role5_entry %in% c("treated", "far") &
               (is.na(g2_5_entry) | ym < g2_5_entry)]
  # a station whose first nearby entry is a suspect 2013 event never treats
  d <- d[!(role5_entry == "treated" & year(g5_entry) < FOCAL_FROM)]
  # control pool excludes the entrants themselves: their own post-opening path
  # is not a valid counterfactual
  d <- d[role5_entry == "treated" | base_2012 == TRUE]
  if (!is.null(dmax)) d <- d[role5_entry == "far" | dist_entry5 <= dmax]

  d[, `:=`(cohort_q = fifelse(role5_entry == "treated", trim_i(g5_entry),
                              COHORTE_NUNCA),
           lp = log(get(fv)))]
  d[, rel_q := period_q - cohort_q]
  d <- d[(role5_entry == "far" | (rel_q >= -4 & rel_q <= 4)) & is.finite(lp)]

  dq <- d[, .(lp = mean(lp), role5 = first(role5_entry),
              ring = first(ring_entry), dist = first(dist_entry5),
              cohort_q = first(cohort_q), region = first(region),
              distribuidor = first(distribuidor), comuna = first(comuna),
              year = first(year)),
          by = .(station_key, period_q)]
  dq[, post := as.integer(role5 == "treated" & period_q >= cohort_q)]
  if (!is.null(brk)) {
    dq[, ring := cut(dist, breaks = brk, labels = lbl, include.lowest = TRUE)]
  }
  dq[, tpr := factor(fifelse(post == 1L, as.character(ring), "none"),
                     levels = c("none", if (is.null(lbl)) anillo_lbl else lbl))]
  dq[]
}

att_por_anillo <- function(dq, lbl) {
  m <- feols(lp ~ i(tpr, ref = "none") |
               station_key + region^period_q + distribuidor^year,
             data = dq, cluster = ~comuna)
  r <- tidy_coefs(m)[grepl("tpr::", term)]
  r[, ring := factor(sub(".*tpr::", "", term), levels = lbl)]
  r[, .(ring, att = estimate * 100, se = se * 100, p,
        ci_low = (estimate - 1.96 * se) * 100,
        ci_high = (estimate + 1.96 * se) * 100)]
}

anillos <- rbindlist(lapply(names(precios), function(fv) {
  att_por_anillo(armar_anillos(fv), anillo_lbl)[, outcome := fv]
}))

# finer 500 m bands inside 0-3 km, where the action is
fino_brk <- seq(0, 3, by = 0.5)
fino_lbl <- paste0(head(fino_brk, -1) * 1000, "-", fino_brk[-1] * 1000, " m")

anillos_fino <- rbindlist(lapply(names(precios), function(fv) {
  dq <- armar_anillos(fv, dmax = 3, brk = fino_brk, lbl = fino_lbl)
  att_por_anillo(dq, fino_lbl)[, outcome := fv]
}))

fwrite(anillos,      file.path(TAB, "atenuacion_distancia.csv"))
fwrite(anillos_fino, file.path(TAB, "atenuacion_distancia_fina.csv"))

cat("\n=== ATT POR ANILLO DE DISTANCIA (entrada, % del precio) ===\n")
cat("si el efecto decae de 0-1 km hacia cero, el control anillo TIENE que atenuar\n")
print(dcast(anillos, ring ~ outcome, value.var = "att")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

cat("\n=== ATT POR ANILLO DE 500 m (0-3 km) ===\n")
print(dcast(anillos_fino, ring ~ outcome, value.var = "att")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

plot_anillos <- function(d, subt) {
  d <- copy(d)[, combustible := factor(outcome, levels = names(precios),
                                       labels = precios)]
  ggplot(d, aes(ring, att)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                  color = AZUL) +
    geom_line(aes(group = 1), color = AZUL) +
    geom_point(color = AZUL, size = 1.6) +
    facet_wrap(~combustible, scales = "free_y") +
    labs(x = "Distancia de la incumbente a la entrada",
         y = "Efecto sobre el precio (%)",
         title = "Atenuación del efecto de la entrada con la distancia",
         subtitle = subt) +
    tema()
}

guardar(plot_anillos(anillos,
                     "Cada estación se asigna al anillo de su primera entrada cercana"),
        file.path(FIG, "atenuacion_distancia.pdf"), 9, 5.4)
guardar(plot_anillos(anillos_fino,
                     "Bandas de 500 m; tratadas restringidas a menos de 3 km") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)),
        file.path(FIG, "atenuacion_distancia_fina.pdf"), 9, 5.8)

# ==============================================================================
# B. ENTRY-TIMING HAZARD
# ==============================================================================

# --- local market conditions: neighbours within RTREAT km ---------------------
sloc <- unique(panel[!is.na(lat), .(station_key, lat, lon)])
a1 <- sloc$lat * pi / 180; o1 <- sloc$lon * pi / 180
dl <- outer(a1, a1, "-"); do_ <- outer(o1, o1, "-")
h <- sin(dl / 2)^2 + outer(cos(a1), cos(a1)) * sin(do_ / 2)^2
h[h > 1] <- 1
D <- 2 * R_EARTH * asin(sqrt(h))
idx <- which(D <= RTREAT, arr.ind = TRUE)   # includes the diagonal: own market
edges <- data.table(station_key = sloc$station_key[idx[, "row"]],
                    vecina      = sloc$station_key[idx[, "col"]])

# station-quarter outcomes, then averaged over the local market
pq <- panel[, .(lp93 = mean(log(p93), na.rm = TRUE),
                mar93 = mean(m93, na.rm = TRUE)),
            by = .(station_key, period_q)]
pq <- pq[is.finite(lp93)]

loc <- merge(edges, pq, by.x = "vecina", by.y = "station_key",
             allow.cartesian = TRUE)
loc <- loc[, .(mar_local = mean(mar93, na.rm = TRUE),
               lp_local  = mean(lp93,  na.rm = TRUE),
               ncomp     = .N - 1L),          # exclude self from the count
           by = .(station_key, period_q)]

hz <- merge(pq, loc, by = c("station_key", "period_q"))
setorder(hz, station_key, period_q)
# lagged one quarter: conditions the entrant could have observed before opening
hz[, `:=`(l_mar_local = shift(mar_local), l_lp_local = shift(lp_local),
          l_ncomp = shift(ncomp), l_mar_own = shift(mar93)),
   by = station_key]
# four-quarter change in local margin: does entry time to IMPROVING markets?
hz[, l_dmar_local := l_mar_local - shift(mar_local, 5L), by = station_key]

# --- risk set -----------------------------------------------------------------
st <- unique(panel[, .(station_key, role_entry, g_entry, base_2012, region,
                       comuna, distribuidor)])
st <- st[, .SD[1], by = station_key]
st[, g_q := trim_i(g_entry)]
# a suspect 2013 first entry is neither an event nor a valid censoring point
st <- st[is.na(g_entry) | year(g_entry) >= FOCAL_FROM]

hz <- merge(hz, st, by = "station_key")
hz[, first_q := min(period_q), by = station_key]
# at risk from the quarter after the station itself starts reporting, until it
# experiences its first nearby entry (inclusive) or the panel ends
hz <- hz[period_q > first_q & (is.na(g_q) | period_q <= g_q)]
hz[, entrada := as.integer(!is.na(g_q) & period_q == g_q)]
hz[, year := period_q %/% 4L]

covs <- c("l_mar_local", "l_dmar_local", "l_ncomp", "l_lp_local")
hz <- hz[complete.cases(hz[, ..covs])]
# standardised so the coefficients are comparable across covariates and samples
hz[, (covs) := lapply(.SD, function(x) x / sd(x)), .SDcols = covs]

muestras <- list(
  incidencia = hz,                                   # every incumbent at risk
  timing     = hz[!is.na(g_q)]                       # only those that do get entry
)

fml_hz <- as.formula(paste("entrada ~", paste(covs, collapse = " + "),
                           "| region^year"))

# the joint test is split into two groups, because they are different stories.
# RENTABILIDAD is the classic endogeneity worry: entry timed to how profitable
# the local market currently is. ESTRUCTURA is densification: markets that are
# already adding stations keep adding them, which shifts timing without saying
# anything about margins. only the first threatens the price estimate directly
grupos <- list(rentabilidad = c("l_mar_local", "l_dmar_local", "l_lp_local"),
               estructura   = "l_ncomp",
               conjunto     = covs)

modelos <- lapply(muestras, function(d)
  feglm(fml_hz, data = d, family = binomial(), cluster = ~comuna))

res_hz <- rbindlist(lapply(names(modelos), function(nm) {
  d <- muestras[[nm]]
  tidy_coefs(modelos[[nm]])[, .(muestra = nm, term, estimate, se, p,
                                or = exp(estimate), n_obs = nobs(modelos[[nm]]),
                                n_est = uniqueN(d$station_key),
                                n_eventos = sum(d$entrada))]
}))

tests_hz <- rbindlist(lapply(names(modelos), function(nm) {
  rbindlist(lapply(names(grupos), function(g) {
    w <- wald(modelos[[nm]], paste(grupos[[g]], collapse = "|"), print = FALSE)
    data.table(muestra = nm, grupo = g, F = w$stat, df = w$df1, p = w$p)
  }))
}))

fwrite(res_hz,   file.path(TAB, "hazard_timing.csv"))
fwrite(tests_hz, file.path(TAB, "hazard_timing_tests.csv"))

cat("\n=== HAZARD DE ENTRADA: INCIDENCIA vs TIMING ===\n")
cat("coeficientes logit por desviacion estandar del regresor; ee cluster comuna\n")
print(res_hz[, .(muestra, term, coef = round(estimate, 3),
                 se = round(se, 3), p = round(p, 4),
                 odds_ratio = round(or, 3))])

cat("\n=== TESTS CONJUNTOS (H0: el grupo no predice la entrada) ===\n")
cat("el supuesto exige que RENTABILIDAD no prediga el timing.\n")
cat("que ESTRUCTURA si lo haga es densificacion, no seleccion sobre margenes\n")
print(tests_hz[, .(muestra, grupo, F = round(F, 2), df, p = round(p, 4))])
print(unique(res_hz[, .(muestra, n_obs, n_est, n_eventos)]))

lbl_cov <- c(l_mar_local = "Margen local (t-1)",
             l_dmar_local = "Cambio del margen local (4 trim.)",
             l_ncomp = "N° de competidores (t-1)",
             l_lp_local = "Log precio local (t-1)")
lbl_mue <- c(incidencia = "Incidencia: todos los incumbentes en riesgo",
             timing = "Timing: solo los que sí reciben entrada")

d_hz <- copy(res_hz)
d_hz[, `:=`(cov = factor(lbl_cov[term], levels = rev(lbl_cov)),
            muestra = factor(lbl_mue[muestra], levels = lbl_mue),
            ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

p_hz <- ggplot(d_hz, aes(estimate, cov)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.15, color = AZUL) +
  geom_point(color = AZUL, size = 1.8) +
  facet_wrap(~muestra) +
  labs(x = "Efecto sobre el log-odds de entrada, por desviación estándar",
       y = NULL,
       title = "Las condiciones locales predicen dónde entra un competidor, no cuándo",
       subtitle = "Logit de tiempo discreto, trimestral") +
  tema()
guardar(p_hz, file.path(FIG, "hazard_timing.pdf"), 9, 3.6)

message("12_atenuacion_y_timing.R: figuras en ", FIG, " y tablas en ", TAB, ".")
