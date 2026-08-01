
# objective: the effect of entry along the price distribution, through the
# recentred influence function of \textcite{FIRPO2009}. this is the route
# fischer, martin & schmidt-dengler (2025) take in their appendix, and the one
# that remains available here after the panel-data quantile estimators of
# callaway & li proved underpowered on this sample: rif imposes no balanced
# panel requirement, so it uses all 734 treated stations rather than the 204
# that survive balancing.
#
# WHAT RIF DOES. for a quantile tau, the recentred influence function of the
# outcome is
#     RIF(y; q_tau) = q_tau + (tau - 1{y <= q_tau}) / f(q_tau),
# whose conditional expectation is the unconditional quantile. regressing it on
# the treatment therefore recovers the effect of entry on the tau-quantile of
# the UNCONDITIONAL distribution -- the distribution consumers actually face --
# and not the conditional quantile a quantile regression would give.
#
# ITS LIMIT, WHICH MATTERS FOR HOW THE RESULT IS READ. rif is a local linear
# approximation around the observed distribution. it is informative quantile by
# quantile, but it does not licence a statement about the global shift of the
# distribution: fischer et al. need a second method for that, and it is
# precisely that global statement which this sample cannot support.
#
# THE OUTCOME. a note from the earlier version of this project makes a point
# that is easy to miss and decisive here. the pooled distribution of raw chilean
# prices over 2015-2026 reflects mostly TIME -- inflation, exchange rate, crude
# -- and not the cross-sectional dispersion between stations that the quantiles
# are meant to describe; in its tails the density f(q_tau) becomes tiny and the
# variance of the rif explodes. germany over 2015-2020 has no such problem.
# that version subtracted the national monthly mean. here the margin does the
# same job by construction and with an economic meaning rather than an ad hoc
# demeaning, since mepco is a national weekly cost common to every station. both
# are reported: the detrended price keeps fischer's units and comparability, the
# margin is the cleaner object.
#
# TAKES:  data/processed/panel_mensual.csv, entradas.csv
# PRODUCES:
#   results/tables/qte_rif.csv, tab_qte_rif.tex
#   results/figures/qte_rif_precio.pdf, qte_rif_margen.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

RADIOS <- c(1, 2)        # fischer's two radii, the outer one being this
                         # project's own treatment radius
VENT_M <- 24L            # months kept around the entry for treated stations
FOCAL_FROM <- 2014L
TAUS <- seq(0.10, 0.90, by = 0.05)
R_EARTH <- 6371

FE <- "station_key + ym + region^year + distribuidor^year"

radio_lbl <- function(r) sprintf("Entrada a [0; %s] km", sub("\\.", ",", r))

mi <- function(d) year(d) * 12L + (month(d) - 1L)

celda_tex <- function(est, se, p, dig = 2) {
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

# ==============================================================================
# 1. treatment at each radius
# ==============================================================================

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), year = year(as.IDate(ym)))]
entradas <- fread("data/processed/entradas.csv")[, g := as.IDate(g)]

sloc <- unique(panel[!is.na(lat), .(station_key, lat, lon)])[, .SD[1], by = station_key]
span <- panel[, .(first_mi = min(miym), last_mi = max(miym)), by = station_key]
sloc <- merge(sloc, span, by = "station_key")

a1 <- sloc$lat * pi / 180; o1 <- sloc$lon * pi / 180
a2 <- entradas$elat * pi / 180; o2 <- entradas$elon * pi / 180
h <- sin(outer(a1, a2, "-") / 2)^2 +
  outer(cos(a1), cos(a2)) * sin(outer(o1, o2, "-") / 2)^2
h[h > 1] <- 1
D_ev <- 2 * R_EARTH * asin(sqrt(h))

# a station is never treated by its own entry nor by one predating its arrival
own <- match(entradas$station_key, sloc$station_key)
D_ev[cbind(own, seq_len(nrow(entradas)))] <- Inf
inc <- outer(sloc$first_mi, mi(entradas$g), "<") &
  outer(sloc$last_mi, mi(entradas$g), ">=")
D_ev[!inc] <- Inf
g_ev <- mi(entradas$g)

# the analysis sample for a radius. treated stations are held to a window around
# their entry, because a static did should not average horizons of up to a
# decade; the never-treated keep their whole series and anchor the time effects
muestra_R <- function(R) {
  near <- lapply(seq_len(nrow(sloc)), function(i) {
    w <- which(D_ev[i, ] <= R)
    if (!length(w)) NULL else sort(unique(g_ev[w]))
  })
  asg <- data.table(
    station_key = sloc$station_key,
    g1 = vapply(near, function(x) if (is.null(x)) NA_integer_ else x[1], integer(1)),
    g2 = vapply(near, function(x) if (length(x) >= 2L) x[2] else NA_integer_,
                integer(1)))

  d <- merge(panel, asg, by = "station_key")
  d[, trat := as.integer(!is.na(g1) & g1 %/% 12L >= FOCAL_FROM)]
  # a station whose first nearby entry predates 2014 is neither treated nor a
  # control: its status in the window is ambiguous
  d <- d[trat == 1L | is.na(g1)]
  # never-treated controls exclude post-2012 entrants, whose own opening path is
  # no counterfactual
  d <- d[trat == 1L | base_2012 == TRUE]
  d <- d[is.na(g2) | miym < g2]
  d[, dmi := miym - g1]
  d <- d[trat == 0L | (dmi >= -VENT_M & dmi <= VENT_M)]
  d[, post := as.integer(trat == 1L & miym >= g1)]
  d[]
}

# ==============================================================================
# 2. rif regressions
#
# the density is estimated once per outcome on the pooled sample, as in firpo,
# fortin & lemieux. the standard errors treat the quantile and the density as
# known, which is the usual practice and what fischer et al. do
# ==============================================================================

# national monthly mean over the WHOLE panel, not over the analysis sample, so
# the benchmark neither depends on the radius nor is contaminated by the treated
# stations after their entry
media_nac <- function(fv) panel[is.finite(get(fv)), .(mm = mean(get(fv))), by = ym]

curva_rif <- function(d, fv, detrend) {
  d <- d[is.finite(get(fv))]
  if (detrend) {
    d <- merge(d, media_nac(fv), by = "ym")
    gm <- mean(d[[fv]])
    d[, yv := get(fv) - mm + gm]      # keeps $/litre units
  } else {
    d[, yv := get(fv)]
  }
  y <- d$yv
  dens <- density(y)
  rbindlist(lapply(TAUS, function(tau) {
    qt <- as.numeric(quantile(y, tau))
    fhat <- approx(dens$x, dens$y, xout = qt)$y
    d[, rif := qt + (tau - as.numeric(yv <= qt)) / fhat]
    m <- feols(as.formula(paste0("rif ~ post | ", FE)), data = d,
               cluster = ~comuna)
    ct <- coeftable(m)["post", ]
    data.table(tau = tau, est = ct[1], se = ct[2], p = ct[4],
               n_obs = nobs(m), n_trat = uniqueN(d[trat == 1L, station_key]))
  }))
}

salidas <- list(
  precio = list(vars = c(p93 = "Gasolina 93", pdi = "Diésel"), detrend = TRUE,
                ylab = "Efecto sobre el precio ($/litro)",
                sub = "precio neto del promedio nacional del mes"),
  margen = list(vars = c(m93 = "Gasolina 93", mdi = "Diésel"), detrend = FALSE,
                ylab = "Efecto sobre el margen ($/litro)",
                sub = "margen sobre el precio mayorista con MEPCO")
)

curvas <- rbindlist(lapply(RADIOS, function(R) {
  d <- muestra_R(R)
  cat(sprintf("radio %.1f km: %d tratadas | %d nunca tratadas\n", R,
              uniqueN(d[trat == 1L, station_key]),
              uniqueN(d[trat == 0L, station_key])))
  rbindlist(lapply(names(salidas), function(tp) {
    cfg <- salidas[[tp]]
    rbindlist(lapply(names(cfg$vars), function(fv) {
      curva_rif(d, fv, cfg$detrend)[, `:=`(radio = R, tipo = tp, outcome = fv)]
    }))
  }))
}))
curvas[, `:=`(ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
              radio_f = factor(radio_lbl(radio), levels = radio_lbl(RADIOS)))]

fwrite(curvas, file.path(TAB, "qte_rif.csv"))

for (tp in names(salidas)) {
  cat(sprintf("\n=== EFECTO CUANTILICO SOBRE EL %s ($/litro) ===\n", toupper(tp)))
  print(dcast(curvas[tipo == tp & tau %in% c(.1, .25, .5, .75, .9)],
              tau ~ outcome + radio, value.var = "est")[
    , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
}

cat("\n=== SIGNIFICANCIA POR CURVA ===\n")
print(curvas[, .(cuantiles = .N, signif_5 = sum(p < 0.05),
                 est_min = round(min(est), 2), est_max = round(max(est), 2),
                 se_mediana = round(median(se), 2)),
             by = .(tipo, outcome, radio)][order(tipo, outcome, radio)])

# ==============================================================================
# 3. table and figures
# ==============================================================================

sel <- curvas[tau %in% c(0.10, 0.25, 0.50, 0.75, 0.90)]
w <- dcast(sel, tau ~ tipo + outcome + radio, value.var = c("est", "se", "p"))
out <- data.table(Cuantil = sprintf("$q_{%.0f}$", 100 * w$tau))
for (tp in names(salidas)) for (fv in names(salidas[[tp]]$vars)) for (R in RADIOS) {
  k <- paste(tp, fv, R, sep = "_")
  out[[sprintf("%s %s, %s km", if (tp == "precio") "Precio" else "Margen",
               salidas[[tp]]$vars[[fv]], sub("\\.", ",", R))]] <-
    celda_tex(w[[paste0("est_", k)]], w[[paste0("se_", k)]], w[[paste0("p_", k)]])
}
guardar_tabla_tex(out, "tab_qte_rif.tex")

pal <- setNames(c(AZUL, NARANJO_OSC), radio_lbl(RADIOS))
for (tp in names(salidas)) {
  cfg <- salidas[[tp]]
  d <- curvas[tipo == tp]
  d[, combustible := factor(outcome, levels = names(cfg$vars), labels = cfg$vars)]
  p <- ggplot(d, aes(tau, est, colour = radio_f, shape = radio_f,
                     linetype = radio_f)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.012,
                  linetype = "solid",
                  position = position_dodge(width = 0.015)) +
    geom_line(position = position_dodge(width = 0.015)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.015)) +
    facet_wrap(~combustible, scales = "free_y") +
    scale_colour_manual(values = pal) +
    scale_shape_manual(values = c(18, 17)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
    labs(x = "cuantil de la distribución (τ)", y = cfg$ylab,
         colour = NULL, shape = NULL, linetype = NULL,
         title = "efecto de la entrada a lo largo de la distribución",
         subtitle = paste0(cfg$sub, "; regresiones RIF, ee cluster comuna")) +
    tema()
  guardar(p, file.path(FIG, sprintf("qte_rif_%s.pdf", tp)), 9.5, 5)
}

message("22_qte_rif.R: figuras en ", FIG, " y tablas en ", TAB, ".")
