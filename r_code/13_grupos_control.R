
# objective: re-estimate the entry event study under three definitions of the
# control group and pick one on robustness grounds rather than on the size of
# the coefficient. run for two outcomes: the retail price and the gross margin.
#
# what fischer, martin & schmidt-dengler (2025) do: their equation (2) compares
# treated incumbents against ALL non-treated stations, with station and
# state x day fixed effects; they never use a ring design. what they vary in
# robustness (online appendix S.14) is the identification cell -- date,
# county-date, ROR-date instead of state-date -- not the control group. so the
# never-treated pool is the benchmark here, and the other two are the checks.
#
# the three groups:
#   anillo : incumbents 2-5 km from an entry. most comparable (same urban
#            fabric) but partially treated: 12_atenuacion_y_timing.R shows the
#            price effect is still alive at 2-3 km, so this group is
#            contaminated at its inner edge and must attenuate.
#   lejano : stations never within 5 km of any entry. clean of spillover but
#            structurally different -- emptier, more rural markets.
#   nunca  : every station never treated at 2 km (anillo + lejano), i.e. the
#            fischer et al. pool.
#
# the two outcomes:
#   precio : log price, coefficient x100 read as a percentage.
#   margen : price minus the mepco wholesale quote, in $/L and NOT in logs --
#            the margin is a small difference that can sit near zero, so a log
#            transform is fragile. only 93, 97 and diesel: mepco publishes no
#            95 reference, so that panel does not exist for the margin.
#            since mepco is national and weekly, the cost side is common to all
#            stations and the margin isolates the retail pricing decision.
#
# TAKES:  data/processed/panel_mensual.csv  (from 10_build_event_panel.R)
# PRODUCES:
#   results/tables/es_control_coefs.csv, es_control_att.csv,
#     es_control_pretrend.csv, es_control_balance.csv
#   results/tables/tab_att_control{,_pesos,_margen}.tex,
#     tab_pretend_control{,_margen}.tex, tab_balance_control.tex
#   results/figures/es_control_{anillo,lejano,nunca,comparacion}.pdf
#   results/figures/es_control_margen_{anillo,lejano,nunca,comparacion}.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

# six-month bins over a +-4 window, endpoints binned: the reporting format of
# fischer, martin & schmidt-dengler (2025), adopted so the figures are readable
# side by side with theirs
BIN_M <- 6L
NBIN  <- 4L
FOCAL_FROM <- 2014L
RTREAT <- 2
R_EARTH <- 6371

precios  <- c(p93 = "Gasolina 93", p95 = "Gasolina 95",
              p97 = "Gasolina 97", pdi = "Diésel")
margenes <- c(m93 = "Gasolina 93", m97 = "Gasolina 97", mdi = "Diésel")

# outcome definitions: the estimator is the same, only the transform and the
# reporting unit change
salidas <- list(
  precio = list(vars = precios, en_log = TRUE, sufijo = "",
                ylab = "Efecto sobre el precio (%)", dig = 3,
                unidad = "\\%"),
  margen = list(vars = margenes, en_log = FALSE, sufijo = "_margen",
                ylab = "Efecto sobre el margen ($/L)", dig = 2,
                unidad = "\\$/L")
)

grupos <- list(anillo = "ring", lejano = "far", nunca = c("ring", "far"))
grupo_lbl <- c(anillo = "Anillo 2-5 km",
               lejano = "Lejano >5 km",
               nunca  = "Nunca tratadas")
grupo_tex <- c(anillo = "Anillo 2--5 km", lejano = "Lejano $>$5 km",
               nunca = "Nunca tratadas")

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------

tidy_coefs <- function(m) {
  ct <- coeftable(m)
  data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2], p = ct[, 4])
}

tidy_es <- function(m, patron, scale = 100) {
  re <- paste0(patron, "::(-?\\d+)")
  d <- tidy_coefs(m)[grepl(re, term)]
  d[, event_time := as.integer(sub(paste0(".*", re, ".*"), "\\1", term))]
  d <- d[, .(event_time, estimate = estimate * scale, se = se * scale)]
  d <- rbind(d, data.table(event_time = -1L, estimate = 0, se = 0))
  d[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]
  setorder(d, event_time)[]
}

# with ref = -1 dropped, every remaining negative event time is a pre-period
# coefficient, so the regex alone gives the joint pre-trend test.
#
# the omnibus test is reported together with the SLOPE of the pre-period path,
# because on its own it is misleading here: it tests against arbitrary
# alternatives and is close to blind to a smooth linear drift, which is exactly
# the shape that would matter. the slope is the exact linear contrast on the
# same coefficients, so both come from the model that the figure plots
test_pretend <- function(m, patron, scale) {
  w <- wald(m, paste0(patron, "::-"), print = FALSE)

  b <- coef(m); V <- vcov(m)
  nm <- grep(paste0("^", patron, "::-"), names(b), value = TRUE)
  tt <- as.integer(sub(paste0("^", patron, "::(-?\\d+).*"), "\\1", nm))
  o <- order(tt); nm <- nm[o]; tt <- tt[o]
  # tau = -1 is the omitted reference and enters the fit as a hard zero, so it
  # belongs in the centring even though it carries no coefficient
  wt <- c(tt, -1L) - mean(c(tt, -1L))
  wt <- (wt / sum(wt^2))[-(length(tt) + 1L)]
  pend <- sum(wt * b[nm]) * scale
  pend_se <- sqrt(as.numeric(t(wt) %*% V[nm, nm] %*% wt)) * scale

  data.table(pre_F = w$stat, pre_df = w$df1, pre_p = w$p,
             pend = pend, pend_se = pend_se,
             pend_p = 2 * pnorm(-abs(pend / pend_se)))
}

# latex cell "estimate^{stars} (se)"; * p<0,10 ** p<0,05 *** p<0,01
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

# ==============================================================================
# 1. estimation samples
# ==============================================================================

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g2_entry = as.IDate(g2_entry))]
panel[, year := year(ym)]

mi <- function(d) year(d) * 12L + (month(d) - 1L)

# brand x year absorbs chain-level pricing shocks: with three franchises holding
# ~80% of chilean stations, a national move by one of them is a confounder that
# station and region x year effects do not touch. distribuidor is a station-level
# constant by construction of the build, so it cannot be a bad control here
FE <- "station_key + ym + region^year + distribuidor^year"

armar <- function(roles) {
  d <- panel[role_entry %in% c("treated", roles)]
  # observations from a second nearby entry onward are dropped: the estimand is
  # the effect of a FIRST entry
  d <- d[is.na(g2_entry) | ym < g2_entry]
  d <- d[!(role_entry == "treated" & year(g_entry) < FOCAL_FROM)]
  # the control pool excludes post-2012 entrants: their own opening path is not
  # a valid counterfactual
  d <- d[role_entry == "treated" | base_2012 == TRUE]
  d[, treated := as.integer(role_entry == "treated")]
  d[, et := miym - mi(g_entry)]
  d[, rel := fifelse(treated == 1L,
                     pmax(-NBIN, pmin(NBIN, as.integer(floor(et / BIN_M)))),
                     -1L)]
  d[, post := fifelse(treated == 1L & !is.na(g_entry) & ym >= g_entry, 1L, 0L)]
  d[]
}

muestras <- lapply(grupos, armar)

# ==============================================================================
# 2. event study and static att, by outcome x control group x fuel
# ==============================================================================

estimar <- function(tipo, gr, fv) {
  cfg <- salidas[[tipo]]
  d <- muestras[[gr]][!is.na(get(fv))]
  lhs <- if (cfg$en_log) sprintf("log(%s)", fv) else fv
  scale <- if (cfg$en_log) 100 else 1

  m <- feols(as.formula(sprintf("%s ~ i(rel, treated, ref = -1) | %s", lhs, FE)),
             data = d, cluster = ~comuna)
  a <- coeftable(feols(as.formula(sprintf("%s ~ post | %s", lhs, FE)),
                       data = d, cluster = ~comuna))
  # for the price, also keep the effect in $/L: the percentage is comparable
  # across years but the peso figure is what a consumer actually sees
  a_lvl <- if (cfg$en_log) {
    coeftable(feols(as.formula(sprintf("%s ~ post | %s", fv, FE)),
                    data = d, cluster = ~comuna))
  } else NULL
  pt <- test_pretend(m, "rel", scale)

  tidy_es(m, "rel", scale = scale)[, `:=`(
    tipo = tipo, grupo = gr, outcome = fv,
    att = a["post", 1] * scale, att_se = a["post", 2] * scale,
    att_p = a["post", 4],
    att_lvl = if (is.null(a_lvl)) NA_real_ else a_lvl["post", 1],
    se_lvl  = if (is.null(a_lvl)) NA_real_ else a_lvl["post", 2],
    p_lvl   = if (is.null(a_lvl)) NA_real_ else a_lvl["post", 4],
    pre_F = pt$pre_F, pre_df = pt$pre_df, pre_p = pt$pre_p,
    pend = pt$pend, pend_se = pt$pend_se, pend_p = pt$pend_p,
    n_obs = nobs(m), n_trat = uniqueN(d[treated == 1L, station_key]),
    n_ctrl = uniqueN(d[treated == 0L, station_key]))]
}

res <- rbindlist(lapply(names(salidas), function(tipo) {
  rbindlist(lapply(names(grupos), function(gr) {
    rbindlist(lapply(names(salidas[[tipo]]$vars), estimar,
                     tipo = tipo, gr = gr))
  }))
}))

fwrite(res[, .(tipo, grupo, outcome, event_time, estimate, se, ci_low, ci_high)],
       file.path(TAB, "es_control_coefs.csv"))

att <- unique(res[, .(tipo, grupo, outcome, att, att_se, att_p,
                      att_lvl, se_lvl, p_lvl, n_obs, n_trat, n_ctrl)])
pretend <- unique(res[, .(tipo, grupo, outcome, pre_F, pre_df, pre_p,
                          pend, pend_se, pend_p)])
# 12 months of drift extrapolated forward, as a share of the estimated effect.
# this is the number that decides whether a pre-trend matters: a slope can be
# statistically indistinguishable from zero and still be large enough to
# account for most of the post-entry effect
pretend <- merge(pretend, att[, .(tipo, grupo, outcome, att)],
                 by = c("tipo", "grupo", "outcome"))
pretend[, `:=`(deriva_12m = pend * 12, share_att = pend * 12 / att)]
fwrite(att,     file.path(TAB, "es_control_att.csv"))
fwrite(pretend, file.path(TAB, "es_control_pretrend.csv"))

redondear <- function(d) d[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 3) else x)]

for (tp in names(salidas)) {
  cat("\n=== ATT DE LA ENTRADA -", toupper(tp), "(",
      if (tp == "precio") "% del precio" else "$/L", ") ===\n")
  print(redondear(dcast(att[tipo == tp], grupo ~ outcome, value.var = "att")))
  cat("--- pre-tendencias: p del test omnibus ---\n")
  print(redondear(dcast(pretend[tipo == tp], grupo ~ outcome,
                        value.var = "pre_p")))
  cat("--- pre-tendencias: p de la pendiente lineal (el test con potencia) ---\n")
  print(redondear(dcast(pretend[tipo == tp], grupo ~ outcome,
                        value.var = "pend_p")))
  cat("--- deriva de 12 meses como fraccion del att (|.| alto = problema) ---\n")
  print(redondear(dcast(pretend[tipo == tp], grupo ~ outcome,
                        value.var = "share_att")))
  cat("--- precision (ee del att) ---\n")
  print(redondear(dcast(att[tipo == tp], grupo ~ outcome,
                        value.var = "att_se")))
}

cat("\n=== TAMANO DE MUESTRA ===\n")
print(unique(att[tipo == "precio", .(grupo, n_trat = max(n_trat),
                                     n_ctrl = max(n_ctrl)), by = grupo][
  , .(grupo, n_trat, n_ctrl)]))

# ==============================================================================
# 3. balance: are the control groups comparable to the treated?
#
# measured over each station's UNTREATED months, and on the margin rather than
# the price, since the margin nets out the national mepco cost and is therefore
# comparable across calendar time
# ==============================================================================

sloc <- unique(panel[!is.na(lat), .(station_key, lat, lon)])
a1 <- sloc$lat * pi / 180; o1 <- sloc$lon * pi / 180
h <- sin(outer(a1, a1, "-") / 2)^2 +
  outer(cos(a1), cos(a1)) * sin(outer(o1, o1, "-") / 2)^2
h[h > 1] <- 1
D <- 2 * R_EARTH * asin(sqrt(h))
diag(D) <- Inf
ncomp <- data.table(station_key = sloc$station_key,
                    ncomp = rowSums(D <= RTREAT))

bal <- panel[is.na(g_entry) | ym < g_entry][
  , .(margen = mean(m93, na.rm = TRUE), precio = mean(p93, na.rm = TRUE),
      franquicia = mean(is_franchise, na.rm = TRUE),
      metropolitana = mean(region == "metropolitana"),
      role_entry = first(role_entry)),
  by = station_key]
bal <- merge(bal, ncomp, by = "station_key")

balance <- rbindlist(lapply(c("treated", names(grupos)), function(gr) {
  b <- bal[role_entry %in% (if (gr == "treated") "treated" else grupos[[gr]])]
  data.table(grupo = gr, n = nrow(b), margen = mean(b$margen, na.rm = TRUE),
             precio = mean(b$precio, na.rm = TRUE), ncomp = mean(b$ncomp),
             franquicia = mean(b$franquicia),
             metropolitana = mean(b$metropolitana))
}))
fwrite(balance, file.path(TAB, "es_control_balance.csv"))

cat("\n=== BALANCE: tratadas vs cada grupo de control (medias pre-tratamiento) ===\n")
print(balance[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ==============================================================================
# 4. latex tables
# ==============================================================================

tabla_att <- function(tp, value_est, value_se, value_p, dig, archivo) {
  w <- dcast(att[tipo == tp], outcome ~ grupo,
             value.var = c(value_est, value_se, value_p))
  out <- data.table(Combustible = salidas[[tp]]$vars[w$outcome])
  for (gr in names(grupos)) {
    out[[grupo_tex[gr]]] <- celda_tex(w[[paste0(value_est, "_", gr)]],
                                      w[[paste0(value_se, "_", gr)]],
                                      w[[paste0(value_p, "_", gr)]], dig)
  }
  guardar_tabla_tex(out, archivo)
}

tabla_att("precio", "att", "att_se", "att_p", 3, "tab_att_control.tex")
tabla_att("precio", "att_lvl", "se_lvl", "p_lvl", 2, "tab_att_control_pesos.tex")
tabla_att("margen", "att", "att_se", "att_p", 2, "tab_att_control_margen.tex")

tabla_pretend <- function(tp, archivo) {
  w <- dcast(pretend[tipo == tp], outcome ~ grupo,
             value.var = c("pre_p", "pend_p", "share_att"))
  out <- data.table(Combustible = salidas[[tp]]$vars[w$outcome])
  for (gr in names(grupos)) {
    out[[paste0(grupo_tex[gr], ": $p$ conj.")]] <-
      sprintf("$%.3f$", w[[paste0("pre_p_", gr)]])
    out[[paste0(grupo_tex[gr], ": $p$ pend.")]] <-
      sprintf("$%.3f$", w[[paste0("pend_p_", gr)]])
    out[[paste0(grupo_tex[gr], ": deriva/ATT")]] <-
      sprintf("$%.2f$", w[[paste0("share_att_", gr)]])
  }
  guardar_tabla_tex(out, archivo)
}

tabla_pretend("precio", "tab_pretend_control.tex")
tabla_pretend("margen", "tab_pretend_control_margen.tex")

tab_bal <- balance[, .(
  Grupo = c(treated = "Tratadas", anillo = grupo_tex[["anillo"]],
            lejano = grupo_tex[["lejano"]], nunca = grupo_tex[["nunca"]])[grupo],
  `$N$` = n,
  `Margen (\\$/L)` = sprintf("$%.1f$", margen),
  `Precio (\\$/L)` = sprintf("$%.0f$", precio),
  `Competidores $<$2 km` = sprintf("$%.2f$", ncomp),
  `Franquicia` = sprintf("$%.2f$", franquicia),
  `RM` = sprintf("$%.2f$", metropolitana))]
guardar_tabla_tex(tab_bal, "tab_balance_control.tex")

# ==============================================================================
# 5. figures
# ==============================================================================

AZUL_ES <- AZUL

# intervalos como barras de error y no como banda: con nueve bins discretos la
# banda sugiere una continuidad que el estimador no tiene
plot_grupo <- function(d, ylab, subt, color_var = NULL) {
  p <- if (is.null(color_var)) {
    ggplot(d, aes(event_time, estimate)) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.16,
                    colour = AZUL_ES) +
      geom_line(colour = AZUL_ES) +
      geom_point(colour = AZUL_ES, size = 1.4)
  } else {
    ggplot(d, aes(event_time, estimate, colour = .data[[color_var]])) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.16,
                    position = position_dodge(width = 0.3)) +
      geom_line(position = position_dodge(width = 0.3)) +
      geom_point(size = 1.3, position = position_dodge(width = 0.3))
  }
  p +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey55") +
    facet_wrap(~combustible, scales = "free_y") +
    scale_x_continuous(breaks = -NBIN:NBIN,
                       labels = c(sprintf("\u2264\u2212%d", NBIN),
                                  as.character(-(NBIN - 1):(NBIN - 1)),
                                  sprintf("\u2265%d", NBIN))) +
    labs(x = "Bins de seis meses desde la entrada (extremos agrupados)",
         y = ylab, colour = NULL, fill = NULL,
         title = "Efecto de la entrada de un competidor",
         subtitle = subt) +
    tema() +
    # marco por panel, igual que en 17_formato_fischer.R: cierra la grilla y
    # separa visualmente las facetas, que llevan escalas libres
    theme(panel.border = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.7))
}

for (tp in names(salidas)) {
  cfg <- salidas[[tp]]
  d <- res[tipo == tp]
  d[, combustible := factor(outcome, levels = names(cfg$vars),
                            labels = cfg$vars)]
  d[, grupo_f := factor(grupo_lbl[grupo], levels = grupo_lbl)]
  # 95 has no mepco reference, so the margin figures carry three panels
  alto <- if (tp == "precio") 5 else 3.4

  for (gr in names(grupos)) {
    guardar(plot_grupo(d[grupo == gr], cfg$ylab,
                      paste0("Control: estaciones ", tolower(grupo_lbl[gr]))),
           file.path(FIG, sprintf("es_control%s_%s.pdf", cfg$sufijo, gr)), 9, alto)
  }
  guardar(plot_grupo(d, cfg$ylab,
                     "Los tres grupos sobre la misma muestra de tratadas",
                     "grupo_f"),
          file.path(FIG, sprintf("es_control%s_comparacion.pdf", cfg$sufijo)),
          9.5, alto + 0.6)
}

# ==============================================================================
# 6. which control group to keep
#
# the criteria are robustness, not effect size: does the group survive the
# pre-trend test, how precise is it, how comparable is it to the treated, and
# is it contaminated by the treatment itself
# ==============================================================================

resumen <- merge(
  att[, .(att = mean(att), ee = mean(att_se), n_ctrl = max(n_ctrl)),
      by = .(tipo, grupo)],
  pretend[, .(pre_p_min = min(pre_p)), by = .(tipo, grupo)],
  by = c("tipo", "grupo"))
resumen <- merge(resumen, balance[, .(grupo, ncomp_ctrl = ncomp)], by = "grupo")
resumen[, dif_ncomp := ncomp_ctrl - balance[grupo == "treated", ncomp]]
setorder(resumen, tipo, grupo)

cat("\n=== RESUMEN PARA ELEGIR ===\n")
print(resumen[, .(tipo, grupo, att_medio = round(att, 3),
                  ee_medio = round(ee, 3), pre_p_min = round(pre_p_min, 3),
                  n_ctrl, dif_ncomp = round(dif_ncomp, 2))])

message("13_grupos_control.R: figuras en ", FIG, " y tablas en ", TAB, ".")
