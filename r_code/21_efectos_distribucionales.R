
# objective: the distributional effect of entry. the average effect answers how
# much prices fall; it does not answer WHERE in the price distribution they
# fall, and that is what determines who gains. if entry compresses the market
# from below, the gain concentrates on consumers who compare prices before
# buying; if it shifts the whole distribution, everyone gains alike.
#
# the market-level evidence in 16_mercado_max_media_min.R already points that
# way -- the market minimum falls several times more than the maximum -- but the
# max and the min are two order statistics, not the distribution. this script
# estimates the quantile treatment effect on the treated over the whole range.
#
# two estimators, deliberately:
#
#   MAIN. callaway & li (2019), via qte::panel_qtt, with the copula base
#   anchored at treatment onset (pre_copula = "short"). the default anchors the
#   copula window to the event horizon, which for the longer horizons of this
#   design falls outside the estimation window and leaves cells without a base. identifies the qtt under a
#   DISTRIBUTIONAL parallel-trends assumption -- the extension to the whole
#   distribution of the mean assumption this thesis already defends -- plus a
#   copula stability assumption: the dependence between the change in untreated
#   outcomes and their initial level is constant over time. it returns uniform
#   confidence bands over the quantile grid, which is what allows a statement
#   about the whole distribution rather than quantile by quantile.
#
#   ROBUSTNESS. change in changes of athey & imbens (2006), via qte::cic. it
#   needs no copula stability assumption, at the price of assuming a monotone
#   scalar-unobservable model. the two rest on different assumptions, so
#   agreement between them is informative: it makes the copula assumption a
#   weaker point of attack.
#
# fischer, martin & schmidt-dengler (2025) reach the same object through rif
# regressions plus the distribution regression of chernozhukov, fernandez-val &
# melly. the route taken here is better matched to this design, because its
# identifying assumption is the distributional version of the one the rest of
# the thesis defends rather than a separate modelling assumption.
#
# DEPENDENCIA DE VERSIONES. qte 2.0.0 llama a BMisc::combine_ecdfs con el
# argumento `ecdflist`, mientras que la version de CRAN de BMisc (1.4.8) lo
# nombra `dflist`. hay que instalar BMisc desde github para que la firma
# coincida:
#     remotes::install_github("bcallaway11/BMisc")
# si se reinstalan paquetes desde CRAN, este script vuelve a fallar con
# "argument dflist is missing, with no default".
#
# TAKES:  data/processed/panel_mensual.csv
# PRODUCES:
#   results/tables/qtt_coefs.csv, qtt_att.csv, tab_qtt.tex
#   results/figures/qtt_distribucional.pdf

library(data.table)
library(qte)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

FOCAL_FROM <- 2014L
# cohorts are ANNUAL rather than semestral. the balance requirement is what
# binds here: with semestral cohorts over 2016-2023 the design ends up with 196
# treated stations spread over fourteen cohorts, i.e. fourteen units per
# (g,t) cell, which is far too few to estimate an empirical cdf and invert
# seventeen quantiles off it. widening the cohort to a year nearly doubles the
# units per cell without losing treated stations
MESES_PER  <- 12L
ANIO_INI   <- 2015L   # window: balance is required over every period in it,
ANIO_FIN   <- 2024L   # so it trades sample size against horizon
PROBS  <- seq(0.10, 0.90, by = 0.05)
BITERS <- 200L
# the bootstrap recomputes every (g,t) cell on each iteration -- with fourteen
# cohorts that is around 196 cells per estimation, four estimations, and 200
# iterations, which is heavy enough to be worth parallelising. one core is left
# free for the rest of the machine
NCORES <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

fuels <- c(p93 = "Gasolina 93", pdi = "Diésel")

mi <- function(d) year(d) * 12L + (month(d) - 1L)

# vectorised: it is called on whole columns
celda_tex <- function(est, li, ls, dig = 3) {
  sig <- fifelse(is.na(li) | is.na(ls), "",
                 fifelse(li > 0 | ls < 0, "^{*}", ""))
  sprintf("$%.*f%s$ $[%.*f;\\,%.*f]$", dig, est, sig, dig, li, dig, ls)
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
# 1. semester panel, balanced
#
# both estimators need a balanced panel: callaway & li recover the missing
# copula from repeated observations of the same unit, which is only possible if
# the unit is observed throughout. the window and the period width are both
# choices with the same trade-off -- a longer window or a narrower period buys
# horizon and resolution but costs units per cell, and it is units per cell that
# distributional inference needs
# ==============================================================================

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g2_entry = as.IDate(g2_entry))]

d <- panel[role_entry %in% c("treated", "ring", "far")]
d <- d[is.na(g2_entry) | ym < g2_entry]
d <- d[!(role_entry == "treated" & year(g_entry) < FOCAL_FROM)]
d <- d[role_entry == "treated" | base_2012 == TRUE]
d[, sem := miym %/% MESES_PER]
d[, trat := as.integer(role_entry == "treated")]
# gname: first treated period, zero for the never treated, as ptetools expects
d[, g := fifelse(trat == 1L, as.integer(mi(g_entry) %/% MESES_PER), 0L)]

PER_ANIO <- 12L %/% MESES_PER          # periods per calendar year
LO <- ANIO_INI * PER_ANIO
HI <- ANIO_FIN * PER_ANIO + (PER_ANIO - 1L)
NSEM <- HI - LO + 1L

armar <- function(fv) {
  x <- d[!is.na(get(fv)) & get(fv) > 0 & sem >= LO & sem <= HI]
  x <- x[, .(y = mean(log(get(fv))) * 100), by = .(station_key, sem, g)]
  # balance: keep only stations observed in every semester of the window
  completos <- x[, .N, by = station_key][N == NSEM, station_key]
  x <- x[station_key %in% completos]
  # a cohort must fall inside the window and leave two pre periods for the
  # copula; earlier-treated stations are dropped rather than used as controls
  x[, g := fifelse(g > 0L & (g < LO + 2L | g > HI), NA_integer_, g)]
  x <- x[!is.na(g)]
  # ptetools expects numeric period and group columns, and is easier to reason
  # about when periods are consecutive integers starting at one
  per <- data.table(sem = sort(unique(x$sem)))
  per[, t := seq_len(.N)]
  x <- merge(x, per, by = "sem")
  x[, gt := fifelse(g == 0L, 0, as.numeric(per$t[match(g, per$sem)]))]
  x[, id := .GRP, by = station_key]
  x[, `:=`(t = as.numeric(t), gt = as.numeric(gt), id = as.numeric(id))]
  setorder(x, id, t)
  x[, .(id, t, gt, y, station_key)]
}

# ==============================================================================
# 2. estimation
# ==============================================================================

corre <- function(fv) {
  x <- as.data.frame(armar(fv))
  n_tr <- length(unique(x$id[x$gt > 0]))
  n_nt <- length(unique(x$id[x$gt == 0]))
  cat(sprintf("\n%s: %d estaciones (%d tratadas, %d nunca tratadas), %d cohortes\n",
              fv, length(unique(x$id)), n_tr, n_nt,
              length(unique(x$gt[x$gt > 0]))))

  # the object carries both pointwise and uniform bands; the uniform ones are
  # what licence a statement about the whole distribution rather than about one
  # quantile at a time
  saca <- function(obj, etiqueta) {
    q <- as.data.table(obj$overall)
    data.table(outcome = fv, estimador = etiqueta, prob = q$probs,
               qtt = q$qtt, se = q$se,
               li_pw = q$lower_pw, ls_pw = q$upper_pw,
               li = q$lower_ub, ls = q$upper_ub)
  }

  set.seed(20260731)
  m_cl <- panel_qtt(yname = "y", gname = "gt", tname = "t", idname = "id",
                    data = x, control_group = "nevertreated",
                    gt_type = "qtt", probs = PROBS, cband = TRUE,
                    biters = BITERS, pre_copula = "short", cl = NCORES)
  set.seed(20260731)
  m_ci <- cic(yname = "y", gname = "gt", tname = "t", idname = "id",
              data = x, panel = TRUE, control_group = "nevertreated",
              gt_type = "qtt", probs = PROBS, cband = TRUE, biters = BITERS,
              cl = NCORES)

  list(q = rbind(saca(m_cl, "Callaway-Li"), saca(m_ci, "Change-in-Changes")),
       n = data.table(outcome = fv, n_trat = n_tr, n_nunca = n_nt))
}

res <- lapply(names(fuels), corre)
qtt <- rbindlist(lapply(res, `[[`, "q"))
tam <- rbindlist(lapply(res, `[[`, "n"))

fwrite(qtt, file.path(TAB, "qtt_coefs.csv"))
fwrite(tam, file.path(TAB, "qtt_muestra.csv"))

cat("\n=== TAMANO DE MUESTRA ===\n"); print(tam)

cat("\n=== QTT POR CUANTIL (%) ===\n")
print(dcast(qtt[prob %in% c(0.1, 0.25, 0.5, 0.75, 0.9)],
            prob ~ outcome + estimador, value.var = "qtt")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

# whether the uniform band lies entirely below zero at every quantile is the
# statement about the whole distribution, and the reason for asking for cband
cat("\n=== SIGNIFICANCIA POR CUANTIL ===\n")
cat("pw = banda puntual bajo cero; ub = banda uniforme bajo cero\n")
print(qtt[, .(cuantiles = .N, pw_bajo_cero = sum(ls_pw < 0),
              ub_bajo_cero = sum(ls < 0),
              qtt_min = round(min(qtt), 3), qtt_max = round(max(qtt), 3),
              se_mediana = round(median(se), 3)),
          by = .(outcome, estimador)])

# ==============================================================================
# 3. table and figure
# ==============================================================================

sel <- qtt[prob %in% c(0.10, 0.25, 0.50, 0.75, 0.90)]
w <- dcast(sel, prob ~ outcome + estimador,
           value.var = c("qtt", "li_pw", "ls_pw"))
out <- data.table(Cuantil = sprintf("$q_{%.0f}$", 100 * w$prob))
for (fv in names(fuels)) for (es in c("Callaway-Li", "Change-in-Changes")) {
  k <- paste0(fv, "_", es)
  # the table reports pointwise bands, which are the readable ones per cell;
  # the uniform bands are in the figure and in the csv
  out[[sprintf("%s: %s", fuels[[fv]], es)]] <-
    celda_tex(w[[paste0("qtt_", k)]], w[[paste0("li_pw_", k)]],
              w[[paste0("ls_pw_", k)]])
}
guardar_tabla_tex(out, "tab_qtt.tex")

qtt[, combustible := factor(outcome, levels = names(fuels), labels = fuels)]
p <- ggplot(qtt, aes(prob, qtt, colour = estimador, fill = estimador)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_ribbon(aes(ymin = li_pw, ymax = ls_pw), alpha = 0.18, colour = NA) +
  geom_ribbon(aes(ymin = li, ymax = ls), alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  facet_wrap(~combustible, scales = "free_y") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("Callaway-Li" = AZUL,
                                 "Change-in-Changes" = NARANJO_OSC)) +
  scale_fill_manual(values = c("Callaway-Li" = AZUL,
                               "Change-in-Changes" = NARANJO_OSC)) +
  labs(x = "cuantil de la distribución de precios", y = "efecto sobre el precio (%)",
       colour = NULL, fill = NULL,
       title = "efecto de la entrada a lo largo de la distribución de precios",
       subtitle = paste("banda interior puntual, banda exterior uniforme;",
                        "control nunca tratadas")) +
  tema()
guardar(p, file.path(FIG, "qtt_distribucional.pdf"), 9.5, 5)

message("21_efectos_distribucionales.R: figuras en ", FIG, " y tablas en ", TAB, ".")
