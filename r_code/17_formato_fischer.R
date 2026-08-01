
# objective: the main specification of this project, presented in the format of
# figure 3 of fischer, martin & schmidt-dengler (2025): station-level event
# study with six-month bins, four bins before and five after, endpoints binned,
# with prices in levels in one panel and logged margins in the other, and the
# two market definitions as the two series.
#
# what is kept from the main specification of this project:
#   - never-treated stations as the control pool (13_grupos_control.R)
#   - station + month + region x year + brand x year fixed effects
#   - only the FIRST entry per incumbent; observations from a second one are
#     dropped
#   - standard errors clustered at comuna level
#
# what follows fischer et al. instead:
#   - treatment radii of 1 km and 2 km rather than the 2 km used elsewhere
#     here, so these numbers are not interchangeable with the other tables
#   - six-month bins over a +-4 window instead of monthly leads and lags
#   - prices in LEVELS ($/L) and margins in LOGS (effect read as % of the
#     margin), which is the opposite transform to the rest of this project and
#     is what makes the panel comparable to their figure
#
# their county-level time-varying controls have no counterpart here: chile
# publishes nothing at comuna level at monthly frequency that matches. region x
# year absorbs the slow-moving part of it.
#
# TAKES:  data/processed/panel_mensual.csv, entradas.csv
# PRODUCES:
#   results/tables/formato_fischer_coefs.csv, formato_fischer_att.csv,
#     tab_formato_fischer.tex
#   results/figures/formato_fischer.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

RADIOS <- c(1, 2)
BIN_M  <- 6L
NBIN   <- 4L
FOCAL_FROM <- 2014L
R_EARTH <- 6371

FE <- "station_key + ym_f + region^year + distribuidor^year"

fuels <- list(p93 = list(precio = "p93", margen = "m93", lbl = "Gasolina 93"),
              p97 = list(precio = "p97", margen = "m97", lbl = "Gasolina 97"),
              pdi = list(precio = "pdi", margen = "mdi", lbl = "Diésel"))

mi <- function(d) year(d) * 12L + (month(d) - 1L)

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
# 1. treatment at each of fischer's radii
# ==============================================================================

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), year = year(as.IDate(ym)))]

entradas <- fread("data/processed/entradas.csv")
entradas[, g := as.IDate(g)]

sloc <- unique(panel[!is.na(lat), .(station_key, lat, lon)])
sloc <- sloc[, .SD[1], by = station_key]
span <- panel[, .(first_mi = min(miym), last_mi = max(miym)), by = station_key]
sloc <- merge(sloc, span, by = "station_key")

a1 <- sloc$lat * pi / 180; o1 <- sloc$lon * pi / 180
a2 <- entradas$elat * pi / 180; o2 <- entradas$elon * pi / 180
h <- sin(outer(a1, a2, "-") / 2)^2 +
  outer(cos(a1), cos(a2)) * sin(outer(o1, o2, "-") / 2)^2
h[h > 1] <- 1
D_ev <- 2 * R_EARTH * asin(sqrt(h))

# never treated by its own entry, nor by one that predates its own arrival
own <- match(entradas$station_key, sloc$station_key)
D_ev[cbind(own, seq_len(nrow(entradas)))] <- Inf
inc <- outer(sloc$first_mi, mi(entradas$g), "<") &
  outer(sloc$last_mi, mi(entradas$g), ">=")
D_ev[!inc] <- Inf
g_ev <- mi(entradas$g)

armar <- function(r) {
  near <- lapply(seq_len(nrow(sloc)), function(i) {
    w <- which(D_ev[i, ] <= r)
    if (!length(w)) NULL else sort(unique(g_ev[w]))
  })
  asg <- data.table(
    station_key = sloc$station_key,
    g  = vapply(near, function(x) if (is.null(x)) NA_integer_ else x[1], integer(1)),
    g2 = vapply(near, function(x) if (length(x) >= 2L) x[2] else NA_integer_,
                integer(1)))

  d <- merge(panel, asg, by = "station_key")
  d[, treated := as.integer(!is.na(g))]
  d <- d[is.na(g2) | miym < g2]
  d <- d[!(treated == 1L & g %/% 12L < FOCAL_FROM)]
  # never-treated control pool: no entry within r km at any point. entrants are
  # excluded from it, their own opening path being no counterfactual
  d <- d[treated == 1L | base_2012 == TRUE]

  d[, et := miym - g]
  d[, bin := fifelse(treated == 1L,
                     pmax(-NBIN, pmin(NBIN, as.integer(floor(et / BIN_M)))),
                     -1L)]
  d[, post := fifelse(treated == 1L & !is.na(et) & et >= 0L, 1L, 0L)]
  d[, ym_f := miym]
  d[, radio := sprintf("Entrada en [0, %d] km", r)]
  d[]
}

muestras <- lapply(RADIOS, armar)
names(muestras) <- sprintf("Entrada en [0, %d] km", RADIOS)

cat("=== MUESTRAS ===\n")
for (nm in names(muestras)) {
  d <- muestras[[nm]]
  cat(nm, ": tratadas", uniqueN(d[treated == 1L, station_key]),
      "| control", uniqueN(d[treated == 0L, station_key]), "\n")
}

# ==============================================================================
# 2. estimation
#
# panel (a) is the price in LEVELS; panel (b) is log(margin) x 100, so the
# coefficient reads as a percentage of the margin. margins at or below zero
# cannot be logged and are dropped -- they are 0.1% of observations, and a
# negative gross margin against the national mepco quote is almost surely a
# reporting or matching error rather than a station selling below cost
# ==============================================================================

estimar <- function(nm, y, tipo) {
  d <- copy(muestras[[nm]])
  if (tipo == "margen") {
    n0 <- nrow(d[!is.na(get(y))])
    d <- d[!is.na(get(y)) & get(y) > 0]
    d[, yv := log(get(y)) * 100]
    perdidas <- n0 - nrow(d)
  } else {
    d <- d[!is.na(get(y))]
    d[, yv := get(y)]
    perdidas <- 0L
  }

  m <- feols(as.formula(sprintf("yv ~ i(bin, treated, ref = -1) | %s", FE)),
             data = d, cluster = ~comuna)
  ct <- coeftable(m)
  es <- data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2])
  es <- es[grepl("^bin::", term)]
  es[, bin := as.integer(sub("^bin::(-?[0-9]+).*", "\\1", term))]
  es <- rbind(es[, .(bin, estimate, se)],
              data.table(bin = -1L, estimate = 0, se = 0))

  a <- coeftable(feols(as.formula(sprintf("yv ~ post | %s", FE)),
                       data = d, cluster = ~comuna))
  # joint test on the pre-period bins, the same check used elsewhere here
  w <- wald(m, "bin::-", print = FALSE)

  es[, `:=`(radio = nm, outcome = y, tipo = tipo,
            ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se,
            att = a["post", 1], att_se = a["post", 2], att_p = a["post", 4],
            pre_F = w$stat, pre_p = w$p, n_obs = nobs(m),
            n_trat = uniqueN(d[treated == 1L, station_key]),
            n_drop = perdidas)]
  setorder(es, bin)[]
}

res <- rbindlist(lapply(names(muestras), function(nm)
  rbindlist(lapply(names(fuels), function(fk) rbind(
    estimar(nm, fuels[[fk]]$precio, "precio"),
    estimar(nm, fuels[[fk]]$margen, "margen"))))))

fwrite(res[, .(radio, outcome, tipo, bin, estimate, se, ci_low, ci_high)],
       file.path(TAB, "formato_fischer_coefs.csv"))

att <- unique(res[, .(radio, outcome, tipo, att, att_se, att_p, pre_F, pre_p,
                      n_trat, n_drop)])
fwrite(att, file.path(TAB, "formato_fischer_att.csv"))

cat("\n=== ATT (precio en $/L, margen en % del margen) ===\n")
print(att[, .(radio, outcome, tipo, att = round(att, 3),
              se = round(att_se, 3), p = round(att_p, 4), n_trat)])
cat("\n=== TEST CONJUNTO DE PRE-TENDENCIAS SOBRE LOS BINS PREVIOS ===\n")
print(att[, .(radio, outcome, F = round(pre_F, 2), p = round(pre_p, 4))])

# ==============================================================================
# 3. table and figure
# ==============================================================================

for (fk in names(fuels)) {
  f <- fuels[[fk]]
  w <- dcast(att[outcome %in% c(f$precio, f$margen)], tipo ~ radio,
             value.var = c("att", "att_se", "att_p"))
  out <- data.table(Resultado = c(margen = "$\\ln$(margen), \\%",
                                  precio = "Precio, \\$/L")[w$tipo])
  for (nm in names(muestras)) {
    out[[nm]] <-
      celda_tex(w[[paste0("att_", nm)]], w[[paste0("att_se_", nm)]],
                w[[paste0("att_p_", nm)]])
  }
  guardar_tabla_tex(out, sprintf("tab_formato_fischer_%s.tex", fk))
}

etiquetas <- c(sprintf("\u2264\u2212%d", NBIN),
               as.character(-(NBIN - 1):(NBIN - 1)),
               sprintf("\u2265%d", NBIN))
pal <- c(AZUL, NARANJO)

panel_plot <- function(d) {
  ggplot(d, aes(bin, estimate, color = radio, shape = radio, linetype = radio)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = -0.5, linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15,
                  linewidth = 0.6, linetype = "solid",
                  position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.3), size = 2.1) +
    scale_x_continuous(breaks = -NBIN:NBIN, labels = etiquetas) +
    scale_color_manual(values = pal) +
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    # filas por transformacion y columnas por combustible: las unidades difieren
    # entre filas ($/L y %), de modo que la escala se libera por fila y se
    # comparte entre combustibles, que es la comparacion de interes
    facet_grid(tipo_f ~ combustible, scales = "free_y", switch = "y") +
    labs(x = "Bins de seis meses respecto de la entrada", y = NULL,
         color = NULL, shape = NULL, linetype = NULL,
         title = "Efecto de la entrada sobre el precio y el margen",
         subtitle = "Control: estaciones nunca tratadas") +
    tema() +
    theme(strip.placement = "outside")
}

lbl_comb <- vapply(fuels, `[[`, character(1), "lbl")
res[, combustible := factor(
  lbl_comb[match(outcome, unlist(lapply(fuels, function(f) c(f$precio, f$margen))) )],
  levels = lbl_comb)]
# el match anterior recorre precio y margen de cada combustible en orden
res[, combustible := factor(
  fifelse(outcome %in% c("p93", "m93"), "Gasolina 93",
          fifelse(outcome %in% c("p97", "m97"), "Gasolina 97", "Diésel")),
  levels = unname(lbl_comb))]
res[, tipo_f := factor(fifelse(tipo == "precio", "Precio ($/L)",
                               "ln(margen) (%)"),
                       levels = c("Precio ($/L)", "ln(margen) (%)"))]

guardar(panel_plot(res), file.path(FIG, "formato_fischer.pdf"), 11, 6)

message("17_formato_fischer.R: figuras en ", FIG, " y tablas en ", TAB, ".")
