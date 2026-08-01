
# objective: replicate the market-level figure of fischer, martin &
# schmidt-dengler (2025): the effect of entry on the MAXIMUM, MEAN and MINIMUM
# price of the local market, in six-month bins, for two market definitions.
#
# what changes relative to every other script here: the unit of observation is
# the MARKET, not the station. a market is the circle of radius r around a focal
# incumbent, and the outcome is the max/mean/min price across all stations
# active inside it. the entrant's own price IS included once it opens -- fischer
# et al. exclude entrants from the station-level analysis but state in their
# footnote 16 that entrant prices "still matter for the consumer's choice and
# rents in the end, so that we will take them into account in our market-level
# analysis on market-level minimum, mean, and maximum prices".
#
# reading the figure: the mean is the average price effect already estimated
# elsewhere. the gap between the three series is the distributional content --
# if entry compresses the market from below, the minimum falls by more than the
# maximum.
#
# radii are fischer's 1 km and 2 km, not the 2 km treatment radius
# in the rest of this project, so the two sets of numbers are not interchangeable.
#
# TAKES:  data/processed/panel_mensual.csv, entradas.csv
# PRODUCES:
#   results/tables/mercado_coefs.csv, mercado_att.csv, tab_mercado.tex
#   results/figures/mercado_max_media_min_{p93,pdi}.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

RADIOS <- c(1, 2)          # km, as in fischer et al.
BIN_M  <- 6L               # months per bin
NBIN   <- 4L               # bins kept each side; endpoints are binned
FOCAL_FROM <- 2014L
MIN_FIRMS <- 2L            # a market needs two firms for a max/min to differ
R_EARTH <- 6371

fuels <- c(p93 = "gasolina 93", pdi = "diesel")
series <- c(max = "precio maximo", media = "precio medio", min = "precio minimo")

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
# 1. inputs
# ==============================================================================

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), year = year(as.IDate(ym)))]

entradas <- fread("data/processed/entradas.csv")
entradas[, g := as.IDate(g)]

sloc <- unique(panel[!is.na(lat), .(station_key, lat, lon, region, comuna,
                                    base_2012)])
sloc <- sloc[, .SD[1], by = station_key]

span <- panel[, .(first_mi = min(miym), last_mi = max(miym)), by = station_key]
sloc <- merge(sloc, span, by = "station_key")

dist_km <- function(la1, lo1, la2, lo2) {
  a1 <- la1 * pi / 180; o1 <- lo1 * pi / 180
  a2 <- la2 * pi / 180; o2 <- lo2 * pi / 180
  h <- sin(outer(a1, a2, "-") / 2)^2 +
    outer(cos(a1), cos(a2)) * sin(outer(o1, o2, "-") / 2)^2
  h[h > 1] <- 1
  2 * R_EARTH * asin(sqrt(h))
}

D_st <- dist_km(sloc$lat, sloc$lon, sloc$lat, sloc$lon)          # station-station
D_ev <- dist_km(sloc$lat, sloc$lon, entradas$elat, entradas$elon) # station-entry

# a station is never treated by its own entry, nor by one that predates it
own <- match(entradas$station_key, sloc$station_key)
D_ev[cbind(own, seq_len(nrow(entradas)))] <- Inf
inc <- outer(sloc$first_mi, mi(entradas$g), "<") &
  outer(sloc$last_mi, mi(entradas$g), ">=")
D_ev[!inc] <- Inf
g_ev <- mi(entradas$g)

# ==============================================================================
# 2. one estimation sample per radius
# ==============================================================================

armar_mercado <- function(r) {
  # --- market aggregates: every active station within r km of the focal one,
  # the focal included, and the entrant included from the month it opens
  ix <- which(D_st <= r, arr.ind = TRUE)          # diagonal kept: own market
  edges <- data.table(market = sloc$station_key[ix[, "row"]],
                      miembro = sloc$station_key[ix[, "col"]])

  px <- panel[, .(miembro = station_key, miym, p93, pdi)]
  mk <- merge(edges, px, by = "miembro", allow.cartesian = TRUE)
  mk <- mk[, .(n_firms = .N,
               max_p93 = max(p93, na.rm = TRUE), media_p93 = mean(p93, na.rm = TRUE),
               min_p93 = min(p93, na.rm = TRUE),
               max_pdi = max(pdi, na.rm = TRUE), media_pdi = mean(pdi, na.rm = TRUE),
               min_pdi = min(pdi, na.rm = TRUE)),
           by = .(market, miym)]
  for (cl in names(mk)) set(mk, which(is.infinite(mk[[cl]])), cl, NA_real_)

  # --- treatment: first and second entry within r km while already incumbent
  near <- lapply(seq_len(nrow(sloc)), function(i) {
    w <- which(D_ev[i, ] <= r)
    if (!length(w)) NULL else sort(unique(g_ev[w]))
  })
  asg <- data.table(
    market = sloc$station_key,
    g  = vapply(near, function(x) if (is.null(x)) NA_integer_ else x[1], integer(1)),
    g2 = vapply(near, function(x) if (length(x) >= 2L) x[2] else NA_integer_,
                integer(1)))

  d <- merge(mk, asg, by = "market")
  d <- merge(d, sloc[, .(market = station_key, region, comuna, base_2012)],
             by = "market")
  d[, year := miym %/% 12L]

  d[, treated := as.integer(!is.na(g))]
  # drop from the second entry on: the estimand is the effect of a FIRST entry
  d <- d[is.na(g2) | miym < g2]
  # a suspect 2013 first entry neither treats nor serves as control
  d <- d[!(treated == 1L & g %/% 12L < FOCAL_FROM)]
  # control markets are centred on stations present since 2012, so an entrant's
  # own market is never used as a counterfactual
  d <- d[treated == 1L | base_2012 == TRUE]

  # --- at least two firms, measured BEFORE treatment so the restriction cannot
  # select on the entry itself
  pre_n <- d[treated == 0L | miym < g, .(n_pre = mean(n_firms)), by = market]
  d <- merge(d, pre_n, by = "market")
  d <- d[n_pre >= MIN_FIRMS]

  # --- six-month bins, endpoints binned
  d[, et := miym - g]
  d[, bin := fifelse(treated == 1L,
                     pmax(-NBIN, pmin(NBIN, as.integer(floor(et / BIN_M)))),
                     -1L)]
  d[, post := fifelse(treated == 1L & !is.na(et) & et >= 0L, 1L, 0L)]
  d[, radio := r]
  d[]
}

muestras <- lapply(RADIOS, armar_mercado)
names(muestras) <- paste0(RADIOS, " km")

cat("=== MUESTRAS DE MERCADO ===\n")
for (nm in names(muestras)) {
  d <- muestras[[nm]]
  cat(nm, ": mercados", uniqueN(d$market),
      "| tratados", uniqueN(d[treated == 1L, market]),
      "| control", uniqueN(d[treated == 0L, market]),
      "| firmas pre (mediana)", round(median(unique(d[, .(market, n_pre)])$n_pre), 2),
      "\n")
}

# ==============================================================================
# 3. estimation
#
# market fixed effect plus month and region x year, mirroring the station-level
# design. brand x year is dropped here on purpose: the outcome belongs to the
# market, not to the focal station, so its brand is not the right control
# ==============================================================================

FE <- "market + ym_f + region^year"

estimar <- function(nm, fv, s) {
  d <- muestras[[nm]]
  y <- paste0(s, "_", fv)
  dd <- d[!is.na(get(y))]
  dd[, ym_f := miym]

  m <- feols(as.formula(sprintf("%s ~ i(bin, treated, ref = -1) | %s", y, FE)),
             data = dd, cluster = ~comuna)
  ct <- coeftable(m)
  es <- data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2])
  es <- es[grepl("^bin::", term)]
  es[, bin := as.integer(sub("^bin::(-?[0-9]+).*", "\\1", term))]
  es <- rbind(es[, .(bin, estimate, se)],
              data.table(bin = -1L, estimate = 0, se = 0))
  a <- coeftable(feols(as.formula(sprintf("%s ~ treated:post | %s", y, FE)),
                       data = dd, cluster = ~comuna))
  es[, `:=`(radio = nm, fuel = fv, serie = s,
            ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se,
            att = a[1, 1], att_se = a[1, 2], att_p = a[1, 4],
            n_obs = nobs(m), n_mercados = uniqueN(dd$market))]
  setorder(es, bin)[]
}

res <- rbindlist(lapply(names(muestras), function(nm)
  rbindlist(lapply(names(fuels), function(fv)
    rbindlist(lapply(names(series), function(s) estimar(nm, fv, s)))))))

fwrite(res[, .(radio, fuel, serie, bin, estimate, se, ci_low, ci_high)],
       file.path(TAB, "mercado_coefs.csv"))

att <- unique(res[, .(radio, fuel, serie, att, att_se, att_p, n_mercados)])
fwrite(att, file.path(TAB, "mercado_att.csv"))

cat("\n=== ATT SOBRE EL PRECIO DE MERCADO ($/L) ===\n")
print(dcast(att, radio + serie ~ fuel, value.var = "att")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
cat("\n--- ee ---\n")
print(dcast(att, radio + serie ~ fuel, value.var = "att_se")[
  , lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

# ==============================================================================
# 4. table and figures
# ==============================================================================

for (fv in names(fuels)) {
  w <- dcast(att[fuel == fv], serie ~ radio,
             value.var = c("att", "att_se", "att_p"))
  out <- data.table(Serie = series[w$serie])
  for (nm in names(muestras)) {
    out[[paste0("Entrada en [0,", sub(" km", "", nm), "] km")]] <-
      celda_tex(w[[paste0("att_", nm)]], w[[paste0("att_se_", nm)]],
                w[[paste0("att_p_", nm)]])
  }
  guardar_tabla_tex(out, sprintf("tab_mercado_%s.tex", fv))
}

etiquetas <- c(sprintf("<=-%d", NBIN), as.character(-(NBIN - 1):(NBIN - 1)),
               sprintf(">=%d", NBIN))
pal <- c(max = AZUL, media = NARANJO_OSC, min = NARANJO)
lty <- c(max = "solid", media = "dashed", min = "dashed")
shp <- c(max = 16, media = 17, min = 15)

for (fv in names(fuels)) {
  d <- res[fuel == fv]
  d[, `:=`(serie_f = factor(series[serie], levels = series),
           panel = factor(sprintf("entrada en [0,%s] km", sub(" km", "", radio))))]
  p <- ggplot(d, aes(bin, estimate, color = serie_f, shape = serie_f,
                     linetype = serie_f)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = -0.5, linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18,
                  linetype = "solid",
                  position = position_dodge(width = 0.35)) +
    geom_line(position = position_dodge(width = 0.35), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.35), size = 1.7) +
    facet_wrap(~panel) +
    scale_x_continuous(breaks = -NBIN:NBIN, labels = etiquetas) +
    scale_color_manual(values = unname(pal[names(series)])) +
    scale_shape_manual(values = unname(shp[names(series)])) +
    scale_linetype_manual(values = unname(lty[names(series)])) +
    labs(x = "bins de seis meses respecto de la entrada",
         y = "efecto sobre el precio ($/L)",
         color = NULL, shape = NULL, linetype = NULL,
         title = sprintf("efecto de la entrada sobre el precio de mercado: %s",
                         fuels[[fv]]),
         subtitle = paste("mercado = circulo alrededor de la incumbente focal,",
                          "incluye al entrante; ee cluster comuna")) +
    tema()
  guardar(p, file.path(FIG, sprintf("mercado_max_media_min_%s.pdf", fv)),
          9.5, 5)
}

message("16_mercado_max_media_min.R: figuras en ", FIG, " y tablas en ", TAB, ".")
