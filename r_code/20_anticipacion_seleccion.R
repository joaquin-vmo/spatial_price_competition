
# objective: the two exercises that speak most directly to the maintained
# assumption -- that the TIMING of a rival's entry is exogenous to the
# incumbent, conditional on the fixed effects.
#
#   A. anticipation, at weekly resolution. building a service station is visible
#      for months before it opens, so the sharpest objection to the design is
#      that incumbents move first and the estimated jump is only the tail of a
#      response that started earlier. the six-month bins of the main
#      specification cannot see that. the weekly panel can: it puts about
#      twenty-six observations on each side of the opening week, where the main
#      figure has one bin. note that anticipation, if present, ATTENUATES the
#      estimate -- prices would already be lower at the reference period -- so
#      this exercise either closes the objection or turns the headline number
#      into a lower bound.
#
#   B. where entrants locate, in levels. the station fixed effects of the main
#      specification are what the design removes: each station's price level
#      net of time effects. regressing them on the distance to the nearest
#      entry recovers the selection the design neutralises, and characterises
#      it. this is the exercise of table 7 of \textcite{ARCIDIACONO2020}, who
#      use it to show that walmart locates next to low-priced supermarkets. it
#      does not test the timing assumption; it measures the LOCATION
#      endogeneity that the assumption explicitly allows, and makes visible how
#      much of the raw association the fixed effects are absorbing.
#
# TAKES:  data/processed/panel_semanal.csv.gz, panel_mensual.csv, entradas.csv
# PRODUCES:
#   results/tables/anticipacion_coefs.csv, anticipacion_test.csv
#   results/tables/seleccion_ubicacion.csv, tab_seleccion.tex
#   results/figures/anticipacion.pdf, seleccion_ubicacion.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

RTREAT <- 2
FOCAL_FROM <- 2014L
VENT_W  <- 26L      # weeks kept each side of the opening; endpoints binned
R_EARTH <- 6371

fuels <- c(p93 = "gasolina 93", pdi = "diesel")

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

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g2_entry = as.IDate(g2_entry))]
panel[, year := year(ym)]
entradas <- fread("data/processed/entradas.csv")[, g := as.IDate(g)]

# ==============================================================================
# A. ANTICIPATION AT WEEKLY RESOLUTION
# ==============================================================================

sem <- fread("data/processed/panel_semanal.csv.gz",
             select = c("station_key", "fuel", "wi", "ym", "price", "comuna",
                        "region", "distribuidor", "lat", "lon", "base_2012"))
sem[, ym := as.IDate(ym)]
sem[, year := year(ym)]

# the opening WEEK of each entrant, which the monthly build only knows to the
# month. weekly precision is the whole point of this exercise
primera <- sem[, .(wi_ini = min(wi)), by = station_key]
ent <- merge(entradas[, .(station_key, g, elat, elon)], primera, by = "station_key")
ent[, g_mi := mi(g)]

sloc <- unique(sem[!is.na(lat), .(station_key, lat, lon)])[, .SD[1], by = station_key]
a1 <- sloc$lat * pi / 180; o1 <- sloc$lon * pi / 180
a2 <- ent$elat * pi / 180; o2 <- ent$elon * pi / 180
h <- sin(outer(a1, a2, "-") / 2)^2 +
  outer(cos(a1), cos(a2)) * sin(outer(o1, o2, "-") / 2)^2
h[h > 1] <- 1
D_ev <- 2 * R_EARTH * asin(sqrt(h))

# for each treated focal, the opening week of the entry that treats it
trat <- unique(panel[role_entry == "treated" & !is.na(g_entry) &
                       year(g_entry) >= FOCAL_FROM,
                     .(station_key, g_mi = mi(g_entry))])
idx <- match(trat$station_key, sloc$station_key)
trat[, wi_evento := vapply(seq_len(.N), function(k) {
  i <- idx[k]
  if (is.na(i)) return(NA_integer_)
  w <- which(D_ev[i, ] <= RTREAT & ent$g_mi == trat$g_mi[k])
  if (!length(w)) NA_integer_ else min(ent$wi_ini[w])
}, integer(1))]
trat <- trat[!is.na(wi_evento)]
cat("tratadas con semana de evento identificada:", nrow(trat), "\n")

# second nearby entry, to cut the window as in the main specification
g2 <- unique(panel[role_entry == "treated" & !is.na(g2_entry),
                   .(station_key, g2_mi = mi(g2_entry))])

d <- merge(sem, trat[, .(station_key, wi_evento)], by = "station_key", all.x = TRUE)
d <- merge(d, g2, by = "station_key", all.x = TRUE)
d[, treated := as.integer(!is.na(wi_evento))]
# controls: never treated at RTREAT and present since 2012
ctrl_ok <- unique(panel[role_entry %in% c("ring", "far") & base_2012 == TRUE,
                        station_key])
d <- d[treated == 1L | station_key %in% ctrl_ok]
d <- d[is.na(g2_mi) | mi(ym) < g2_mi]
d[, etw := fifelse(treated == 1L,
                   pmax(-VENT_W, pmin(VENT_W, wi - wi_evento)), -1L)]

pend_lin <- function(m, sc) {
  b <- coef(m); V <- vcov(m)
  nm <- grep("^etw::-", names(b), value = TRUE)
  tt <- as.integer(sub("^etw::(-?[0-9]+).*", "\\1", nm))
  o <- order(tt); nm <- nm[o]; tt <- tt[o]
  w <- c(tt, -1L) - mean(c(tt, -1L)); w <- (w / sum(w^2))[-(length(tt) + 1L)]
  e <- sum(w * b[nm]) * sc
  s <- sqrt(as.numeric(t(w) %*% V[nm, nm] %*% w)) * sc
  data.table(pend = e, pend_se = s, pend_p = 2 * pnorm(-abs(e / s)))
}

FE_W <- "station_key + wi + region^year + distribuidor^year"

anti <- rbindlist(lapply(names(fuels), function(fv) {
  dd <- d[fuel == sub("^p", "", fv) | (fv == "pdi" & fuel == "di")]
  dd <- dd[!is.na(price) & price > 0]
  m <- feols(as.formula(sprintf("log(price) ~ i(etw, treated, ref = -1) | %s",
                                FE_W)), data = dd, cluster = ~comuna)
  ct <- coeftable(m)
  r <- data.table(term = rownames(ct), estimate = ct[, 1] * 100, se = ct[, 2] * 100)
  r <- r[grepl("^etw::", term)]
  r[, etw := as.integer(sub("^etw::(-?[0-9]+).*", "\\1", term))]
  r <- rbind(r[, .(etw, estimate, se)], data.table(etw = -1L, estimate = 0, se = 0))
  w <- wald(m, "etw::-", print = FALSE)
  pl <- pend_lin(m, 100)
  r[, `:=`(outcome = fv, ci_low = estimate - 1.96 * se,
           ci_high = estimate + 1.96 * se,
           pre_F = w$stat, pre_p = w$p, pend = pl$pend, pend_p = pl$pend_p,
           n_obs = nobs(m), n_trat = uniqueN(dd[treated == 1L, station_key]))]
  setorder(r, etw)[]
}))

fwrite(anti[, .(outcome, etw, estimate, se, ci_low, ci_high)],
       file.path(TAB, "anticipacion_coefs.csv"))
test_anti <- unique(anti[, .(outcome, pre_F, pre_p, pend, pend_p, n_trat)])
fwrite(test_anti, file.path(TAB, "anticipacion_test.csv"))

cat("\n=== ANTICIPACION: test sobre las 25 semanas previas ===\n")
print(test_anti[, .(outcome, F = round(pre_F, 2), p_conj = round(pre_p, 4),
                    pend_semanal = round(pend, 4), p_pend = round(pend_p, 4),
                    n_trat)])

anti[, combustible := factor(outcome, levels = names(fuels), labels = fuels)]
p_anti <- ggplot(anti, aes(etw, estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey55") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, fill = AZUL) +
  geom_line(colour = AZUL) +
  facet_wrap(~combustible, scales = "free_y") +
  labs(x = "semanas desde la apertura de la entrante (extremos agrupados)",
       y = "efecto sobre el precio (%)",
       title = "¿anticipan las incumbentes la entrada?",
       subtitle = paste("resolución semanal alrededor de la apertura;",
                        "control nunca tratadas")) +
  tema()
guardar(p_anti, file.path(FIG, "anticipacion.pdf"), 9.5, 4.6)

rm(sem, d); invisible(gc())

# ==============================================================================
# B. WHERE ENTRANTS LOCATE: THE STATION FIXED EFFECTS
# ==============================================================================

dm <- panel[role_entry %in% c("treated", "ring", "far")]
dm <- dm[is.na(g2_entry) | ym < g2_entry]
dm <- dm[!(role_entry == "treated" & year(g_entry) < FOCAL_FROM)]
dm <- dm[role_entry == "treated" | base_2012 == TRUE]
dm[, post := fifelse(role_entry == "treated" & !is.na(g_entry) & ym >= g_entry,
                     1L, 0L)]

# distance from every station to the nearest entry event, at any date
idx2 <- match(sloc$station_key, sloc$station_key)
dmin <- data.table(station_key = sloc$station_key,
                   dist_min = apply(D_ev, 1, min))
bandas <- c(0, 1, 2, 3, 5, 10, Inf)
lbl_b <- c("$<$1 km", "1--2 km", "2--3 km", "3--5 km", "5--10 km", "$>$10 km")
dmin[, banda := cut(dist_min, breaks = bandas, labels = lbl_b,
                    include.lowest = TRUE)]

sel <- rbindlist(lapply(names(fuels), function(fv) {
  dd <- dm[!is.na(get(fv))]
  m <- feols(as.formula(sprintf(
    "log(%s) ~ post | station_key + ym + region^year + distribuidor^year", fv)),
    data = dd, cluster = ~comuna)
  fe <- fixef(m)$station_key
  fx <- data.table(station_key = names(fe), ef = as.numeric(fe) * 100)
  fx <- merge(fx, dmin, by = "station_key")
  fx <- merge(fx, unique(dd[, .(station_key, region, comuna, distribuidor)]),
              by = "station_key")
  # the fixed effects are a price level net of time effects; regressing them on
  # the distance band recovers where entrants chose to locate
  mm <- feols(ef ~ i(banda, ref = "$>$10 km") | region + distribuidor,
              data = fx, cluster = ~comuna)
  ct <- coeftable(mm)
  r <- data.table(term = rownames(ct), est = ct[, 1], se = ct[, 2], p = ct[, 4])
  r <- r[grepl("banda::", term)]
  r[, banda := factor(sub(".*banda::", "", term), levels = lbl_b)]
  r[, `:=`(outcome = fv, n = nrow(fx))]
  r[, .(outcome, banda, est, se, p, n)]
}))
fwrite(sel, file.path(TAB, "seleccion_ubicacion.csv"))

cat("\n=== EFECTO FIJO DE ESTACION SEGUN DISTANCIA A LA ENTRADA MAS CERCANA ===\n")
cat("(% respecto de las estaciones a mas de 10 km de toda entrada)\n")
print(dcast(sel, banda ~ outcome, value.var = "est")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

w <- dcast(sel, banda ~ outcome, value.var = c("est", "se", "p"))
tab_sel <- data.table(`Distancia a la entrada más cercana` = as.character(w$banda))
for (fv in names(fuels)) {
  tab_sel[[unname(fuels[fv])]] <- celda_tex(w[[paste0("est_", fv)]],
                                            w[[paste0("se_", fv)]],
                                            w[[paste0("p_", fv)]])
}
guardar_tabla_tex(tab_sel, "tab_seleccion.tex")

sel[, combustible := factor(outcome, levels = names(fuels), labels = fuels)]
sel[, `:=`(ci_low = est - 1.96 * se, ci_high = est + 1.96 * se)]
p_sel <- ggplot(sel, aes(banda, est)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, colour = AZUL) +
  geom_point(colour = AZUL, size = 1.8) +
  facet_wrap(~combustible) +
  labs(x = "distancia a la entrada más cercana",
       y = "nivel de precio de la estación (%)",
       title = "dónde se instalan las entrantes",
       subtitle = paste("efectos fijos de estación de la especificación",
                        "principal; base: estaciones a más de 10 km")) +
  tema() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
guardar(p_sel, file.path(FIG, "seleccion_ubicacion.pdf"), 9.5, 4.6)

message("20_anticipacion_seleccion.R: figuras en ", FIG, " y tablas en ", TAB, ".")
