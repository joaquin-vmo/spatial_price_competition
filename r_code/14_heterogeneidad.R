
# objective: does the entry effect depend on how much margin the incumbent was
# already taking, relative to its own local market?
#
# the moderator is the station's price rank inside its comuna. within a
# comuna-month the mepco wholesale cost is identical for every station, so the
# price rank IS the margin rank -- ranking the price is what the code does,
# which keeps the months before mepco starts (2014-08) usable.
#
# the split is binary: the top tercile of its comuna against the rest. a
# station pricing at the top has room to cut that a station already at the
# bottom does not, and the two lower terciles behave alike, so pooling them
# gives a cleaner contrast than three curves on top of each other.
#
# two design points that decide whether this is credible:
#
#   the moderator is measured over the station's PRE-treatment months, never
#   contemporaneously. a margin rank measured after entry is an outcome of the
#   entry, and interacting the treatment with it would be a bad control.
#
#   selecting stations on a high margin and following them afterwards produces
#   a decline through mean reversion alone. two things guard against reading
#   that as an entry effect: the moderator is a pre-window AVERAGE rather than
#   a single month, and the same cutoff is applied to control stations, so
#   whatever reversion exists is present in both arms and differences out.
#   the per-group pre-trend slope reported below is the check on this.
#
# TAKES:  data/processed/panel_mensual.csv  (from 10_build_event_panel.R)
# PRODUCES:
#   results/tables/het_att.csv, het_diferencia.csv, het_pretendencia.csv,
#     het_composicion.csv, het_robustez_ncomp.csv, het_robustez_minn.csv,
#     tab_het_rank.tex
#   results/figures/het_rank_precio.pdf, het_rank_margen.pdf

library(data.table)
library(fixest)
library(ggplot2)
source("r_code/00_estilo.R")

FIG <- "results/figures"
TAB <- "results/tables"

PRE <- 12L
POST <- 12L
FOCAL_FROM <- 2014L
RTREAT <- 2       # treatment radius (km), same as the main specification
MIN_N <- 4L      # min stations in a comuna-month for its rank to mean anything
CORTE <- 2 / 3   # top tercile of the treated distribution is the "high" group

precios  <- c(p93 = "Gasolina 93", p95 = "Gasolina 95",
              p97 = "Gasolina 97", pdi = "Diésel")
margenes <- c(m93 = "Gasolina 93", m97 = "Gasolina 97", mdi = "Diésel")

# brand x year is in the main specification: with three franchises holding ~80%
# of chilean stations, a chain-level national move is a confounder that station
# and region x year effects do not absorb
FE <- "station_key + ym + region^year + distribuidor^year"

# etiquetas cortas para la tabla y la leyenda; el texto explica que el tercil se
# define dentro de la comuna, de modo que repetirlo en cada fila solo desborda
lbl <- c("Margen bajo o medio", "Margen alto")

tidy_coefs <- function(m) {
  ct <- coeftable(m)
  data.table(term = rownames(ct), estimate = ct[, 1], se = ct[, 2], p = ct[, 4])
}

celda_tex <- function(est, se, p, dig = 3) {
  st <- fifelse(p < 0.01, "^{***}",
                fifelse(p < 0.05, "^{**}", fifelse(p < 0.10, "^{*}", "")))
  sprintf("$%.*f%s$ ($%.*f$)", dig, est, st, dig, se)
}

# la envoltura en \small es necesaria: con cinco columnas de estimador y error
# estandar la tabla no cabe a tamano normal en el ancho de texto
guardar_tabla_tex <- function(d, archivo, regla_antes_de = NULL) {
  cuerpo <- apply(d, 1, function(r) paste(paste(r, collapse = " & "), "\\\\"))
  if (!is.null(regla_antes_de)) {
    k <- which(d[[1]] == regla_antes_de)
    if (length(k)) cuerpo[k] <- paste("\\midrule", cuerpo[k])
  }
  writeLines(
    c("\\small", sprintf("\\begin{tabular}{l%s}", strrep("c", ncol(d) - 1L)),
      "\\toprule",
      paste(paste(names(d), collapse = " & "), "\\\\"), "\\midrule",
      cuerpo, "\\bottomrule", "\\end{tabular}"),
    file.path(TAB, archivo)
  )
}

# ==============================================================================
# 1. estimation sample: same as the chosen design in 13_grupos_control.R
#    (never-treated control pool, first entry only, entrants excluded)
# ==============================================================================

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g2_entry = as.IDate(g2_entry))]
panel[, year := year(ym)]

mi <- function(d) year(d) * 12L + (month(d) - 1L)

d <- panel[role_entry %in% c("treated", "ring", "far")]
d <- d[is.na(g2_entry) | ym < g2_entry]
d <- d[!(role_entry == "treated" & year(g_entry) < FOCAL_FROM)]
d <- d[role_entry == "treated" | base_2012 == TRUE]
d[, treated := as.integer(role_entry == "treated")]
d[, et := miym - mi(g_entry)]
d[, rel := fifelse(treated == 1L, pmax(-PRE, pmin(POST, et)), -1L)]
d[, post := fifelse(treated == 1L & !is.na(g_entry) & ym >= g_entry, 1L, 0L)]

# ==============================================================================
# 2. moderator: price rank inside the comuna, frozen on the pre-entry window
# ==============================================================================

rank_local <- function(geo) {
  r <- panel[!is.na(p93), .(station_key, miym, p93, geo = get(geo))]
  r[, n_geo := .N, by = .(geo, miym)]
  r <- r[n_geo >= MIN_N]
  # percentile rank: 0 = cheapest of its market, 1 = most expensive
  r[, rk := (frank(p93, ties.method = "average") - 1) / (.N - 1),
    by = .(geo, miym)]
  r[, .(station_key, miym, rk)]
}

rk_com <- rank_local("comuna")
rk_reg <- rank_local("region")   # kept as a robustness moderator in the csv

pre <- d[post == 0L, .(station_key, miym, treated)]
pre <- merge(pre, rk_com, by = c("station_key", "miym"), all.x = TRUE)
setnames(pre, "rk", "rk_com")
pre <- merge(pre, rk_reg, by = c("station_key", "miym"), all.x = TRUE)
setnames(pre, "rk", "rk_reg")

mod <- pre[, .(rk_com = mean(rk_com, na.rm = TRUE),
               rk_reg = mean(rk_reg, na.rm = TRUE),
               treated = first(treated)), by = station_key]
mod <- mod[!is.na(rk_com)]

# the cutoff comes from the TREATED distribution, because the object of interest
# is heterogeneity among treated stations; controls are placed by the same
# cutoff so "high" means the same thing in both arms
q <- quantile(mod[treated == 1L, rk_com], CORTE, na.rm = TRUE)
mod[, alto := as.integer(rk_com > q)]
mod[, grupo := factor(lbl[alto + 1L], levels = lbl)]

cat("\n=== CORTE ===\n")
cat("percentil", CORTE, "del rank comunal entre tratadas:", round(q, 3), "\n")

cat("\n=== COMPOSICION: estaciones por grupo ===\n")
comp <- mod[, .N, by = .(grupo, treated)]
comp <- dcast(comp, grupo ~ treated, value.var = "N")
setnames(comp, c("grupo", "control", "tratadas"))
print(comp)
fwrite(comp, file.path(TAB, "het_composicion.csv"))

d <- merge(d, mod[, .(station_key, rk_com, rk_reg, alto, grupo)],
           by = "station_key", all.x = TRUE)

# ==============================================================================
# 3. att by group, the difference between them, and the event study
#
# explicit dummies rather than i(): post is zero for every control, so
# post x group gives one clean att per group with the common time path still
# pinned down by the whole control pool
# ==============================================================================

estimar <- function(fv) {
  en_log <- substr(fv, 1, 1) == "p"
  lhs <- if (en_log) sprintf("log(%s)", fv) else fv
  sc  <- if (en_log) 100 else 1

  dd <- d[!is.na(get(fv)) & !is.na(alto)]
  dd[, `:=`(post_bajo = post * (alto == 0L), post_alto = post * alto,
            trat_bajo = treated * (alto == 0L), trat_alto = treated * alto)]

  # one att per group
  m_g <- feols(as.formula(sprintf("%s ~ post_bajo + post_alto | %s", lhs, FE)),
               data = dd, cluster = ~comuna)
  att <- tidy_coefs(m_g)[, .(
    outcome = fv, grupo = lbl[(term == "post_alto") + 1L],
    att = estimate * sc, se = se * sc, p)]

  # same model reparameterised: post is the effect for the low/medium group and
  # post:alto is the DIFFERENCE, so its p-value tests the heterogeneity directly
  m_d <- feols(as.formula(sprintf("%s ~ post + post:alto | %s", lhs, FE)),
               data = dd, cluster = ~comuna)
  ct <- coeftable(m_d)
  dif <- data.table(outcome = fv, base = ct["post", 1] * sc,
                    dif = ct["post:alto", 1] * sc,
                    se = ct["post:alto", 2] * sc, p = ct["post:alto", 4])

  # event study, one path per group
  m_es <- feols(as.formula(sprintf(
    "%s ~ i(rel, trat_bajo, ref = -1) + i(rel, trat_alto, ref = -1) | %s",
    lhs, FE)), data = dd, cluster = ~comuna)
  ce <- tidy_coefs(m_es)[grepl("^rel::", term)]
  ce[, alto_k := as.integer(grepl("trat_alto$", term))]
  ce[, event_time := as.integer(sub("^rel::(-?\\d+).*", "\\1", term))]
  es <- ce[, .(outcome = fv, grupo = lbl[alto_k + 1L], event_time,
               estimate = estimate * sc, se = se * sc)]
  es <- rbind(es, CJ(outcome = fv, grupo = lbl, event_time = -1L,
                     estimate = 0, se = 0, unique = TRUE))
  es[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

  # pre-trend slope per group: mean reversion would drift the high group down
  # before the entry, and that is exactly what must not be there
  b <- coef(m_es); V <- vcov(m_es)
  pend <- rbindlist(lapply(0:1, function(k) {
    sfx <- if (k == 1L) "trat_alto" else "trat_bajo"
    nm <- grep(sprintf("^rel::-.*%s$", sfx), names(b), value = TRUE)
    tt <- as.integer(sub("^rel::(-?\\d+).*", "\\1", nm))
    o <- order(tt); nm <- nm[o]; tt <- tt[o]
    wt <- c(tt, -1L) - mean(c(tt, -1L))
    wt <- (wt / sum(wt^2))[-(length(tt) + 1L)]
    e <- sum(wt * b[nm]) * sc
    s <- sqrt(as.numeric(t(wt) %*% V[nm, nm] %*% wt)) * sc
    data.table(outcome = fv, grupo = lbl[k + 1L], pend = e, pend_se = s,
               pend_p = 2 * pnorm(-abs(e / s)))
  }))

  list(att = att, dif = dif, es = es, pend = pend)
}

res  <- lapply(c(names(precios), names(margenes)), estimar)
att  <- rbindlist(lapply(res, `[[`, "att"))
dif  <- rbindlist(lapply(res, `[[`, "dif"))
es   <- rbindlist(lapply(res, `[[`, "es"))
pend <- rbindlist(lapply(res, `[[`, "pend"))

pend <- merge(pend, att, by = c("outcome", "grupo"), suffixes = c("", "_att"))
pend[, share_att := pend * 12 / att]

fwrite(att,  file.path(TAB, "het_att.csv"))
fwrite(dif,  file.path(TAB, "het_diferencia.csv"))
fwrite(pend, file.path(TAB, "het_pretendencia.csv"))

redondear <- function(x) x[, lapply(.SD, function(z)
  if (is.numeric(z)) round(z, 3) else z)]

cat("\n=== ATT POR GRUPO ===\n")
print(redondear(dcast(att, grupo ~ outcome, value.var = "att")))
cat("\n--- ee ---\n")
print(redondear(dcast(att, grupo ~ outcome, value.var = "se")))

cat("\n=== DIFERENCIA (alto menos bajo/medio) y su p ===\n")
cat("este p es el test de heterogeneidad: si es chico, el efecto SI depende\n")
cat("de cuanto margen tenia la incumbente para ceder\n")
print(redondear(dif[, .(outcome, base, dif, se, p)]))

cat("\n=== PRE-TENDENCIA POR GRUPO (reversion a la media) ===\n")
cat("un p alto en el grupo alto descarta que la caida sea reversion\n")
print(redondear(dcast(pend, grupo ~ outcome, value.var = "pend_p")))
cat("--- deriva 12m / att (no interpretable si el att es ~0) ---\n")
print(redondear(dcast(pend, grupo ~ outcome, value.var = "share_att")))

# ==============================================================================
# 4. robustness: is this really about margin position?
#
#   (a) competitor count as a SECOND moderator. the pre-entry number of rivals
#       is orthogonal to the margin rank in this sample, so the margin gradient
#       should survive letting the effect vary with it too. note it enters
#       interacted with post and measured pre-entry: the contemporaneous count
#       is what entry changes, so as a plain regressor it would be a bad
#       control that absorbs the treatment itself.
#
#   (b) sensitivity to MIN_N. the rank is computed inside comuna-month, so
#       comuna density differences out of the moderator by construction, but
#       "top tercile" of a comuna with four stations is a noisier object than
#       of one with twenty.
# ==============================================================================

sloc <- unique(panel[!is.na(lat), .(station_key, lat, lon)])
a1 <- sloc$lat * pi / 180; o1 <- sloc$lon * pi / 180
hv <- sin(outer(a1, a1, "-") / 2)^2 +
  outer(cos(a1), cos(a1)) * sin(outer(o1, o1, "-") / 2)^2
hv[hv > 1] <- 1
Dm <- 2 * 6371 * asin(sqrt(hv))
diag(Dm) <- Inf
ix <- which(Dm <= RTREAT, arr.ind = TRUE)
edges <- data.table(station_key = sloc$station_key[ix[, "row"]],
                    vecina = sloc$station_key[ix[, "col"]])
presencia <- unique(panel[, .(vecina = station_key, miym)])
ncomp_mes <- merge(edges, presencia, by = "vecina",
                   allow.cartesian = TRUE)[, .N, by = .(station_key, miym)]
setnames(ncomp_mes, "N", "ncomp")

ncomp_pre <- merge(d[post == 0L, .(station_key, miym)], ncomp_mes,
                   by = c("station_key", "miym"))[
  , .(ncomp = mean(ncomp)), by = station_key]

cat("\n=== CORRELACION rank comunal vs n competidores (tratadas, pre-entrada) ===\n")
cc <- merge(mod[treated == 1L, .(station_key, rk_com)], ncomp_pre,
            by = "station_key")
cat(round(cc[, cor(rk_com, ncomp)], 3), "\n")

d <- merge(d, ncomp_pre, by = "station_key", all.x = TRUE)
d[, ncomp_z := (ncomp - mean(ncomp, na.rm = TRUE)) / sd(ncomp, na.rm = TRUE)]

rob_ncomp <- rbindlist(lapply(c(names(precios), names(margenes)), function(fv) {
  en_log <- substr(fv, 1, 1) == "p"
  lhs <- if (en_log) sprintf("log(%s)", fv) else fv
  sc  <- if (en_log) 100 else 1
  dd <- d[!is.na(get(fv)) & !is.na(alto) & !is.na(ncomp_z)]
  m <- feols(as.formula(sprintf("%s ~ post + post:alto + post:ncomp_z | %s",
                                lhs, FE)), data = dd, cluster = ~comuna)
  ct <- coeftable(m)
  data.table(outcome = fv,
             dif_alto = ct["post:alto", 1] * sc, se_alto = ct["post:alto", 2] * sc,
             p_alto = ct["post:alto", 4],
             grad_ncomp = ct["post:ncomp_z", 1] * sc,
             se_ncomp = ct["post:ncomp_z", 2] * sc,
             p_ncomp = ct["post:ncomp_z", 4])
}))
fwrite(rob_ncomp, file.path(TAB, "het_robustez_ncomp.csv"))

cat("\n=== ROBUSTEZ (a): diferencia por margen, controlando el gradiente por competidores ===\n")
print(merge(dif[, .(outcome, dif_base = dif, p_base = p)],
            rob_ncomp[, .(outcome, dif_ctrl = dif_alto, p_ctrl = p_alto,
                          grad_ncomp, p_ncomp)],
            by = "outcome")[, lapply(.SD, function(z)
              if (is.numeric(z)) round(z, 3) else z)])

rob_minn <- rbindlist(lapply(c(4L, 8L, 12L), function(mn) {
  r <- panel[!is.na(p93), .(station_key, miym, p93, geo = comuna)]
  r[, n_geo := .N, by = .(geo, miym)]
  r <- r[n_geo >= mn]
  r[, rk := (frank(p93, ties.method = "average") - 1) / (.N - 1),
    by = .(geo, miym)]
  mm <- merge(d[post == 0L, .(station_key, miym, treated)],
              r[, .(station_key, miym, rk)], by = c("station_key", "miym"))[
    , .(rk = mean(rk), treated = first(treated)), by = station_key]
  qq <- quantile(mm[treated == 1L, rk], CORTE, na.rm = TRUE)
  mm[, alto2 := as.integer(rk > qq)]
  dd0 <- merge(d, mm[, .(station_key, alto2)], by = "station_key")
  rbindlist(lapply(c("p93", "pdi"), function(fv) {
    dd <- dd0[!is.na(get(fv)) & !is.na(alto2)]
    ct <- coeftable(feols(as.formula(sprintf("log(%s) ~ post + post:alto2 | %s",
                                             fv, FE)),
                          data = dd, cluster = ~comuna))
    data.table(min_n = mn, outcome = fv, n_est = uniqueN(dd$station_key),
               dif = ct["post:alto2", 1] * 100, se = ct["post:alto2", 2] * 100,
               p = ct["post:alto2", 4])
  }))
}))
fwrite(rob_minn, file.path(TAB, "het_robustez_minn.csv"))

cat("\n=== ROBUSTEZ (b): sensibilidad al minimo de estaciones por comuna-mes ===\n")
print(rob_minn[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

# ==============================================================================
# 5. latex table and figures
# ==============================================================================

cols <- c(p93 = "Precio 93 (\\%)", p95 = "Precio 95 (\\%)",
          pdi = "Precio diésel (\\%)", m93 = "Margen 93 (\\$/L)",
          mdi = "Margen diésel (\\$/L)")

w <- dcast(att, grupo ~ outcome, value.var = c("att", "se", "p"))
setorder(w, grupo)
out <- data.table(Grupo = as.character(w$grupo))
for (fv in names(cols)) {
  out[[cols[[fv]]]] <- celda_tex(w[[paste0("att_", fv)]],
                                 w[[paste0("se_", fv)]],
                                 w[[paste0("p_", fv)]],
                                 if (substr(fv, 1, 1) == "p") 3 else 2)
}
wd <- dcast(dif, . ~ outcome, value.var = c("dif", "se", "p"))
out <- rbind(out, c(list(Grupo = "Diferencia"),
                    setNames(lapply(names(cols), function(fv)
                      celda_tex(wd[[paste0("dif_", fv)]],
                                wd[[paste0("se_", fv)]],
                                wd[[paste0("p_", fv)]],
                                if (substr(fv, 1, 1) == "p") 3 else 2)),
                      unname(cols))))
guardar_tabla_tex(out, "tab_het_rank.tex", regla_antes_de = "Diferencia")

pal <- c(AZUL_CL, NARANJO_OSC)

for (tipo in c("precio", "margen")) {
  vars <- if (tipo == "precio") precios else margenes
  dd <- es[outcome %in% names(vars)]
  dd[, combustible := factor(outcome, levels = names(vars), labels = vars)]
  dd[, grupo := factor(grupo, levels = lbl)]
  p <- ggplot(dd, aes(event_time, estimate, color = grupo, fill = grupo)) +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.55) + geom_point(size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey55") +
    facet_wrap(~combustible, scales = "free_y") +
    scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
    labs(x = "Meses desde la entrada (extremos agrupados en ±12)",
         y = if (tipo == "precio") "Efecto sobre el precio (%)"
             else "Efecto sobre el margen ($/L)",
         color = NULL, fill = NULL,
         title = "El efecto de la entrada cae sobre quien tenía margen que ceder",
         subtitle = paste("Grupos según la posición de la estación en la",
                          "distribución de precios de su comuna antes de la entrada")) +
    tema()
  guardar(p, file.path(FIG, sprintf("het_rank_%s.pdf", tipo)),
          9.5, if (tipo == "precio") 5.8 else 4.2)
}

# the competitor-count moderator was dropped: it showed no gradient in prices
# and the wrong sign in margins, and its outputs would be stale here
unlink(file.path(c(FIG, FIG, TAB), c("het_ncomp_precio.pdf",
                                     "het_ncomp_margen.pdf",
                                     "tab_het_ncomp.tex")))
unlink(file.path(TAB, "het_gradiente.csv"))

message("14_heterogeneidad.R: figuras en ", FIG, " y tablas en ", TAB, ".")
