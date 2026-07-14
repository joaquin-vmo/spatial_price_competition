
# objective: "rockets and feathers" -- is cost pass-through faster/more
# complete when the mepco wholesale price rises than when it falls? the
# central question in this literature (bacon 1991; borenstein-cameron-
# gilbert 1997; bachmeier-griffin 2003; peltzman 2000). same local
# projections machinery as 06_local_projections.R (pooled over every mepco
# update 2014-2026, no date fe -- see that script's header note on why),
# just split the sample by the sign of dlog_wholesale instead of by brand
# or competition.

library(data.table)
library(fixest)
library(ggplot2)

db <- fread(file.path("data/processed", "database.csv"),
            encoding = "UTF-8", colClasses = list(character = "id"))
mepco <- fread(file.path("data/processed", "whosale_prices_mepco.csv"))

H <- 12
fuels <- c("93", "97", "di")

dbf <- db[fuel %in% fuels & !is.na(wholesale_w)]

#####

# resample onto the mepco weekly calendar (locf), same construction as 06
mepco_dates <- sort(unique(as.Date(mepco$date)))
id_fuel_bounds <- dbf[, .(min_date = min(date), max_date = max(date)), by = .(id, fuel)]
grid <- id_fuel_bounds[, .(date = mepco_dates[mepco_dates >= min_date & mepco_dates <= max_date]),
                        by = .(id, fuel)]

price_obs <- dbf[, .(id, fuel, date, price)]
setkey(price_obs, id, fuel, date)
setkey(grid, id, fuel, date)
grid <- price_obs[grid, roll = TRUE]

mepco_long <- rbindlist(lapply(fuels, function(f) {
  data.table(fuel = f, date = mepco$date, wholesale_w = mepco[[paste0(f, "_w/")]])
}))
setorder(mepco_long, fuel, date)
mepco_long[, log_w := log(wholesale_w)]
mepco_long[, dlog_wholesale := log_w - shift(log_w), by = fuel]

# cumulative wholesale change over [t-1, t+h], same window as y_h below.
# using only the origin week's dlog_wholesale as the regressor at every
# horizon is the wrong specification here: it's correlated ~0.40-0.44 with
# how much wholesale keeps moving over the following weeks (mepco doesn't
# mean-revert), so that single-week shock acts as a proxy for a larger
# cumulative move and inflates beta_h at longer horizons (checked: gave
# beta_h > 2 for "baja" at h=8, not economically sensible). regressing on
# the cumulative change over the same window as the price response is the
# standard fix (as in exchange-rate pass-through regressions).
for (h in 0:H) {
  mepco_long[, paste0("cumw_h", h) := shift(log_w, -h) - shift(log_w, 1), by = fuel]
}

# decomposition for the robustness section at the bottom: split each weekly
# wholesale step into its positive and negative part, then accumulate each
# separately over the window [t, t+h]. cumw_up_h = sum of the weeks that rose,
# cumw_dn_h = sum of the weeks that fell (<= 0). by construction
# cumw_up_h + cumw_dn_h == cumw_h. this is the borenstein-cameron-gilbert
# (1997) way to test asymmetry in one regression, without splitting the sample
# on the origin-week sign (which conditions on a systematically different
# cost path and can manufacture a spurious asymmetry)
mepco_long[, dpos := pmax(dlog_wholesale, 0)]
mepco_long[, dneg := pmin(dlog_wholesale, 0)]
mepco_long[, cumw_up_h0 := dpos, by = fuel]
mepco_long[, cumw_dn_h0 := dneg, by = fuel]
for (h in 1:H) {
  mepco_long[, paste0("cumw_up_h", h) := get(paste0("cumw_up_h", h - 1)) + shift(dpos, -h), by = fuel]
  mepco_long[, paste0("cumw_dn_h", h) := get(paste0("cumw_dn_h", h - 1)) + shift(dneg, -h), by = fuel]
}

grid <- merge(grid, mepco_long, by = c("fuel", "date"), all.x = TRUE)

setorder(grid, id, fuel, date)
grid[, log_price := log(price)]
for (h in 0:H) {
  grid[, paste0("y_h", h) := shift(log_price, -h) - shift(log_price, 1), by = .(id, fuel)]
}

# sign of that week's shock; drop exact-zero weeks (uninformative either way)
grid[, sentido := fifelse(dlog_wholesale > 0, "sube", fifelse(dlog_wholesale < 0, "baja", NA_character_))]

# heterogeneity dimensions to interact with the sign of the shock, same
# construction as 06_distributed_lags.R / 05_local_projection.R. the point is
# to ask whether the rockets-and-feathers asymmetry (faster up than down) is
# disciplined by local competition or differs across brands.

# brand: modal distributor over the station's history, big-3-vs-rest
brand <- dbf[, .N, by = .(id, distributor)]
setorder(brand, id, -N)
brand <- brand[, .SD[1], by = id]
brand[, brand_group := fifelse(distributor %in% c("copec", "shell", "aramco_petrobras"),
                                distributor, "independiente")]
grid <- merge(grid, brand[, .(id, brand_group)], by = "id", all.x = TRUE)

# n° competitors within 2km, resampled onto the grid (locf), then terciles.
# labelled poca/media/mucha (not baja/media/alta) so "baja" doesn't collide
# with sentido == "baja" (the cost falling) once both appear in this script
comp_hist <- unique(dbf[, .(id, date, competitors_2km)])
setkey(comp_hist, id, date)
setkey(grid, id, date)
grid <- comp_hist[grid, roll = TRUE]
comp_breaks <- quantile(grid$competitors_2km, probs = 0:3 / 3, na.rm = TRUE)
grid[, comp_nivel := cut(competitors_2km, breaks = unique(comp_breaks),
                          include.lowest = TRUE, labels = c("poca", "media", "mucha"))]
setorder(grid, id, fuel, date)

#####

run_lp <- function(data, group_col) {
  data <- data[!is.na(get(group_col))]
  groups <- split(data, data[[group_col]])
  rbindlist(lapply(names(groups), function(g) {
    d <- groups[[g]]
    rbindlist(lapply(0:H, function(h) {
      x_h <- paste0("cumw_h", h)
      m <- feols(as.formula(paste0("y_h", h, " ~ ", x_h, " | id^fuel")), data = d, cluster = ~id)
      ct <- coeftable(m)
      data.table(horizon = h, group = g,
                 estimate = ct[x_h, "Estimate"],
                 se = ct[x_h, "Std. Error"],
                 n = nrow(d))
    }))
  }))
}

irf <- run_lp(grid, "sentido")
irf[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

fwrite(irf, file.path("results/tables", "rockets_feathers_irf.csv"))
print(irf[horizon %in% c(0, 1, 4, 8, 12)])

# same asymmetry, but now within each level of a heterogeneity dimension
# (brand or competition): one lp per sentido x het_col cell. sentido stays as
# its own column (mapped to colour) and het_col as another (mapped to facet),
# so the plot can show up-vs-down asymmetry within each subgroup
run_lp_het <- function(data, het_col) {
  d0 <- data[!is.na(sentido) & !is.na(get(het_col))]
  combos <- unique(d0[, .(sentido, het = get(het_col))])
  rbindlist(lapply(seq_len(nrow(combos)), function(i) {
    s <- combos$sentido[i]
    hv <- combos$het[i]
    d <- d0[sentido == s & get(het_col) == hv]
    rbindlist(lapply(0:H, function(h) {
      x_h <- paste0("cumw_h", h)
      m <- feols(as.formula(paste0("y_h", h, " ~ ", x_h, " | id^fuel")), data = d, cluster = ~id)
      ct <- coeftable(m)
      data.table(horizon = h, sentido = s, het = hv,
                 estimate = ct[x_h, "Estimate"],
                 se = ct[x_h, "Std. Error"],
                 n = nrow(d))
    }))
  }))
}

#####

p <- ggplot(irf, aes(horizon, estimate, color = group, fill = group)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c(sube = "firebrick", baja = "steelblue")) +
  scale_fill_manual(values = c(sube = "firebrick", baja = "steelblue")) +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "sentido del\ncambio mayorista", fill = "sentido del\ncambio mayorista",
    title = "rockets and feathers: pass-through segun el signo del shock de costo",
    subtitle = "linea punteada en 1 = pass-through completo; usa todas las actualizaciones mepco 2014-2026"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "rockets_feathers_irf.pdf"), p, width = 7.5, height = 4.5)

#####

# asymmetry by competition and by brand

irf_comp <- run_lp_het(grid, "comp_nivel")[, dimension := "competencia"]
irf_brand <- run_lp_het(grid, "brand_group")[, dimension := "marca"]

irf_het <- rbind(irf_comp, irf_brand)
irf_het[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)]

fwrite(irf_het, file.path("results/tables", "rockets_feathers_het_irf.csv"))

sentido_scale <- list(
  scale_color_manual(values = c(sube = "firebrick", baja = "steelblue")),
  scale_fill_manual(values = c(sube = "firebrick", baja = "steelblue"))
)

# competition: the rockets-and-feathers story predicts the up-vs-down gap
# (sube above baja) should shrink where competition is stronger ("mucha")
p_comp <- ggplot(irf_het[dimension == "competencia"],
                 aes(horizon, estimate, color = sentido, fill = sentido)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~factor(het, levels = c("poca", "media", "mucha"))) +
  sentido_scale +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "sentido del\ncambio mayorista", fill = "sentido del\ncambio mayorista",
    title = "rockets and feathers por nivel de competencia (competidores a 2km)",
    subtitle = "asimetria = brecha entre sube (rojo) y baja (azul); linea punteada en 1 = pass-through completo"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "rockets_feathers_competencia.pdf"), p_comp, width = 9.5, height = 4.5)

# brand
p_brand <- ggplot(irf_het[dimension == "marca"],
                  aes(horizon, estimate, color = sentido, fill = sentido)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~het) +
  sentido_scale +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": " * Delta * log(precio)[t+h] - Delta * log(precio)[t-1] * " por " * Delta * log(mayorista)),
    color = "sentido del\ncambio mayorista", fill = "sentido del\ncambio mayorista",
    title = "rockets and feathers por marca",
    subtitle = "asimetria = brecha entre sube (rojo) y baja (azul); linea punteada en 1 = pass-through completo"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "rockets_feathers_marca.pdf"), p_brand, width = 9.5, height = 4.5)

#####

# ROBUSTNESS: is the reverse-asymmetry sign above (pass-through larger when
# cost FALLS than when it rises) real, or an artifact of splitting the sample
# on the origin-week sign? re-estimate without any sample split, regressing
# y_h on the cumulative up-moves and cumulative down-moves of the cost jointly
# (see the cumw_up_h / cumw_dn_h construction above):
#
#   y_h ~ beta_up * cumw_up_h + beta_down * cumw_dn_h | id^fuel
#
# beta_up is the pass-through of cost *increases*, beta_down of *decreases*.
# textbook rockets-and-feathers => beta_up > beta_down (asym > 0). the reverse
# sign found above => asym < 0. symmetric (artifact) => asym ~ 0. the asymmetry
# and its se come from the joint clustered vcov, so it's a real test on the
# difference, not an eyeballed gap between two separate subsample fits

run_asym_decomp <- function(data, group_col = NULL) {
  if (is.null(group_col)) {
    groups <- list(pooled = data)
  } else {
    data <- data[!is.na(get(group_col))]
    groups <- split(data, data[[group_col]])
  }
  rbindlist(lapply(names(groups), function(g) {
    d <- groups[[g]]
    rbindlist(lapply(0:H, function(h) {
      xu <- paste0("cumw_up_h", h)
      xd <- paste0("cumw_dn_h", h)
      m <- feols(as.formula(paste0("y_h", h, " ~ ", xu, " + ", xd, " | id^fuel")),
                 data = d, cluster = ~id)
      b <- coef(m)
      V <- vcov(m)
      asym <- b[[xu]] - b[[xd]]
      se_asym <- sqrt(V[xu, xu] + V[xd, xd] - 2 * V[xu, xd])
      data.table(horizon = h, group = g,
                 beta_up = b[[xu]], se_up = sqrt(V[xu, xu]),
                 beta_down = b[[xd]], se_down = sqrt(V[xd, xd]),
                 asym = asym, se_asym = se_asym,
                 t_asym = asym / se_asym, n = nrow(d))
    }))
  }))
}

asym_pooled <- run_asym_decomp(grid)[, dimension := "pooled"]
asym_comp <- run_asym_decomp(grid, "comp_nivel")[, dimension := "competencia"]
asym_brand <- run_asym_decomp(grid, "brand_group")[, dimension := "marca"]

asym_all <- rbind(asym_pooled, asym_comp, asym_brand)
asym_all[, `:=`(asym_low = asym - 1.96 * se_asym, asym_high = asym + 1.96 * se_asym)]

fwrite(asym_all, file.path("results/tables", "rockets_feathers_robustness_irf.csv"))
cat("\n=== decomposition asymmetry (asym = beta_up - beta_down; >0 = textbook rocket) ===\n")
print(asym_all[group == "pooled" & horizon %in% c(0, 1, 4, 8, 12),
               .(horizon, beta_up, beta_down, asym, se_asym, t_asym)])

#####

# graphs

# pooled: beta_up (response to cost increases) vs beta_down (to decreases),
# on the same axis. textbook rockets => red (up) above blue (down)
p_rob_pooled <- ggplot(rbind(
    asym_pooled[, .(horizon, sentido = "sube", estimate = beta_up, se = se_up)],
    asym_pooled[, .(horizon, sentido = "baja", estimate = beta_down, se = se_down)]
  )[, `:=`(ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se)],
  aes(horizon, estimate, color = sentido, fill = sentido)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  sentido_scale +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression(beta[h] * ": respuesta del precio al costo acumulado (subas vs bajas)"),
    color = "sentido del\ncambio mayorista", fill = "sentido del\ncambio mayorista",
    title = "robustez rockets and feathers: descomposicion subas/bajas (sin partir la muestra)",
    subtitle = "beta_up (rojo) vs beta_down (azul) en una sola regresion; textbook => rojo por encima de azul"
  ) +
  theme_minimal()

ggsave(file.path("results/figures", "rockets_feathers_robustness.pdf"), p_rob_pooled, width = 7.5, height = 4.5)

# the asymmetry itself (beta_up - beta_down) with a proper 95% ci, pooled and
# by subgroup. ci crossing 0 => no significant asymmetry at that horizon
asym_all[, group_lab := factor(group, levels = c("pooled", "poca", "media", "mucha",
                                                 "aramco_petrobras", "copec", "independiente", "shell"))]

p_rob_asym <- ggplot(asym_all, aes(horizon, asym, color = dimension, fill = dimension)) +
  geom_ribbon(aes(ymin = asym_low, ymax = asym_high), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~group_lab, nrow = 2) +
  labs(
    x = "semanas desde el cambio de costo mayorista",
    y = expression("asimetria " * beta[up] - beta[down] * " (>0 = textbook rocket)"),
    color = "dimension", fill = "dimension",
    title = "robustez: asimetria estimada con su ci 95%, pooled y por subgrupo",
    subtitle = "banda que cruza 0 = asimetria no significativa a ese horizonte"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path("results/figures", "rockets_feathers_robustness_asym.pdf"), p_rob_asym, width = 10, height = 5.5)

#####

# SUMMARY OF FINDINGS (robustness run 2026-07-14)
#
# 1. the reverse-asymmetry sign is REAL, not an artifact of splitting the
#    sample on the origin-week shock sign. the borenstein-cameron-gilbert
#    one-regression decomposition (y_h on cumw_up_h + cumw_dn_h, no split)
#    gives asym = beta_up - beta_down < 0 at every horizon, pooled and in
#    every subgroup, with cis excluding 0 by a wide margin (pooled t_asym
#    ranges from -30 at h=0 to -7 at h=12). so pass-through responds MORE to
#    cost DECREASES than to increases -- the opposite of textbook rockets-and-
#    feathers (which would give asym > 0).
#
# 2. it is mostly an impact / short-run phenomenon: |asym| is largest at h=0-1
#    and fades toward h=8-12, where pass-through converges to roughly the same
#    level whether cost rose or fell. reportable as a clean result.
#
# 3. competition AMPLIFIES the (reverse) asymmetry rather than dampening it:
#    at h=0 asym is -0.15 (poca) / -0.25 (media) / -0.31 (mucha), and by h=8
#    the low-competition tercile has converged to symmetry (asym ~ 0, t ~ 0)
#    while high-competition stays asymmetric (-0.08, t = -7). note this is the
#    OPPOSITE of the search/tacit-collusion prediction (competition -> less
#    asymmetry); mechanism is not obvious and shouldn't be asserted without
#    more structure -- flagged, not resolved.
#
# 4. by brand: shell most asymmetric (h=0 asym = -0.35, t = -47), aramco/copec
#    intermediate, and independientes the only group not distinguishable from
#    symmetry (h=0 asym = -0.03, ci crosses 0) -- consistent with the sample-
#    split result where independientes were the only textbook-signed group.
