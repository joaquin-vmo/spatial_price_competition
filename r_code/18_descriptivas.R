
# objective: the descriptive tables and figures of the thesis -- how many
# stations there are, how the treated and control groups break down, entry and
# exit flows, the price series, and what the local market looks like at the
# moment of entry.
#
# a caveat that governs everything about exits below. the reporting regime
# changes between the legacy files (2012-2022) and the current ones (2023+),
# and the change leaves a scar: 120 stations stop being reported in december
# 2022 alone, against one or two a year afterwards. part of that is a data
# artifact already fixed upstream (the duplicate "a"-suffixed ids folded into
# their base in 10_build_event_panel.R), but the rest are stations that had
# gone dormant at some unobserved earlier date and were purged at the
# migration. exit DATES are therefore not trustworthy around the boundary, and
# every exit table here reports december 2022 separately instead of hiding it
# inside an annual total. entry dates do not have this problem: a station that
# appears for the first time is genuinely new.
#
# TAKES:  data/processed/panel_mensual.csv, entradas.csv
# PRODUCES:
#   results/tables/desc_*.csv and tab_desc_*.tex
#   results/figures/desc_*.pdf

library(data.table)
library(ggplot2)

FIG <- "results/figures"
TAB <- "results/tables"

RTREAT <- 2
RCTRL  <- 5
# a replacement AT THE SAME SITE cannot show up as an exit here:
# 10_build_event_panel.R already folds a station that succeeds another within
# 50 m into the same station_key, precisely so that a change of operator is not
# counted as an exit plus an entry. what is left to ask is whether the local
# market was replenished, so replacement is measured on a wider radius and
# within a bounded window
REEMP_M  <- 500    # metres
REEMP_MESES <- 24L # months after the exit
CENSURA <- 3L      # months before the panel end within which an ending spell
                   # is right-censored rather than an exit
FOCAL_FROM <- 2014L
R_EARTH <- 6371

fuente <- "Fuente: elaboración propia a partir de la base de precios de la CNE."

precios <- c(p93 = "Gasolina 93", p95 = "Gasolina 95",
             p97 = "Gasolina 97", pdi = "Diésel")

guardar_tabla_tex <- function(d, archivo, align = NULL) {
  if (is.null(align)) align <- paste0("l", strrep("r", ncol(d) - 1L))
  cuerpo <- apply(d, 1, function(r) paste(paste(r, collapse = " & "), "\\\\"))
  writeLines(
    c(sprintf("\\begin{tabular}{%s}", align), "\\toprule",
      paste(paste(names(d), collapse = " & "), "\\\\"), "\\midrule",
      cuerpo, "\\bottomrule", "\\end{tabular}"),
    file.path(TAB, archivo)
  )
}

# palatino es la tipografia del documento. el dispositivo pdf() por defecto no
# la embebe y descarta los simbolos unicode de las etiquetas; quartz, el motor
# nativo de macos, si lo hace. fuera de macos se cae de vuelta a ggsave
TIPO <- "Palatino"

tema <- function(base = 14) {
  theme_minimal(base_size = base, base_family = TIPO) +
    theme(text = element_text(family = TIPO),
          panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.caption = element_text(colour = "grey40", size = rel(0.75)))
}

guardar <- function(fig, archivo, width, height) {
  if (capabilities("aqua")) {
    grDevices::quartz(file = file.path(FIG, archivo), type = "pdf",
                      width = width, height = height)
    print(fig)
    grDevices::dev.off()
  } else {
    ggsave(file.path(FIG, archivo), fig, width = width, height = height)
  }
}

# ==============================================================================
# 1. inputs
# ==============================================================================

panel <- fread("data/processed/panel_mensual.csv")
panel[, `:=`(ym = as.IDate(ym), g_entry = as.IDate(g_entry),
             g2_entry = as.IDate(g2_entry))]
entradas <- fread("data/processed/entradas.csv")[, g := as.IDate(g)]
entradas[, grupo := fcase(
  edist == "copec", "Copec",
  edist == "shell", "Shell",
  edist == "aramco_petrobras", "Aramco/Petrobras",
  default = "Pequeña/sin bandera")]
grupos_orden <- c("Copec", "Shell", "Aramco/Petrobras", "Pequeña/sin bandera")

fin_panel <- max(panel$miym)

st <- unique(panel[, .(station_key, region, comuna, distribuidor, is_franchise,
                       base_2012, role_entry, g_entry, lat, lon)])
st <- st[, .SD[1], by = station_key]

# ==============================================================================
# 2. spells, entries and exits
#
# the monthly panel already drops months a station did not report within the
# carry-forward window, so a gap in it is a gap in operation (or in reporting).
# a spell is a maximal run of consecutive months
# ==============================================================================

setorder(panel, station_key, miym)
panel[, gap := miym - shift(miym) - 1L, by = station_key]
panel[, spell := 1L + cumsum(!is.na(gap) & gap > 0L), by = station_key]

sp <- panel[, .(ini = min(miym), fin = max(miym)), by = .(station_key, spell)]
sp[, n_sp := max(spell), by = station_key]
sp[, ultimo := spell == n_sp]

# a spell that ends within CENSURA months of the panel end is censored, not an
# exit: the station may well still be operating
salidas <- sp[!(ultimo & fin >= fin_panel - CENSURA)]
# a spell that ends and is followed by another one is a reporting interruption,
# not a closure: the station comes back. the two are kept apart throughout
salidas[, tipo := fifelse(ultimo, "Salida definitiva",
                          "Interrupción de reporte")]
salidas[, `:=`(anio = fin %/% 12L, mes = fin %% 12L + 1L)]
salidas <- merge(salidas, st[, .(station_key, region, distribuidor, lat, lon)],
                 by = "station_key")

# replacement: the station itself comes back (temporary), or another station
# starts up within REEMP_M metres afterwards
ini_st <- panel[, .(ini = min(miym)), by = station_key]
loc <- merge(unique(panel[!is.na(lat), .(station_key, lat, lon)])[
  , .SD[1], by = station_key], ini_st, by = "station_key")

a1 <- loc$lat * pi / 180; o1 <- loc$lon * pi / 180
h <- sin(outer(a1, a1, "-") / 2)^2 +
  outer(cos(a1), cos(a1)) * sin(outer(o1, o1, "-") / 2)^2
h[h > 1] <- 1
D <- 2 * R_EARTH * asin(sqrt(h)) * 1000
diag(D) <- Inf

idx <- match(salidas$station_key, loc$station_key)
salidas[, geo_reemp := vapply(seq_len(.N), function(k) {
  i <- idx[k]
  if (is.na(i)) return(FALSE)
  any(D[i, ] < REEMP_M & loc$ini > salidas$fin[k] &
        loc$ini <= salidas$fin[k] + REEMP_MESES)
}, logical(1))]
salidas[, reemplazo := fifelse(geo_reemp, "Con reemplazo", "Sin reemplazo")]
# the december 2022 wave is the migration scar, not a market event
salidas[, artefacto := anio == 2022L & mes == 12L]

# ==============================================================================
# 3. summary table: what the sample looks like
# ==============================================================================

resumen <- data.table(
  Concepto = c(
    "Estaciones observadas (2012–2026)",
    "\\quad presentes desde 2012 (censura izquierda)",
    "\\quad entrantes posteriores a 2012",
    "Entradas registradas",
    "\\quad focales (2014 en adelante)",
    "Incumbentes tratadas (entrada a $<$2 km)",
    "Controles anillo (2--5 km)",
    "Controles lejanos ($>$5 km)",
    "\\quad total nunca tratadas",
    "Salidas definitivas (excl. dic-2022)",
    "Interrupciones de reporte ($>$3 meses)",
    "\\quad salidas en dic-2022 (artefacto de fuente)"),
  N = c(
    uniqueN(panel$station_key),
    st[base_2012 == TRUE, .N],
    st[base_2012 == FALSE, .N],
    nrow(entradas),
    entradas[year(g) >= FOCAL_FROM, .N],
    st[role_entry == "treated", .N],
    st[role_entry == "ring", .N],
    st[role_entry == "far", .N],
    st[role_entry %in% c("ring", "far"), .N],
    salidas[tipo == "Salida definitiva" & artefacto == FALSE, .N],
    salidas[tipo == "Interrupción de reporte", .N],
    salidas[artefacto == TRUE, .N]))

fwrite(resumen, file.path(TAB, "desc_resumen.csv"))
guardar_tabla_tex(resumen, "tab_desc_resumen.tex")
cat("\n=== RESUMEN DE LA MUESTRA ===\n"); print(resumen)

# ==============================================================================
# 4. entries and exits by year
# ==============================================================================

flujo_e <- entradas[, .(N = .N), by = .(anio = year(g))][, tipo := "Entradas"]
flujo_s <- salidas[artefacto == FALSE, .N, by = .(anio, tipo)]
flujo_a <- salidas[artefacto == TRUE, .(N = .N), by = anio][
  , tipo := "Salidas dic-2022 (artefacto)"]

flujo <- rbind(flujo_e, flujo_s, flujo_a)
flujo_w <- dcast(flujo, anio ~ tipo, value.var = "N", fill = 0L)
setorder(flujo_w, anio)

# active stations at the end of each year, which makes the source break visible
parque <- panel[, .(activas = uniqueN(station_key)), by = .(anio = year(ym), ym)]
parque <- parque[parque[, .I[which.max(ym)], by = anio]$V1, .(anio, activas)]
flujo_w <- merge(flujo_w, parque, by = "anio", all.x = TRUE)

fwrite(flujo_w, file.path(TAB, "desc_flujos_por_anio.csv"))
flujo_tex <- copy(flujo_w)
setnames(flujo_tex, c("anio", "activas"), c("Año", "Activas"))
guardar_tabla_tex(flujo_tex[, lapply(.SD, as.character)],
                  "tab_desc_flujos_por_anio.tex")
cat("\n=== ENTRADAS Y SALIDAS POR ANIO ===\n"); print(flujo_w)

pal_flujo <- c("Entradas" = "#1a9850",
               "Salida definitiva" = "#d73027",
               "Interrupción de reporte" = "#4575b4",
               "Salidas dic-2022 (artefacto)" = "grey60")
fl <- copy(flujo)
fl[tipo != "Entradas", N := -N]
fl[, tipo := factor(tipo, levels = names(pal_flujo))]

p_flujo <- ggplot(fl, aes(factor(anio), N, fill = tipo)) +
  geom_col() +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = pal_flujo) +
  scale_y_continuous(labels = abs) +
  labs(x = "Año", y = "Número de estaciones", fill = NULL,
       title = "Entradas y salidas de estaciones de servicio, 2012–2026",
       subtitle = "Entradas hacia arriba, salidas hacia abajo",
       caption = paste(fuente,
                       "La barra gris de 2022 corresponde al quiebre de fuente entre los archivos legacy y actuales.")) +
  tema() + theme(panel.grid.major.x = element_blank())
guardar(p_flujo, "desc_flujos_por_anio.pdf", 9.5, 5)

# ==============================================================================
# 4b. entries by year and by region, broken down by the entrant's brand
#
# the 2013 cohort is reported but flagged: 79 entries against 20-50 in a normal
# year is a coverage artifact at the start of the sample, so those events
# contaminate and cut in the design but are never used as treatment
# ==============================================================================

ent_anio <- dcast(entradas[, .N, by = .(anio = year(g), grupo)],
                  anio ~ grupo, value.var = "N", fill = 0L)
setcolorder(ent_anio, c("anio", intersect(grupos_orden, names(ent_anio))))
ent_anio[, Total := rowSums(.SD), .SDcols = -"anio"]
ent_anio[, Focal := fifelse(anio >= FOCAL_FROM, "Sí", "No")]
setorder(ent_anio, anio)
fwrite(ent_anio, file.path(TAB, "desc_entradas_por_anio.csv"))

setnames(ent_anio, "anio", "Año")
tot <- as.list(c(`Año` = "Total",
                 lapply(ent_anio[, -c("Año", "Focal")], sum), Focal = ""))
guardar_tabla_tex(rbind(ent_anio[, lapply(.SD, as.character)],
                        lapply(tot, as.character)),
                  "tab_desc_entradas_por_anio.tex")
setnames(ent_anio, "Año", "anio")
cat("\n=== ENTRADAS POR ANIO Y MARCA DEL ENTRANTE ===\n"); print(ent_anio)

ent_reg <- entradas[year(g) >= FOCAL_FROM & !is.na(eregion), .N,
                    by = .(eregion, grupo)]
ent_reg[, region_lbl := tools::toTitleCase(gsub("_", " ", eregion))]
ent_reg[, total := sum(N), by = region_lbl]
fwrite(dcast(ent_reg, region_lbl + total ~ grupo, value.var = "N", fill = 0L),
       file.path(TAB, "desc_entradas_por_region.csv"))

# colores corporativos de cada cadena, para que la figura se lea sin leyenda;
# el resto va en gris neutro, que no compite con los tres saturados
pal_marca <- c("Copec" = "#156dfd", "Shell" = "#fbce06",
               "Aramco/Petrobras" = "#01af4a", "Pequeña/sin bandera" = "#7f7f7f")

# la figura arranca en 2014: la cohorte 2013 es un artefacto de cobertura y
# ademas el mepco solo rige desde agosto de ese ano, de modo que el periodo
# graficado coincide con el que el diseno efectivamente usa
p_ent_anio <- ggplot(entradas[year(g) >= FOCAL_FROM, .N, by = .(anio = year(g), grupo)],
                     aes(factor(anio), N,
                         fill = factor(grupo, levels = grupos_orden))) +
  geom_col() +
  scale_fill_manual(values = pal_marca) +
  labs(x = "Año", y = "Número de entradas", fill = NULL,
       title = "Entradas de estaciones de servicio por año y marca, 2014–2026",
       caption = fuente) +
  tema() + theme(panel.grid.major.x = element_blank())
guardar(p_ent_anio, "desc_entradas_por_anio.pdf", 9.5, 5)

p_ent_reg <- ggplot(ent_reg, aes(reorder(region_lbl, total), N,
                                 fill = factor(grupo, levels = grupos_orden))) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = pal_marca) +
  labs(x = NULL, y = "Número de entradas", fill = NULL,
       title = "Entradas de estaciones de servicio por región, 2014–2026",
       subtitle = "Según la marca de la estación entrante", caption = fuente) +
  tema() + theme(panel.grid.major.y = element_blank())
guardar(p_ent_reg, "desc_entradas_por_region.pdf", 9, 6.5)

# ==============================================================================
# 5. average price and margin over time
# ==============================================================================

# a month is kept only if the panel spans it end to end. the panel stops on
# 2026-07-12, so july 2026 covers barely two weeks and lands in the middle of an
# unfinished pass-through: mepco fell 95 $/L on 18 june and retail was still
# adjusting, which drags the average margin down to a quarter of its usual level.
# that is a truncation artifact, not a market outcome
ult_dia <- as.IDate("2026-07-12")   # ultima fecha observada en la base de precios
meses_ok <- panel[, unique(ym)]
meses_ok <- meses_ok[meses_ok + 31L <= ult_dia | month(meses_ok) != month(ult_dia) |
                       year(meses_ok) != year(ult_dia)]
completo <- function(d) d[ym %in% meses_ok & ym > min(panel$ym)]

serie <- melt(panel, id.vars = "ym", measure.vars = names(precios),
              variable.name = "comb", value.name = "precio",
              variable.factor = FALSE)
serie <- serie[!is.na(precio), .(precio = mean(precio)), by = .(ym, comb)]
serie <- completo(serie)
serie[, combustible := factor(comb, levels = names(precios), labels = precios)]
fwrite(dcast(serie, ym ~ combustible, value.var = "precio"),
       file.path(TAB, "desc_precio_promedio.csv"))

p_precio <- ggplot(serie, aes(ym, precio, colour = combustible)) +
  geom_line(linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = NULL, y = "Precio promedio ($/litro)", colour = NULL,
       title = "Precio promedio nacional por combustible, 2012–2026",
       subtitle = "Promedio simple entre estaciones activas cada mes",
       caption = fuente) +
  tema()
guardar(p_precio, "desc_precio_promedio.pdf", 9, 5)

# margin against the national mepco quote, which only exists from 2014-08
marg <- melt(panel, id.vars = "ym", measure.vars = c("m93", "m97", "mdi"),
             variable.name = "comb", value.name = "margen",
             variable.factor = FALSE)
marg <- marg[!is.na(margen), .(margen = mean(margen)), by = .(ym, comb)]
marg <- completo(marg)
marg[, combustible := factor(comb, levels = c("m93", "m97", "mdi"),
                             labels = c("Gasolina 93", "Gasolina 97", "Diésel"))]
fwrite(dcast(marg, ym ~ combustible, value.var = "margen"),
       file.path(TAB, "desc_margen_promedio.csv"))

p_margen <- ggplot(marg, aes(ym, margen, colour = combustible)) +
  geom_line(linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = NULL, y = "Margen promedio ($/litro)", colour = NULL,
       title = "Margen bruto promedio sobre el precio mayorista MEPCO, 2014–2026",
       caption = fuente) +
  tema()
guardar(p_margen, "desc_margen_promedio.pdf", 9, 5)

# ==============================================================================
# 6. exits with and without replacement
# ==============================================================================

reemp_tab <- dcast(salidas[artefacto == FALSE], tipo ~ reemplazo, value.var = "station_key",
                   fun.aggregate = length)
reemp_tab[, Total := rowSums(.SD), .SDcols = -"tipo"]
fwrite(reemp_tab, file.path(TAB, "desc_salidas_reemplazo.csv"))
cat("\n=== SALIDAS SEGUN REEMPLAZO (excl. dic-2022) ===\n"); print(reemp_tab)

reemp_reg <- salidas[artefacto == FALSE & !is.na(region), .N, by = .(region, reemplazo)]
reemp_reg[, region_lbl := tools::toTitleCase(gsub("_", " ", region))]
reemp_reg[, total := sum(N), by = region_lbl]
fwrite(dcast(reemp_reg, region_lbl + total ~ reemplazo, value.var = "N", fill = 0L),
       file.path(TAB, "desc_salidas_reemplazo_region.csv"))

p_reemp <- ggplot(reemp_reg, aes(reorder(region_lbl, total), N, fill = reemplazo)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("Con reemplazo" = "#4575b4",
                               "Sin reemplazo" = "#d73027")) +
  labs(x = NULL, y = "Número de salidas", fill = NULL,
       title = "Salidas de estaciones por región, 2012–2026",
       subtitle = sprintf(paste("Con reemplazo: otra estación abre a menos de %d m dentro de %d meses.",
                               "Los reemplazos en el mismo sitio ya están fusionados. Excluye dic-2022"),
                          REEMP_M, REEMP_MESES),
       caption = fuente) +
  tema() + theme(panel.grid.major.y = element_blank())
guardar(p_reemp, "desc_salidas_reemplazo.pdf", 9, 6.5)

# ==============================================================================
# 7. the local market at the moment of entry
#
# for each entry, the stations already active that month within each radius,
# split by whether they carry the entrant's brand
# ==============================================================================

act <- panel[!is.na(lat), .(station_key, miym, lat, lon, distribuidor)]
ent <- merge(entradas[, .(station_key, g, edist, grupo)],
             unique(panel[, .(station_key, lat, lon)])[, .SD[1], by = station_key],
             by = "station_key")
ent[, gi := year(g) * 12L + (month(g) - 1L)]

setkey(act, miym)
vecinos <- rbindlist(lapply(seq_len(nrow(ent)), function(k) {
  a <- act[.(ent$gi[k]), nomatch = NULL][station_key != ent$station_key[k]]
  if (!nrow(a)) return(data.table(k = k, n25 = 0L, s25 = 0L, n50 = 0L, s50 = 0L))
  la <- a$lat * pi / 180; lo <- a$lon * pi / 180
  l0 <- ent$lat[k] * pi / 180; o0 <- ent$lon[k] * pi / 180
  hh <- sin((la - l0) / 2)^2 + cos(l0) * cos(la) * sin((lo - o0) / 2)^2
  hh[hh > 1] <- 1
  d <- 2 * R_EARTH * asin(sqrt(hh))
  misma <- a$distribuidor == ent$edist[k]
  data.table(k = k, n25 = sum(d <= RTREAT), s25 = sum(d <= RTREAT & misma),
             n50 = sum(d <= RCTRL), s50 = sum(d <= RCTRL & misma))
}))
ent <- cbind(ent, vecinos[, !"k"])
fwrite(ent[, .(station_key, g, edist, n25, s25, n50, s50)],
       file.path(TAB, "desc_mercado_local_detalle.csv"))

mercado <- data.table(
  Radio = c("2 km", "5 km"),
  `Entradas` = nrow(ent),
  `\\% con al menos un vecino` = round(100 * c(mean(ent$n25 > 0), mean(ent$n50 > 0)), 1),
  `Vecinos (media)` = round(c(mean(ent$n25), mean(ent$n50)), 2),
  `Vecinos (mediana)` = c(median(ent$n25), median(ent$n50)),
  `Vecinos (máx.)` = c(max(ent$n25), max(ent$n50)),
  `\\% misma marca` = round(100 * c(sum(ent$s25) / sum(ent$n25),
                                    sum(ent$s50) / sum(ent$n50)), 1))
fwrite(mercado, file.path(TAB, "desc_mercado_local.csv"))
guardar_tabla_tex(mercado[, lapply(.SD, as.character)], "tab_desc_mercado_local.tex")
cat("\n=== MERCADO LOCAL AL ENTRAR ===\n"); print(mercado)

# composition by the entrant's own brand
comp <- rbind(
  ent[, .(radio = "2 km", tot = sum(n25), same = sum(s25), n = .N), by = grupo],
  ent[, .(radio = "5 km",   tot = sum(n50), same = sum(s50), n = .N), by = grupo])
comp[, `:=`(`Misma marca` = 100 * same / tot,
            `Distinta marca` = 100 * (1 - same / tot))]
fwrite(comp, file.path(TAB, "desc_mercado_local_por_marca.csv"))

comp_l <- melt(comp, id.vars = c("grupo", "radio", "n"),
               measure.vars = c("Misma marca", "Distinta marca"),
               variable.name = "marca", value.name = "pct")
comp_l[, lbl := sprintf("%s (n=%d)", grupo, n)]
orden <- grupos_orden
comp_l[, lbl := factor(lbl, levels = unique(lbl[order(match(grupo, orden))]))]

p_merc <- ggplot(comp_l, aes(radio, pct, fill = marca)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_stack(vjust = 0.5), colour = "white", size = 3.1) +
  facet_wrap(~lbl, nrow = 1) +
  scale_fill_manual(values = c("Misma marca" = "#4575b4",
                               "Distinta marca" = "#d73027")) +
  labs(x = "Radio", y = "% de estaciones vecinas", fill = NULL,
       title = "Composición de marca del mercado local al momento de la entrada",
       subtitle = sprintf("Vecinos promedio: %.1f a 2 km y %.1f a 5 km; %.0f%% de las entradas tiene al menos un vecino a 2 km",
                          mean(ent$n25), mean(ent$n50), 100 * mean(ent$n25 > 0)),
       caption = fuente) +
  tema() + theme(panel.grid.major.x = element_blank())
guardar(p_merc, "desc_mercado_local.pdf", 11, 5)

message("18_descriptivas.R: figuras en ", FIG, " y tablas en ", TAB, ".")
