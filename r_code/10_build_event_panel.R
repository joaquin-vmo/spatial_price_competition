
# objective: build the balanced panels and the entry-event design used by the
# event studies. only ENTRY is treated as an event here: exits are not
# classified (the 2022/2023 spike in last_date is a legacy->current format
# artifact, and the 2026 mass is right-censoring), so nothing downstream
# depends on telling a real closure apart from a reporting gap.
#
# frequency: WEEKLY on the mepco thursday grid, plus a monthly collapse.
# a daily panel is not built on purpose: reporting cadence is ~weekly (median
# 50 obs per station-fuel-year, 99.8% of rows are price changes), so a daily
# grid would be ~85% carry-forward with no extra information, and mepco --the
# cost benchmark-- is itself weekly.
#
# TAKES:
#   data/processed/database.csv          (from 04_build_database.R)
#   data/processed/whosale_prices_mepco.csv
#
# PRODUCES (all in data/processed/):
#   panel_semanal.csv.gz  station_key x fuel x week, price + margin, locf-filled
#   panel_mensual.csv     station_key x month, wide prices/margins + treatment
#   entradas.csv          one row per entry event
#   stack_entrada.csv     stacked did, ring control (RTREAT-RCTRL km)
#   stack_entrada_far.csv stacked did, far control (never-treated, same region)

library(data.table)

# ------------------------------------------------------------------------------
# parameters
# ------------------------------------------------------------------------------
RTREAT   <- 2     # treatment radius (km): incumbent is treated if entry within
RCTRL    <- 5     # ring-control radius (km)
PRE      <- 12L   # event-study pre window (months)
POST     <- 12L   # event-study post window (months)
MAXGAP_W <- 13L   # max weeks of carry-forward (~3 months); beyond it, no price
WEEK0    <- as.IDate("2012-01-05")  # thursday anchor; mepco starts 2014-08-07,
                                    # exactly 135 weeks later, so it sits on the grid
FOCAL_FROM <- 2014L  # entries from this year are focal treatments. 2013 entries
                     # (81 vs ~30-50 in a normal year) are suspect as a coverage
                     # artifact at the start of the sample: they contaminate and
                     # cut, but are never used as treatment
LINK_M   <- 50    # metres: a new id succeeding another one this close is the
                  # same physical station (recoding), not an entry

R_EARTH <- 6371   # km

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------

modal <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mi <- function(d) year(d) * 12L + (month(d) - 1L)   # date -> continuous month

mi2date <- function(m) {
  out <- rep(NA_integer_, length(m))
  ok <- !is.na(m)
  out[ok] <- as.integer(as.IDate(sprintf("%d-%02d-01",
                                         m[ok] %/% 12L, m[ok] %% 12L + 1L)))
  as.IDate(out)
}

# great-circle distance matrix (km) between two sets of coordinates
dist_km <- function(lat1, lon1, lat2 = lat1, lon2 = lon1) {
  a1 <- lat1 * pi / 180; o1 <- lon1 * pi / 180
  a2 <- lat2 * pi / 180; o2 <- lon2 * pi / 180
  dlat <- outer(a1, a2, "-")
  dlon <- outer(o1, o2, "-")
  h <- sin(dlat / 2)^2 + outer(cos(a1), cos(a2)) * sin(dlon / 2)^2
  h[h > 1] <- 1
  2 * R_EARTH * asin(sqrt(h))
}

# ==============================================================================
# 1. load
# ==============================================================================

cols <- c("id", "date", "fuel", "price", "wholesale_w", "distributor",
          "municipality", "region", "latitud", "longitud", "is_franchise")
db <- fread(file.path("data/processed", "database.csv"), select = cols,
            encoding = "UTF-8", colClasses = list(character = "id"))
db[, date := as.IDate(date)]
setnames(db, c("municipality", "region", "distributor"),
         c("comuna", "region", "distribuidor"))

# ==============================================================================
# 2. physical station identity
#
# a new id whose first report comes at or after another id's last report, and
# whose coordinates are within LINK_M, is the same site recoded (typically a
# change of operator). counting it as an entry would fabricate treatment, so the
# ids are chained to a single station_key. this is defence #14 (cleaning the
# entry definition) done at the source.
# ==============================================================================

idinfo <- db[, .(first_date = min(date), last_date = max(date),
                 lat = median(latitud), lon = median(longitud)), by = id]
setorder(idinfo, first_date, id)
idinfo[, parent := NA_character_]

cand <- idinfo[!is.na(lat)]
D_link <- dist_km(cand$lat, cand$lon) * 1000  # metres
diag(D_link) <- Inf

for (i in seq_len(nrow(cand))) {
  # potential ancestors: strictly earlier start, already dead when i appears
  prev <- which(cand$last_date <= cand$first_date[i] &
                  cand$first_date < cand$first_date[i])
  if (!length(prev)) next
  d <- D_link[i, prev]
  j <- which.min(d)
  if (d[j] < LINK_M) idinfo[id == cand$id[i], parent := cand$id[prev[j]]]
}

# duplicate records of the same site. the legacy files carry a second id for
# some stations, formed by appending "a" to the base id -- 31 such pairs, every
# one of which has its base id present and coexists with it, and 25 of which
# stop being reported in december 2022 when the source consolidated them into
# the service_type column of the new format. they are the same physical station,
# so the suffixed id is folded into its base. the succession rule below cannot
# catch these because the two records COEXIST rather than follow one another
dup <- idinfo[grepl("[0-9]a$", id) & sub("a$", "", id) %in% idinfo$id]
idinfo[dup, on = "id", parent := sub("a$", "", id)]
message(sprintf("duplicados de sitio (sufijo 'a') fusionados: %d", nrow(dup)))

# follow the chain to the root id
par <- setNames(idinfo$parent, idinfo$id)
root_of <- function(x) {
  while (!is.na(par[[x]])) x <- par[[x]]
  x
}
idinfo[, station_key := vapply(id, root_of, character(1))]
message(sprintf("recodificaciones enlazadas: %d de %d ids",
                idinfo[!is.na(parent), .N], nrow(idinfo)))

db <- merge(db, idinfo[, .(id, station_key)], by = "id", all.x = TRUE)

# ==============================================================================
# 3. weekly panel: last price of each thursday-week, locf within MAXGAP_W
# ==============================================================================

db[, wi := as.integer(date - WEEK0) %/% 7L]        # continuous week index
db <- db[wi >= 0L]

setorder(db, station_key, fuel, date)
obs_w <- db[db[, .I[.N], by = .(station_key, fuel, wi)]$V1,
            .(station_key, fuel, wi, price, wholesale_w)]

grid_w <- obs_w[, .(wi = seq.int(min(wi), max(wi))), by = .(station_key, fuel)]
grid_w <- merge(grid_w, obs_w[, .(station_key, fuel, wi, price, wholesale_w,
                                  obs_wi = wi)],
                by = c("station_key", "fuel", "wi"), all.x = TRUE)
grid_w[, observed := !is.na(price)]
setorder(grid_w, station_key, fuel, wi)
grid_w[, `:=`(price  = nafill(price,  "locf"),
              obs_wi = nafill(obs_wi, "locf")),
       by = .(station_key, fuel)]
# a gap longer than MAXGAP_W is a reporting hole (or an inactive site): the
# price is not carried across it and the station is simply absent those weeks
grid_w[wi - obs_wi > MAXGAP_W, price := NA_real_]
grid_w <- grid_w[!is.na(price)]
grid_w[, obs_wi := NULL]

# wholesale cost: mepco is national and weekly, so it is attached by (fuel, week)
# rather than carried forward per station. it is the vigente quote of that week
mepco <- fread(file.path("data/processed", "whosale_prices_mepco.csv"))
mepco[, date := as.IDate(date)]
mepco_long <- melt(
  mepco[, .(date, `93` = `93_w/`, `97` = `97_w/`, di = `di_w/`)],
  id.vars = "date", variable.name = "fuel", value.name = "cost",
  variable.factor = FALSE
)
mepco_long[, wi := as.integer(date - WEEK0) %/% 7L][, date := NULL]

grid_w[, wholesale_w := NULL]
grid_w <- merge(grid_w, mepco_long, by = c("fuel", "wi"), all.x = TRUE)
grid_w[, margin := price - cost]        # NA for fuel 95: mepco publishes no 95
grid_w[, wk := WEEK0 + wi * 7L]
grid_w[, ym := mi2date(mi(wk))]

# station attributes, one row per station_key
sattr <- db[, .(id = modal(id), distribuidor = modal(distribuidor),
                comuna = modal(comuna), region = modal(region),
                lat = median(latitud, na.rm = TRUE),
                lon = median(longitud, na.rm = TRUE),
                is_franchise = as.logical(modal(as.character(is_franchise)))),
            by = station_key]

panel_w <- merge(grid_w, sattr, by = "station_key", all.x = TRUE)
setcolorder(panel_w, c("station_key", "fuel", "wi", "wk", "ym"))
setorder(panel_w, station_key, fuel, wi)

# ==============================================================================
# 4. monthly panel (wide): last week of each month
# ==============================================================================

setorder(panel_w, station_key, fuel, wi)
last_m <- panel_w[panel_w[, .I[.N], by = .(station_key, fuel, ym)]$V1,
                  .(station_key, fuel, ym, price, margin,
                    observed_m = observed)]

fuel_lbl <- c(`93` = "93", `95` = "95", `97` = "97", di = "di")
px  <- dcast(last_m, station_key + ym ~ fuel, value.var = "price")
setnames(px,  names(fuel_lbl), paste0("p", fuel_lbl), skip_absent = TRUE)
mg  <- dcast(last_m, station_key + ym ~ fuel, value.var = "margin")
setnames(mg,  names(fuel_lbl), paste0("m", fuel_lbl), skip_absent = TRUE)
mg[, m95 := NULL]                       # no mepco reference for 95
obsm <- last_m[, .(observed = any(observed_m)), by = .(station_key, ym)]

panel_m <- Reduce(function(a, b) merge(a, b, by = c("station_key", "ym")),
                  list(px, mg, obsm))
panel_m <- merge(panel_m, sattr, by = "station_key", all.x = TRUE)
panel_m[, miym := mi(ym)]
setorder(panel_m, station_key, ym)

# ==============================================================================
# 5. entry events
#
# an entry is the first month a physical station is ever observed. the 2012
# cohort is left-censored (1.548 of 2.086 stations) and is never an entry.
# stations that go dark for more than MAXGAP_W weeks and come back are not
# re-entries here either: only the very first appearance counts, so a reporting
# hole cannot manufacture an event.
# ==============================================================================

st_span <- panel_m[, .(entry_ym = min(ym), last_ym = max(ym)), by = station_key]
st_span <- merge(st_span, sattr, by = "station_key")
st_span[, base_2012 := year(entry_ym) == 2012L]

entradas_all <- st_span[base_2012 == FALSE & !is.na(lat),
                        .(station_key, g = entry_ym, elat = lat, elon = lon,
                          eregion = region, edist = distribuidor,
                          efranchise = is_franchise)]
setorder(entradas_all, g, station_key)
entradas_all[, event_id := paste0("EN", .I)]
entradas <- entradas_all[year(g) >= FOCAL_FROM]     # focal events only

message(sprintf("entradas: %d totales (>=2013), %d focales (>=%d)",
                nrow(entradas_all), nrow(entradas), FOCAL_FROM))

# ==============================================================================
# 6. treatment assignment at station level
#
# for every station: the months of all entries within RTREAT km that happen
# while it is already an incumbent. g = first such entry (its cohort), g2 =
# second (observations from g2 on are dropped, so the estimate is the effect of
# a first entry, as in fischer et al.). role: treated / ring / far.
# ==============================================================================

sloc <- st_span[!is.na(lat), .(station_key, slat = lat, slon = lon,
                               region, entry_ym, last_ym)]

D_ev <- dist_km(sloc$slat, sloc$slon, entradas_all$elat, entradas_all$elon)
# a station is never treated by its own entry
own <- match(entradas_all$station_key, sloc$station_key)
D_ev[cbind(own, seq_len(nrow(entradas_all)))] <- Inf
# and only by entries that occur while it is already active
inc <- outer(mi(sloc$entry_ym), mi(entradas_all$g), "<") &
  outer(mi(sloc$last_ym), mi(entradas_all$g), ">=")
D_ev[!inc] <- Inf

g_ev <- mi(entradas_all$g)

near_months <- function(dmax) {
  lapply(seq_len(nrow(sloc)), function(i) {
    w <- which(D_ev[i, ] <= dmax)
    if (!length(w)) NULL else sort(unique(g_ev[w]))
  })
}
near_treat <- near_months(RTREAT)

asg <- data.table(
  station_key = sloc$station_key,
  mindist = round(apply(D_ev, 1, min), 3),
  g_mi  = vapply(near_treat, function(x) if (is.null(x)) NA_integer_ else x[1],
                 integer(1)),
  g2_mi = vapply(near_treat,
                 function(x) if (length(x) >= 2L) x[2] else NA_integer_,
                 integer(1))
)
asg[, `:=`(g_entry = mi2date(g_mi), g2_entry = mi2date(g2_mi))]
asg[, role_entry := fifelse(!is.na(g_mi), "treated",
                            fifelse(mindist <= RCTRL, "ring", "far"))]
# distance to the entry that treats the station (its own cohort's entrant)
asg[, dist_entry := NA_real_]
for (i in which(!is.na(asg$g_mi))) {
  w <- which(D_ev[i, ] <= RTREAT & g_ev == asg$g_mi[i])
  asg$dist_entry[i] <- round(min(D_ev[i, w]), 3)
}
asg[, c("g_mi", "g2_mi") := NULL]

# ------------------------------------------------------------------------------
# 6b. distance-ring assignment, for the attenuation test
#
# a separate assignment with a wider net: a station counts as treated by ANY
# entry within RCTRL km and is placed in the ring of its first one, so the
# effect can be traced from 0-1 km outwards. control = stations never within
# RCTRL of any entry. observations from a second entry within RCTRL are cut
# ------------------------------------------------------------------------------

ring_brk <- c(0, 1, 2, 3, 4, 5)
ring_lbl <- c("0-1 km", "1-2 km", "2-3 km", "3-4 km", "4-5 km")

asg5 <- rbindlist(lapply(seq_len(nrow(sloc)), function(i) {
  w <- which(D_ev[i, ] <= RCTRL)
  if (!length(w)) {
    return(data.table(station_key = sloc$station_key[i], g5_mi = NA_integer_,
                      g2_5_mi = NA_integer_, dist5 = NA_real_))
  }
  # earliest entry first; ties within a month broken by proximity
  w <- w[order(g_ev[w], D_ev[i, w])]
  ms <- unique(g_ev[w])
  data.table(station_key = sloc$station_key[i], g5_mi = ms[1],
             g2_5_mi = if (length(ms) >= 2L) ms[2] else NA_integer_,
             dist5 = D_ev[i, w[1]])
}))
asg5[, `:=`(g5_entry = mi2date(g5_mi), g2_5_entry = mi2date(g2_5_mi),
            dist_entry5 = round(dist5, 3),
            ring_entry = cut(dist5, breaks = ring_brk, labels = ring_lbl,
                             include.lowest = TRUE),
            role5_entry = fifelse(!is.na(dist5), "treated", "far"))]
asg5[, c("g5_mi", "g2_5_mi", "dist5") := NULL]

asg <- merge(asg, asg5, by = "station_key")

panel_m <- merge(panel_m, asg, by = "station_key", all.x = TRUE)
panel_m <- merge(panel_m, st_span[, .(station_key, base_2012)],
                 by = "station_key", all.x = TRUE)
panel_w <- merge(panel_w, asg, by = "station_key", all.x = TRUE)
panel_w <- merge(panel_w, st_span[, .(station_key, base_2012)],
                 by = "station_key", all.x = TRUE)

message(sprintf("asignacion: %d treated | %d ring | %d far",
                asg[role_entry == "treated", .N], asg[role_entry == "ring", .N],
                asg[role_entry == "far", .N]))

# ==============================================================================
# 7. stacked did design
#
# one clean experiment per focal entry: treated incumbents inside RTREAT,
# controls either in the RTREAT-RCTRL ring or never-treated in the same region.
# clean-control rules follow fischer et al.:
#   - a treated station enters only if the focal entry is its FIRST nearby
#     entry, and its observations are cut before a second one;
#   - a ring control enters only if it sees no nearby entry inside the window.
# ==============================================================================

pcols <- c("station_key", "ym", "miym", "p93", "p95", "p97", "pdi",
           "m93", "m97", "mdi", "comuna", "region", "distribuidor",
           "is_franchise", "observed")
pm <- panel_m[, ..pcols]
setkey(pm, station_key)

names(near_treat) <- sloc$station_key
first_mi <- setNames(mi(sloc$entry_ym), sloc$station_key)

build_stack <- function(ev, control) {
  far_ids <- asg[role_entry == "far", station_key]
  out <- vector("list", nrow(ev))
  for (i in seq_len(nrow(ev))) {
    e  <- ev[i]
    gi <- mi(e$g)
    lo <- gi - PRE
    hi <- gi + POST

    d <- D_ev[, match(e$event_id, entradas_all$event_id)]
    # treated: inside RTREAT, focal entry is their first, cut before the second
    tw <- which(d <= RTREAT)
    tr <- NULL
    if (length(tw)) {
      keys <- sloc$station_key[tw]
      cut <- rep(NA_integer_, length(keys))
      for (k in seq_along(keys)) {
        nb <- near_treat[[keys[k]]]
        nb <- nb[nb > first_mi[[keys[k]]]]
        if (!length(nb) || nb[1] != gi) next       # focal must be the first
        cut[k] <- if (length(nb) >= 2L) nb[2] - gi - 1L else POST
      }
      ok <- !is.na(cut)
      if (any(ok)) tr <- data.table(station_key = keys[ok],
                                    dist = d[tw][ok], role = "treated",
                                    cutv = pmin(cut[ok], POST))
    }

    # controls
    if (control == "ring") {
      cw <- which(d > RTREAT & d <= RCTRL)
      ct <- NULL
      if (length(cw)) {
        # clean ring control: no entry within RTREAT km inside the window. the
        # check is at RTREAT, not RCTRL: a ring station sits 2-5 km from the
        # focal entry, so testing at RCTRL would disqualify every one of them
        # on account of the focal entry itself
        keys <- sloc$station_key[cw]
        clean <- vapply(keys, function(s) {
          nb <- near_treat[[s]]
          is.null(nb) || !any(nb >= lo & nb <= hi)
        }, logical(1))
        if (any(clean)) ct <- data.table(station_key = keys[clean],
                                         dist = d[cw][clean],
                                         role = "control_ring", cutv = POST)
      }
    } else {
      cw <- which(sloc$station_key %in% far_ids & sloc$region == e$eregion &
                    mi(sloc$entry_ym) < gi & mi(sloc$last_ym) >= gi)
      ct <- if (length(cw))
        data.table(station_key = sloc$station_key[cw], dist = d[cw],
                   role = "control_far", cutv = POST) else NULL
    }

    sel <- rbind(tr, ct)
    if (is.null(sel) || !nrow(sel)) next

    pk <- pm[.(sel$station_key), nomatch = NULL][sel, on = "station_key"]
    pk[, event_time := miym - gi]
    pk <- pk[event_time >= -PRE & event_time <= cutv]
    if (!nrow(pk)) next

    out[[i]] <- pk[, .(
      event_id = e$event_id, g = e$g, src_key = e$station_key,
      src_distribuidor = e$edist, src_franchise = e$efranchise,
      station_key, role, treated = as.integer(role == "treated"),
      dist_km = round(dist, 3), ym, event_time,
      post = as.integer(event_time >= 0),
      p93, p95, p97, pdi, m93, m97, mdi,
      comuna, region, distribuidor, is_franchise, observed
    )]
  }
  rbindlist(out)
}

stack_entrada     <- build_stack(entradas, "ring")
stack_entrada_far <- build_stack(entradas, "far")

# ==============================================================================
# 8. write
# ==============================================================================

fwrite(panel_w,           "data/processed/panel_semanal.csv.gz")
fwrite(panel_m,           "data/processed/panel_mensual.csv")
fwrite(entradas_all,      "data/processed/entradas.csv")
fwrite(stack_entrada,     "data/processed/stack_entrada.csv")
fwrite(stack_entrada_far, "data/processed/stack_entrada_far.csv")

# ------------------------------------------------------------------------------
# summary
# ------------------------------------------------------------------------------

cat("\n=== PANEL SEMANAL ===\n")
cat("filas:", nrow(panel_w), "| estaciones:", uniqueN(panel_w$station_key),
    "| semanas:", uniqueN(panel_w$wi), "\n")
cat("rango:", format(min(panel_w$wk)), "a", format(max(panel_w$wk)), "\n")
cat("% semanas observadas (resto arrastrado):",
    round(100 * mean(panel_w$observed), 1), "\n")

cat("\n=== PANEL MENSUAL ===\n")
cat("filas:", nrow(panel_m), "| estaciones:", uniqueN(panel_m$station_key),
    "| meses:", uniqueN(panel_m$ym), "\n")

cat("\n=== ENTRADAS por anio ===\n")
print(entradas_all[, .N, by = .(anio = year(g))][order(anio)])

cat("\n=== ASIGNACION (estaciones) ===\n")
print(asg[, .N, by = role_entry][order(-N)])
cat("distancia a la entrada tratante (km), tratadas:\n")
print(round(quantile(asg[role_entry == "treated", dist_entry],
                     probs = c(0, .25, .5, .75, 1), na.rm = TRUE), 2))
# coincident coordinates: two station_keys that coexist at the same point are
# either a genuine adjacent pair or a duplicate record the LINK_M rule cannot
# catch (it only links successions, never coexistences). flagged, not dropped
cat("tratadas con entrante a <50 m (coordenadas coincidentes):",
    asg[role_entry == "treated" & dist_entry < 0.05, .N], "\n")

cat("\n=== ASIGNACION POR ANILLOS (radio", RCTRL, "km) ===\n")
print(asg[, .N, by = .(role5_entry, ring_entry)][order(role5_entry, ring_entry)])

cat("\n=== STACKED DiD ===\n")
for (nm in c("stack_entrada", "stack_entrada_far")) {
  d <- get(nm)
  cat(nm, ": filas", nrow(d), "| eventos", uniqueN(d$event_id),
      "| pares evento-estacion", uniqueN(paste(d$event_id, d$station_key)), "\n")
  print(d[, .(pares = uniqueN(paste(event_id, station_key))), by = role])
}

message("10_build_event_panel.R: paneles y diseno de entrada escritos en data/processed/")
