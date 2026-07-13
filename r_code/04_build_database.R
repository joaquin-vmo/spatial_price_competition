
# objective: merge the processed exports into one database

library(data.table)

cne_prices <- fread(file.path("data/processed", "cne_prices.csv"),
                     encoding = "UTF-8", colClasses = list(character = "id"))

cne_stations <- fread(file.path("data/processed", "cne_stations.csv"),
                       encoding = "UTF-8", colClasses = list(character = "codigo"))

mepco <- fread(file.path("data/processed", "whosale_prices_mepco.csv"))

# prices <- stations attributes, matched by station id (codigo)
db <- merge(cne_prices, cne_stations, by.x = "id", by.y = "codigo", all.x = TRUE)

#####

# add mepco wholesale prices: for each retail price, the wholesale price
# vigente at that date, i.e. the most recent mepco quote on or before it.
# mepco has no 95 quote -> those rows stay na.

fuels <- c("93", "97", "di")

mepco_long <- rbindlist(lapply(fuels, function(f) {
  data.table(
    fuel = f,
    date = mepco$date,
    wholesale_w_o = mepco[[paste0(f, "_w/o")]],
    wholesale_w = mepco[[paste0(f, "_w/")]],
    variable_specific_tax_utm_m3 = mepco[[paste0(f, "_variable_specific_tax_utm_m3")]]
  )
}))

# wholesale price change at each mepco update: the size of the shock a
# station is facing when its vigente wholesale price last moved
setorder(mepco_long, fuel, date)
mepco_long[, `:=`(
  d_wholesale_w_o = wholesale_w_o - shift(wholesale_w_o),
  d_wholesale_w   = wholesale_w   - shift(wholesale_w)
), by = fuel]

setkey(mepco_long, fuel, date)
db <- mepco_long[db, on = .(fuel, date), roll = TRUE]

#####

# first/last date each station (id) was observed in the data

db[, `:=`(
  first_date = min(date),
  last_date = max(date)
), by = id]

#####

# event-time relative to the march 2026 mepco shock (2026-03-26, +372 on 93 w/)

event_date <- as.Date("2026-03-26")
db[, `:=`(
  days_to_event = as.integer(date - event_date),
  post_shock = date >= event_date
)]

#####

# price rigidity: days since this station/fuel's price last changed, and
# whether the current row is itself a change point

setorder(db, id, fuel, date)
db[, price_spell_id := rleid(price), by = .(id, fuel)]
db[, `:=`(
  last_price_change_date = min(date),
  price_changed = date == min(date)
), by = .(id, fuel, price_spell_id)]
db[, days_since_last_price_change := as.integer(date - last_price_change_date)]
db[, `:=`(price_spell_id = NULL, last_price_change_date = NULL)]
# first observed spell per station has no known start -> change flag/duration undefined
db[date == first_date, `:=`(price_changed = NA, days_since_last_price_change = NA_integer_)]

#####

# spatial competition: n° of other stations within k km that were active
# (first_date <= date <= last_date) on that same date, based on each
# station's modal reported coordinate (a handful drift slightly across years)

stations_coords <- db[, .N, by = .(id, latitud, longitud)]
setorder(stations_coords, id, -N)
stations_coords <- stations_coords[, .SD[1], by = id][, .(id, latitud, longitud)]

lat_rad <- stations_coords$latitud * pi / 180
lon_rad <- stations_coords$longitud * pi / 180

R_earth <- 6371 # km
dlat <- outer(lat_rad, lat_rad, "-")
dlon <- outer(lon_rad, lon_rad, "-")
a <- sin(dlat / 2)^2 + outer(cos(lat_rad), cos(lat_rad)) * sin(dlon / 2)^2
a[a > 1] <- 1 # guard against floating-point overshoot past asin's domain
dist_km <- 2 * R_earth * asin(sqrt(a))
diag(dist_km) <- Inf # a station is not its own competitor

radii <- 1:5
max_dist_cutoff <- 10 # km; also bounds the continuous intensity measure below

# every (i, j) pair within the cutoff, with the activity window of j
pairs_idx <- which(dist_km <= max_dist_cutoff, arr.ind = TRUE)
pairs <- data.table(
  i_id = stations_coords$id[pairs_idx[, "row"]],
  j_id = stations_coords$id[pairs_idx[, "col"]],
  dist = dist_km[pairs_idx]
)
first_last <- unique(db[, .(id, first_date, last_date)])
pairs <- merge(pairs, first_last, by.x = "j_id", by.y = "id")

# for each radius, turn j's activity window into +1/-1 entry/exit events on
# i's timeline, then a running sum rolled onto i's actual report dates gives
# the count of active competitors on that date (analogous to the mepco join)
id_date <- unique(db[, .(id, date)])
setkey(id_date, id, date)

for (r in radii) {
  in_radius <- pairs[dist <= r]
  events <- rbind(
    in_radius[, .(id = i_id, event_date = first_date, delta = 1L)],
    in_radius[, .(id = i_id, event_date = last_date + 1, delta = -1L)]
  )
  events <- events[, .(delta = sum(delta)), by = .(id, event_date)]
  setorder(events, id, event_date)
  events[, cum_count := cumsum(delta), by = id]
  setkey(events, id, event_date)

  col_name <- paste0("competitors_", r, "km")
  vals <- events[id_date, on = .(id, event_date = date), roll = TRUE]$cum_count
  id_date[[col_name]] <- fifelse(is.na(vals), 0L, vals)
}

# continuous intensity: inverse-distance-weighted sum of active competitors
# within the cutoff, so a rival at 1km counts for more than one at 1.2km.
# distances below the floor (likely geocoding noise between near-duplicate
# points) are clamped so weights don't blow up.
dist_floor <- 0.1 # km
pairs[, weight := 1 / pmax(dist, dist_floor)]

events <- rbind(
  pairs[, .(id = i_id, event_date = first_date, delta = weight)],
  pairs[, .(id = i_id, event_date = last_date + 1, delta = -weight)]
)
events <- events[, .(delta = sum(delta)), by = .(id, event_date)]
setorder(events, id, event_date)
events[, cum_intensity := cumsum(delta), by = id]
setkey(events, id, event_date)

vals <- events[id_date, on = .(id, event_date = date), roll = TRUE]$cum_intensity
id_date$competition_intensity <- fifelse(is.na(vals), 0, vals)

db <- merge(db, id_date, by = c("id", "date"), all.x = TRUE)

rm(stations_coords, lat_rad, lon_rad, dlat, dlon, a, dist_km, R_earth,
   radii, max_dist_cutoff, pairs_idx, pairs, first_last, id_date, r,
   in_radius, events, col_name, vals, dist_floor)

#####

# brand structure: "sin bandera" is cne's own unbranded/independent label;
# among named distributors, a single-station name is an independent operator
# trading under their own name rather than a real multi-station network

distributor_n_stations <- db[, .(distributor_n_stations = uniqueN(id)), by = distributor]
db <- merge(db, distributor_n_stations, by = "distributor", all.x = TRUE)
db[, is_franchise := distributor != "sin bandera" & distributor_n_stations > 1]

rm(distributor_n_stations)

#####

setorder(db, id, date, fuel)

fwrite(db, file.path("data/processed", "database.csv"))
