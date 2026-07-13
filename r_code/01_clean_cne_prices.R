
# objective: load and unify cne prices data into one data frame

# loading
library(data.table)

# cleaning
library(dplyr)
library(stringr)

for (año in 2012:2026) {
  assign(
    paste0("raw_", año),
    fread(
      file.path("data/raw/cne_prices", paste0(año, ".csv")),
      encoding = "UTF-8",
      colClasses = "character"
    )
  )
}


#####

# legacy datase loading and cleaning

raw_legacy <- data.table::rbindlist(mget(paste0("raw_", 2012:2022)))

colnames(raw_legacy) <- c("id", "legal_name", "distributor",
                          "address_street", "address_number",
                          "municipality", "region", "price",
                          "date", "fuel", "latitud", "longitud")

setorder(raw_legacy, id, date)
rm(list = paste0("raw_", 2012:2022))

# some data cleaning

clean_legacy <- raw_legacy |>
  select(-address_street, -address_number, -legal_name)

# id cleaning

clean_legacy <- clean_legacy |>
  filter(!id %in% c("prueba", "CNE 01")) |>
  mutate(
    id = if_else(id == "co730401 co730401", "co730401", id),
    id = tolower(trimws(id))
  )

# distributor cleaning

clean_legacy <- clean_legacy |>
  mutate(
    distributor = tolower(trimws(distributor)),
    distributor = distributor |>
      replace_values(
        "vigu ltda." ~ "vigu",
        "ruta 45" ~ "ruta v45",
        "petrobras" ~ "aramco_petrobras"
      )
  )

# municipality cleaning

clean_legacy <- clean_legacy |>
  mutate(
    municipality = tolower(trimws(municipality)),
    municipality = municipality |>
      replace_values(
        "aisén" ~ "aysen"
      ),
    municipality = chartr("áéíóúüñ", "aeiouun", municipality),
    municipality = str_replace_all(municipality, "['’]", ""),
    municipality = str_replace_all(municipality, "\\s+", "_")
  )

# region cleaning

clean_legacy <- clean_legacy |>
  mutate(
    region = str_replace_all(region, "\\p{Cf}", ""),
    region = tolower(trimws(region)),
    region = region |>
      replace_values(
        "bío bío" ~ "biobio",
        "biobío" ~ "biobio",
        "tarapacá" ~ "tarapaca",
        "los ríos" ~ "los_rios",
        "ñuble" ~ "nuble",
        "valparaíso" ~ "valparaiso",
        "los lagos" ~ "los_lagos",
        "gral. bernardo o'higgins" ~ "ohiggins",
        "araucanía" ~ "araucania",
        "magallanes y la antártida chilena" ~ "magallanes_antartida",
        "aysén gral. c. ibáñez del campo" ~ "aysen",
        "arica y parinacota" ~ "arica"
      )
  )

# fuel and prices cleaning

setDT(clean_legacy)

comb_map <- c(
  "Gasolina 93"     = "93",
  "Gasolina 95"     = "95",
  "Gasolina 97"     = "97",
  "Petroleo Diesel" = "di"
)
comb_interes <- names(comb_map)

# price to numeric; 0 and 99999999 sentinel -> na
clean_legacy[, price := as.numeric(trimws(price))]
clean_legacy[price == 0 | price == 99999999, price := NA_real_]

# temp year for annual medians
clean_legacy[, año := as.integer(substr(date, 1, 4))]
clean_legacy[fuel %in% comb_interes, comb_can := comb_map[fuel]]

med_anual <- clean_legacy[
  fuel %in% comb_interes & price >= 200 & price <= 5000,
  .(med_anual = median(price, na.rm = TRUE)),
  by = .(comb_can, año)
]
clean_legacy <- merge(clean_legacy, med_anual, by = c("comb_can", "año"), all.x = TRUE)

# fix extra-zero errors: /10 or /100 if within +-35% of median, else na
clean_legacy[
  fuel %in% comb_interes & !is.na(price) &
    price > 2 * med_anual,
  price := {
    p10  <- price / 10
    p100 <- price / 100
    d10  <- abs(p10  - med_anual) / med_anual
    d100 <- abs(p100 - med_anual) / med_anual
    fifelse(
      d10 <= 0.35, p10,
      fifelse(d100 <= 0.35, p100, NA_real_)
    )
  }
]

# implausibly cheap -> na
clean_legacy[
  fuel %in% comb_interes & !is.na(price) & price < 200,
  price := NA_real_
]

clean_legacy[, c("comb_can", "año", "med_anual") := NULL]

# keep only fuels of interest with canonical names, drop na prices
clean_legacy <- clean_legacy[fuel %in% comb_interes & !is.na(price)]
clean_legacy[, fuel := comb_map[fuel]]

rm(med_anual)

# location cleaning
clean_legacy[, latitud := as.numeric(gsub(",", ".", latitud))]
clean_legacy[, longitud := as.numeric(gsub(",", ".", longitud))]
clean_legacy[latitud > 0 | latitud < -90, latitud := NA_real_]
clean_legacy[longitud > 0 | longitud < -180, longitud := NA_real_]

clean_legacy <- clean_legacy[!is.na(latitud) & !is.na(longitud)]

# date correction
clean_legacy[, date := as.Date(date)]

#####

# new data structure cleaning

raw_current <- data.table::rbindlist(mget(paste0("raw_", 2023:2026)))

rm(list = paste0("raw_", 2023:2026))

colnames(raw_current) <- c("id", "legal_name", "distributor",
                           "address", "latitud", "longitud",
                           "municipality", "region", "fuel",
                           "price", "pricing_unit", "service_type",
                           "date", "time", "ev_station", "gas_station")

setorder(raw_current, id, date)

# some data cleaning

clean_current <- raw_current |>
  select(-legal_name, -address)

# id cleaning

clean_current <- clean_current |>
  mutate(
    id = if_else(id == "co730401 co730401", "co730401", id),
    id = tolower(trimws(id))
  )

# distributor cleaning

clean_current <- clean_current |>
  mutate(
    distributor = tolower(trimws(distributor)),
    distributor = distributor |>
      replace_values(
        "go abastible" ~ "abastible",
        "gasco autogas" ~ "autogasco",
        "hola" ~ "hola!",
        "esa" ~ "sesa",
        "petrobras" ~ "aramco_petrobras",
        "aramco" ~ "aramco_petrobras"
      )
  )

# location cleaning

clean_current <- clean_current |>
  mutate(
    latitud = as.numeric(gsub(",", ".", latitud)),
    longitud = as.numeric(gsub(",", ".", longitud)),
    latitud = if_else(latitud > 0 | latitud < -90, NA_real_, latitud),
    longitud = if_else(longitud > 0 | longitud < -180, NA_real_, longitud)
  ) |>
  filter(!is.na(latitud) & !is.na(longitud))

# one coordinate per station: keep the modal one (a few sites re-geocoded in 2026)
clean_current <- clean_current |>
  add_count(id, latitud, longitud, name = "coord_n") |>
  mutate(
    latitud = latitud[which.max(coord_n)],
    longitud = longitud[which.max(coord_n)],
    .by = id
  ) |>
  select(-coord_n)

# municipality cleaning

clean_current <- clean_current |>
  mutate(
    municipality = tolower(trimws(municipality)),
    municipality = municipality |>
      replace_values(
        "santiago centro" ~ "santiago"
      ),
    municipality = chartr("áéíóúüñ", "aeiouun", municipality),
    municipality = str_replace_all(municipality, "['’]", ""),
    municipality = str_replace_all(municipality, "\\s+", "_")
  )

# region cleaning

clean_current <- clean_current |>
  mutate(
    region = str_replace_all(region, "\\p{Cf}", ""),
    region = tolower(trimws(region)),
    region = if_else(str_detect(region, "libertador"), "ohiggins", region),
    region = region |>
      replace_values(
        "metropolitana de santiago" ~ "metropolitana",
        "valparaíso" ~ "valparaiso",
        "del biobío" ~ "biobio",
        "del maule" ~ "maule",
        "de la araucanía" ~ "araucania",
        "de los lagos" ~ "los_lagos",
        "ñuble" ~ "nuble",
        "de los ríos" ~ "los_rios",
        "tarapacá" ~ "tarapaca",
        "magallanes y de la antártica chilena" ~ "magallanes_antartida",
        "aysén del gral. carlos ibáñez del campo" ~ "aysen",
        "arica y parinacota" ~ "arica"
      )
  )

# fuel cleaning

setDT(clean_current)

comb_map <- c(
  "93" = "93", "A93" = "93",
  "95" = "95", "A95" = "95",
  "97" = "97", "A97" = "97",
  "DI" = "di", "ADI" = "di"
)
comb_interes <- names(comb_map)

# keep fuels of interest, merge modalities 
# into one grade ("A" prefix = autoservicio)
clean_current <- clean_current[fuel %in% comb_interes]
clean_current[, fuel := comb_map[fuel]]

# service_type: full_service = asistido, self_service = autoservicio
clean_current[, service_type := fifelse(service_type == "Autoservicio", "self_service", "full_service")]

# after filtering, pricing_unit is always $/L -> drop it
clean_current[, pricing_unit := NULL]

# price cleaning

# price to numeric; 0 and 99999999 sentinel -> na
clean_current[, price := as.numeric(trimws(price))]
clean_current[price == 0 | price == 99999999, price := NA_real_]

# temp year for annual medians
clean_current[, año := as.integer(substr(date, 1, 4))]

med_anual <- clean_current[
  price >= 200 & price <= 5000,
  .(med_anual = median(price, na.rm = TRUE)),
  by = .(fuel, año)
]
clean_current <- merge(clean_current, med_anual, by = c("fuel", "año"), all.x = TRUE)

# fix extra-zero errors: /10 or /100 if within +-35% of median, else na
clean_current[
  !is.na(price) & price > 2 * med_anual,
  price := {
    p10  <- price / 10
    p100 <- price / 100
    d10  <- abs(p10  - med_anual) / med_anual
    d100 <- abs(p100 - med_anual) / med_anual
    fifelse(
      d10 <= 0.35, p10,
      fifelse(d100 <= 0.35, p100, NA_real_)
    )
  }
]

# implausibly cheap -> na
clean_current[!is.na(price) & price < 200, price := NA_real_]

clean_current[, c("año", "med_anual") := NULL]

# drop na prices
clean_current <- clean_current[!is.na(price)]

# date and time correction
clean_current[, date := as.Date(date)]
clean_current[, time := as.ITime(time, format = "%H:%M:%S")]

rm(med_anual, año, comb_interes, comb_map)

#####

# unify: legacy for the past (2012-2022), current from 2023 onwards
# legacy has no time column -> time stays NA for those rows (fill = TRUE)

clean_current <- clean_current[date >= as.Date("2023-01-01")]

clean_prices <- rbind(clean_legacy, clean_current, fill = TRUE)
setorder(clean_prices, id, date, fuel)

fwrite(clean_prices, file.path("data/processed/", "cne_prices.csv"))
