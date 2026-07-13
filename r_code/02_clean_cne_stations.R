
# objective: rebuild the station attributes table from a corrupted xlsx

library(readxl)
library(data.table)
library(stringr)
library(dplyr)

# station_data.xlsx is a ';'-delimited, comma-decimal csv that got mangled by an
# excel round-trip: excel split each row at the decimal commas and merged some
# records via embedded newlines. rebuild it: rejoin cells with ",", split on
# every newline, keep the well-formed 37-field records, then parse by ";".

raw <- read_excel("data/raw/station_data.xlsx", sheet = "in",
                  col_names = FALSE, col_types = "text")

lines <- apply(raw, 1, function(r) paste(r[!is.na(r)], collapse = ","))
lines <- unlist(str_split(paste(lines, collapse = "\n"), "\n"))
lines <- str_trim(lines)

keep <- str_count(lines, ";") + 1 == 37

station_data <- fread(text = paste(lines[keep], collapse = "\n"),
                      sep = ";", header = TRUE, quote = "",
                      colClasses = "character", encoding = "UTF-8")

# coordinates: strip stray quotes, comma decimal -> point
station_data[, latitud  := as.numeric(gsub(",", ".", gsub('"', "", latitud)))]
station_data[, longitud := as.numeric(gsub(",", ".", gsub('"', "", longitud)))]

# cleaning unnecessary columns

station_data <- station_data |>
  select(-razon_social, -logo, -marca,
         -combustible, -precio, -unidad_cobro,
         -hora_actualizacion, -tipo_atencion, -horario_atencion,
         -cheque, -tarjeta_grandes_tiendas, -tarjeta_bancaria,
         -direccion, -nombre_region, -nombre_comuna, -latitud, -longitud) |>
  unique()

# keep the last update per station (attributes are constant within a station)

station_data <- station_data |>
  mutate(fecha_actualizacion = as.Date(fecha_actualizacion)) |>
  group_by(codigo) |>
  slice_max(fecha_actualizacion, n = 1, with_ties = FALSE) |>
  ungroup()

fwrite(station_data, file.path("data/processed/", "cne_stations.csv"))

rm(lines, keep, raw)
