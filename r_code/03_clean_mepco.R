
# objective: load and clean mepco weekly data

library(readxl)
library(dplyr)

mepco <- read_excel("data/raw/mepco.xlsx",
                    col_types = c("date", "text", "text",
                                  "text", "text", "text", "text", "skip",
                                  "text", "text", "text"))

# mepco only covers 93, 97 and diesel (no 95); col 8 is an empty separator.
# w/o = precio sin mepco, w/ = precio con mepco
colnames(mepco) <- c("date",
                     "93_w/o", "93_w/",
                     "97_w/o", "97_w/",
                     "di_w/o", "di_w/",
                     "93_variable_specific_tax_utm_m3",
                     "97_variable_specific_tax_utm_m3",
                     "di_variable_specific_tax_utm_m3")

# drop the sub-header row ("Precio sin/con Mepco")
mepco <- mepco[-1, ]

# all numeric columns use '.' as decimal -> plain as.numeric
mepco <- mepco |>
  mutate(across(-date, as.numeric),
         date = as.Date(date))

fwrite(mepco, file.path("data/processed/", "whosale_prices_mepco.csv"))

