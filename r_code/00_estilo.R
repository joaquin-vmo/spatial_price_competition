
# objective: one place for the look of every figure in the thesis, so that the
# scripts do not each carry their own theme and palette.
#
# palatino is the typeface of the document. the default pdf() device neither
# embeds it nor handles the unicode symbols in the axis labels; quartz, the
# native macos engine, does both. outside macos it falls back to ggsave.
#
# the palette is blue and orange rather than the blue and red of fischer, martin
# & schmidt-dengler (2025): red reads as a warning and these series are not
# warnings. blue and orange is also the pair that survives red-green colour
# blindness and keeps its contrast in greyscale printing -- gold and yellow were
# tried first and proved too faint against white, particularly for error bars.

library(ggplot2)

TIPO <- "Palatino"

AZUL     <- "#0037b9"
NARANJO  <- "#e8710a"
GRIS     <- "#7f7f7f"
AZUL_CL  <- "#5b8db8"
NARANJO_OSC <- "#b8560a"

# sequential ramp for ordered categories (distance rings, terciles)
PAL_SEC <- c(AZUL, AZUL_CL, NARANJO, NARANJO_OSC, GRIS)

tema <- function(base = 14) {
  theme_minimal(base_size = base, base_family = TIPO) +
    theme(text = element_text(family = TIPO),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.caption = element_text(colour = "grey40", size = rel(0.75)))
}

guardar <- function(fig, ruta, width, height) {
  if (capabilities("aqua")) {
    grDevices::quartz(file = ruta, type = "pdf", width = width, height = height)
    print(fig)
    grDevices::dev.off()
  } else {
    ggsave(ruta, fig, width = width, height = height)
  }
}
