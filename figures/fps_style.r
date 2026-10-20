# 5PS figure style for this repo's R figure scripts
# (R port of the Python style copied from ASHG2026/old/python_implementation/pc_data.py).
#
# Usage: source(file.path(this.path::here(), "fps_style.r"))
# Render with ragg (ggsave(..., device = ragg::agg_png)) so the system font is picked up.

# Keep only font families the system can find, preserving preference order
.fps_installed <- function(families) {
  have <- unique(systemfonts::system_fonts()$family)
  found <- families[families %in% have]
  if (length(found) == 0) "DejaVu Sans" else found
}

# 5PS design system (9Sep2026 5PS Design System, tokens.css): chart-1..5, chart-null, grid, text
fps <- list(
  navy      = "#264454",
  steel     = "#628096",
  slate     = "#3d6073",
  sand      = "#e4b386",
  sand_mid  = "#edcaa8",
  null      = "#999999",
  grid      = grDevices::rgb(98, 128, 150, alpha = 0.25 * 255, maxColorValue = 255),
  charcoal  = "#2a2a2a",
  white     = "#ffffff",
  off_white = "#f9f8f6"
)

# 5PS core face is Raleway; until real TTFs are installed, the closest available fallback is Source Sans Pro
fps_fonts <- .fps_installed(c("Raleway", "Source Sans Pro", "DejaVu Sans"))
fps_font <- fps_fonts[1]

# Axis label, subscript 10, drawn in the figure font
neg_log10_p <- quote(-log[10](p))

# ggplot2 theme built from the tokens above
theme_fps <- function(base_size = 14, base_family = fps_font) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = fps$charcoal),
      axis.text = ggplot2::element_text(colour = fps$charcoal),
      axis.line = ggplot2::element_line(colour = fps$charcoal),
      axis.ticks = ggplot2::element_line(colour = fps$charcoal),
      panel.grid.major.y = ggplot2::element_line(colour = fps$grid, linewidth = 0.35),
      plot.background = ggplot2::element_rect(fill = fps$white, colour = NA),
      panel.background = ggplot2::element_rect(fill = fps$white, colour = NA),
      legend.background = ggplot2::element_rect(fill = fps$white, colour = NA),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = fps$charcoal, face = "bold")
    )
}

# Draw annotation text in the figure font and text colour by default
ggplot2::update_geom_defaults("text", list(family = fps_font, colour = fps$charcoal))
ggplot2::update_geom_defaults("label", list(family = fps_font, colour = fps$charcoal))
