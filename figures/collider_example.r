library(ggplot2)
library(dplyr)
library(this.path)

set.seed(20260428)
setwd(this.path::here())
source("fps_style.r")

n <- 1200

dat <- tibble(
  cad_prs = rnorm(n),
  response = rnorm(n)
) %>%
  mutate(
    # Deterministic inclusion rule:
    # inclusion if combined CAD PRS + LDL-C liability is high
    inclusion_score = cad_prs + response,
    included = inclusion_score > quantile(inclusion_score, 0.88),
    group = if_else(included,
      "Individuals included in the trial",
      "Individuals excluded from the trial"
    )
  )

p <- ggplot(dat, aes(x = cad_prs, y = response)) +
  geom_point(
    data = filter(dat, !included),
    aes(fill = group, colour = group),
    shape = 21,
    size = 3.0,
    stroke = 0.6,
    alpha = 0.45
  ) +
  geom_point(
    data = filter(dat, included),
    aes(fill = group, colour = group),
    shape = 21,
    size = 3.0,
    stroke = 0.7,
    alpha = 0.65
  ) +
  # Full population line of best fit
  geom_smooth(
    data = dat,
    aes(x = cad_prs, y = response, linetype = "Full population"),
    method = "lm",
    se = FALSE,
    colour = fps$navy,
    linewidth = 1.5
  ) +
  # Trial inclusion line of best fit
  geom_smooth(
    data = filter(dat, included),
    aes(x = cad_prs, y = response, linetype = "Trial participants"),
    method = "lm",
    se = FALSE,
    colour = fps$sand,
    linewidth = 1.5
  ) +
  scale_fill_manual(
    values = c(
      "Individuals excluded from the trial" = fps$steel,
      "Individuals included in the trial"   = fps$sand_mid
    ),
    labels = c(
      "Individuals excluded from the trial" = "Individuals excluded\nfrom the trial",
      "Individuals included in the trial"   = "Individuals included\nin the trial"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Individuals excluded from the trial" = fps$slate,
      "Individuals included in the trial"   = fps$sand
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Full population"   = "dashed",
      "Trial participants" = "dashed"
    ),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-4.3, 4.3), ylim = c(-4.3, 4.3)) +
  labs(
    x = "CAD polygenic risk",
    y = "Other risk factors",
    fill = NULL,
    colour = NULL
  ) +
  guides(
    colour = "none",
    fill = guide_legend(
      order = 1,
      override.aes = list(alpha = 0.75, size = 3.4, shape = 21)
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(
        colour    = c(fps$navy, fps$sand),
        linewidth = 1.5
      )
    )
  ) +
  theme_fps(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(colour = fps$steel, fill = fps$white)
  )

if (interactive()) print(p)
ggsave(
  filename = "collider_example.png",
  plot = p,
  width = 8.5,
  height = 5.5,
  dpi = 320,
  device = ragg::agg_png
)
