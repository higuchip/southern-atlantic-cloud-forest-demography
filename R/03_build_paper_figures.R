suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

required_objects <- c("pca_result", "eigenvalues", "area_final", "tabela_final")
if (!all(vapply(required_objects, exists, logical(1), inherits = TRUE))) {
  stop("Run R/01_run_analysis.R before building the figures, or run R/run_all.R.")
}

area_colors <- c(
  `1` = "#6A3D9A", `2` = "#1F78B4",
  `3` = "#33A02C", `4` = "#E69F00"
)
area_shapes <- c(`1` = 16, `2` = 17, `3` = 15, `4` = 8)

theme_paper <- function(base_size = 9) {
  theme_minimal(base_family = "Arial", base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, color = "#E6E6E6"),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      legend.position = "bottom",
      plot.title = element_blank()
    )
}

save_paper_plot <- function(plot, stem, width_mm, height_mm) {
  png_path <- file.path("results/figures", paste0(stem, ".png"))
  tiff_path <- file.path("results/figures", paste0(stem, ".tiff"))
  ggsave(
    png_path, plot, device = ragg::agg_png,
    width = width_mm, height = height_mm, units = "mm",
    dpi = 600, bg = "white"
  )
  ggsave(
    tiff_path, plot, device = "tiff",
    width = width_mm, height = height_mm, units = "mm",
    dpi = 600, compression = "lzw", bg = "white"
  )
}

# Figure 2: PCA variance, environmental loadings, and plot scores.
variance <- eigenvalues
loadings <- as.data.frame(pca_result$var$coord[, 1:2, drop = FALSE])
scores <- data.frame(
  Area = as.integer(area_final),
  Dim.1 = pca_result$ind$coord[, 1],
  Dim.2 = pca_result$ind$coord[, 2]
)

scree <- data.frame(
  component = seq_len(min(10, nrow(variance))),
  explained = variance$variance.percent[seq_len(min(10, nrow(variance)))]
)

p2a <- ggplot(scree, aes(component, explained)) +
  geom_col(width = 0.72, fill = "#4C78A8", color = "#1D3557", linewidth = 0.25) +
  geom_line(color = "#C23B22", linewidth = 0.4) +
  geom_point(color = "#C23B22", size = 1.3) +
  geom_hline(yintercept = 100 / 14, linetype = "dashed", color = "#777777", linewidth = 0.3) +
  geom_text(
    data = scree[1:2, ], aes(label = sprintf("%.1f%%", explained)),
    vjust = -0.6, size = 2.8, family = "Arial"
  ) +
  scale_x_continuous(breaks = scree$component) +
  coord_cartesian(ylim = c(0, 59), clip = "off") +
  labs(x = "Principal component", y = "Explained variance (%)") +
  theme_paper()

label_map <- c(
  pH = "pH", P = "P", K = "K", MO = "OM", Al = "Al",
  Ca = "Ca", Mg = "Mg", CTC = "CEC", SatBase = "Base sat.",
  Arg = "Clay", Silte = "Silt", Areia = "Sand",
  Desn = "Topo. var.", AberturaDossel = "Canopy openness"
)
loadings$variable <- rownames(loadings)
loadings$label <- unname(label_map[loadings$variable])
label_x <- c(
  pH = 1.03, Ca = 1.03, Mg = 1.03, SatBase = 1.03, CTC = 0.66,
  K = 0.78, MO = -0.62, P = -0.62, Areia = -0.96, Al = -0.96,
  Arg = 0.63, Silte = 0.31, Desn = -0.43, AberturaDossel = -0.42
)
label_y <- c(
  pH = -0.05, Ca = 0.07, Mg = 0.31, SatBase = 0.19, CTC = 0.56,
  K = 0.69, MO = 0.82, P = 0.68, Areia = 0.38, Al = 0.15,
  Arg = -0.59, Silte = 0.14, Desn = 0.54, AberturaDossel = -0.69
)
loadings$label_x <- unname(label_x[loadings$variable])
loadings$label_y <- unname(label_y[loadings$variable])
loadings$label_hjust <- ifelse(loadings$label_x >= 0, 0, 1)
circle <- data.frame(
  x = cos(seq(0, 2 * pi, length.out = 361)),
  y = sin(seq(0, 2 * pi, length.out = 361))
)

p2b <- ggplot() +
  geom_path(data = circle, aes(x, y), color = "#999999", linewidth = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_segment(
    data = loadings,
    aes(
      x = 0, y = 0, xend = Dim.1, yend = Dim.2
    ),
    arrow = arrow(length = grid::unit(1.5, "mm")),
    linewidth = 0.4, color = "#C15A1A", linetype = "solid"
  ) +
  geom_text(
    data = loadings,
    aes(x = label_x, y = label_y, label = label),
    size = 2.8, family = "Arial", color = "#C15A1A",
    hjust = loadings$label_hjust, vjust = 0.5
  ) +
  coord_equal(xlim = c(-1.18, 1.18), ylim = c(-1.15, 1.15), clip = "off") +
  labs(
    x = sprintf("PC1 (%.1f%%)", variance$variance.percent[1]),
    y = sprintf("PC2 (%.1f%%)", variance$variance.percent[2])
  ) +
  theme_paper() +
  theme(legend.position = "none")

scores$Area <- factor(scores$Area, levels = 1:4)
p2c <- ggplot(scores, aes(Dim.1, Dim.2, color = Area, fill = Area, shape = Area)) +
  stat_ellipse(
    geom = "polygon", type = "norm", level = 0.95,
    alpha = 0.10, linewidth = 0.35
  ) +
  geom_point(size = 2.0, stroke = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, color = "#777777") +
  scale_color_manual(values = area_colors, labels = paste("Area", 1:4)) +
  scale_fill_manual(values = area_colors, labels = paste("Area", 1:4)) +
  scale_shape_manual(values = area_shapes, labels = paste("Area", 1:4)) +
  labs(
    x = sprintf("PC1 (%.1f%%)", variance$variance.percent[1]),
    y = sprintf("PC2 (%.1f%%)", variance$variance.percent[2]),
    color = NULL, fill = NULL, shape = NULL
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = "none",
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_paper() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(4, "mm")
  )

figure2 <- (p2a / p2b / p2c) +
  plot_layout(heights = c(0.82, 1.12, 1.18)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 11))
save_paper_plot(figure2, "Fig2_environmental_PCA", 84.0, 170.7)

# Figure 3: mean net demographic changes and normal-approximation 95% CIs.
demography <- tabela_final
summary_area <- demography %>%
  group_by(Area) %>%
  summarise(
    abundance_mean = mean(Mudanca_Liq_Media),
    abundance_se = sd(Mudanca_Liq_Media) / sqrt(n()),
    basal_area_mean = mean(Mudanca_AB_Media),
    basal_area_se = sd(Mudanca_AB_Media) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    area_number = as.integer(sub("PNSJ", "", Area)),
    label = paste("Area", area_number),
    abundance_low = abundance_mean - 1.96 * abundance_se,
    abundance_high = abundance_mean + 1.96 * abundance_se,
    basal_area_low = basal_area_mean - 1.96 * basal_area_se,
    basal_area_high = basal_area_mean + 1.96 * basal_area_se
  )
figure3 <- ggplot(
  summary_area,
  aes(abundance_mean, basal_area_mean, color = factor(area_number), shape = factor(area_number))
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#555555") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#555555") +
  # Draw both uncertainty dimensions explicitly. This avoids orientation-dependent
  # rendering and gives the horizontal and vertical intervals matching end caps.
  geom_segment(
    aes(
      x = abundance_low, xend = abundance_high,
      y = basal_area_mean, yend = basal_area_mean
    ),
    linewidth = 0.6, lineend = "butt"
  ) +
  geom_segment(
    aes(
      x = abundance_low, xend = abundance_low,
      y = basal_area_mean - 0.115, yend = basal_area_mean + 0.115
    ),
    linewidth = 0.6, lineend = "butt"
  ) +
  geom_segment(
    aes(
      x = abundance_high, xend = abundance_high,
      y = basal_area_mean - 0.115, yend = basal_area_mean + 0.115
    ),
    linewidth = 0.6, lineend = "butt"
  ) +
  geom_segment(
    aes(
      x = abundance_mean, xend = abundance_mean,
      y = basal_area_low, yend = basal_area_high
    ),
    linewidth = 0.6, lineend = "butt"
  ) +
  geom_segment(
    aes(
      x = abundance_mean - 0.045, xend = abundance_mean + 0.045,
      y = basal_area_low, yend = basal_area_low
    ),
    linewidth = 0.6, lineend = "butt"
  ) +
  geom_segment(
    aes(
      x = abundance_mean - 0.045, xend = abundance_mean + 0.045,
      y = basal_area_high, yend = basal_area_high
    ),
    linewidth = 0.6, lineend = "butt"
  ) +
  geom_point(size = 3.2, stroke = 0.5, fill = "white") +
  geom_text(
    aes(label = label), nudge_y = 0.28,
    family = "Arial", fontface = "bold", size = 2.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = area_colors) +
  scale_shape_manual(values = area_shapes) +
  labs(
    x = expression("Net abundance change (" * "% " * yr^{-1} * ")"),
    y = expression("Net basal-area change (" * "% " * yr^{-1} * ")")
  ) +
  theme_paper() +
  theme(legend.position = "none")
save_paper_plot(figure3, "Fig3_demographic_change", 174.0, 131.6)

cat("Figures 2 and 3 written to results/figures\n")
