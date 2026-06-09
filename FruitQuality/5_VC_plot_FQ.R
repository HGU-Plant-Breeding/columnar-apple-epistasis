rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggh4x)
  library(svglite)
})

FQ_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
VC_DIR <- file.path(FQ_DIR, "VC_results")

# ── Read ──────────────────────────────────────────────────────────────────────
files <- list.files(VC_DIR, pattern = "^vc_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("No vc_*.csv files found in ", VC_DIR)

vc_raw <- lapply(files, read.csv, check.names = FALSE) |> dplyr::bind_rows()
names(vc_raw) <- tolower(names(vc_raw))

# ── Tidy ──────────────────────────────────────────────────────────────────────
MODEL_KEEP   <- c("A", "A_D", "A_AA", "A_AA_D")
MODEL_LABELS <- c("A", "A+D", "A+AA", "A+AA+D")

TRAIT_ORDER <- c("brix", "acids", "phenols", "area_cm2",
                 "firmness", "h_w_ratio", "circularity", "CIE_a")

vc <- vc_raw |>
  dplyr::filter(model %in% MODEL_KEEP) |>
  dplyr::mutate(
    component = factor(component,
                       levels = c("A", "D", "AA", "Spatial", "Residual")),
    model     = factor(model, levels = MODEL_KEEP, labels = MODEL_LABELS),
    trait     = factor(trait, levels = intersect(TRAIT_ORDER, unique(trait)))
  )

# ── Colours (matching original script) ───────────────────────────────────────
fill_vals <- c(
  "A"        = "#1f78b4",
  "D"        = "#d95f02",
  "AA"       = "#7570b3",
  "Spatial"  = "#cc79a7",
  "Residual" = "#e69f00"
)
fill_labels <- c(
  "A"        = "Additive",
  "D"        = "Dominance",
  "AA"       = "Epistasis (AA)",
  "Spatial"  = "Field spatial",
  "Residual" = "Residual"
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(vc, aes(x = model, y = prop, fill = component)) +
  geom_col(width = 0.7) +
  facet_wrap(~ trait, nrow = 1) +
  scale_fill_manual(values = fill_vals, labels = fill_labels,
                    name = "Variance component") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = c(0, 0), limits = c(0, 1.01)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y        = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(0.4, "lines"),
    strip.background   = element_rect(fill = "gray95", colour = "black", linewidth = 0.4),
    strip.text         = element_text(size = 10, face = "bold"),
    legend.position    = "right",
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 11),
    plot.margin        = margin(t = 6, r = 18, b = 8, l = 8)
  )

print(p)

out_svg <- file.path(FQ_DIR, "vc_FruitQuality.svg")
ggsave(out_svg, plot = p, width = 26, height = 7, device = svglite)
cat(sprintf("Saved: %s\n", out_svg))