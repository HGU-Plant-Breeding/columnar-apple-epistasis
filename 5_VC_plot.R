# ─────────────────────────────────────────────────────────────
# VC stacked barplot per Trait × Year × MutState (pop vs Co)
# Four NOIA models side-by-side in each panel — publication-ready
# ─────────────────────────────────────────────────────────────

rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(forcats)
  library(ggplot2); library(readr); library(scales)
  library(ggh4x); library(svglite)
})

DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
IN_DIR   <- file.path(DATA_DIR, "VC_Plot")

TRAITS_ORDER <- c("in_len", "lsh_ang_1", "tr_brz_prop", "tr_coni", "tr_len")

# ── Read vc_* files ───────────────────────────────────────────
files <- list.files(IN_DIR, pattern = "^vc_.*_(pop|Co)\\.csv$", full.names = TRUE)
if (!length(files)) stop("No vc_*_pop/Co.csv files found.")

vc <- lapply(files, function(f) {
  df <- read.csv(f, check.names = FALSE)
  if (!("Model" %in% names(df))) df$Model <- "GBLUP_popfixed"
  df$MutState <- sub("^vc_.*_(pop|Co)\\.csv$", "\\1", basename(f))
  df
}) |> bind_rows()

MODEL_KEEP <- c("noia_A", "noia_A_D", "noia_A_AA", "noia_A_AA_D")

keep_components <- c("Additive","A","D","AA","AD","DD","Field spatial","Residual")

vc <- vc |>
  filter(Model %in% MODEL_KEEP, Component %in% keep_components) |>
  select(Trait, Year, MutState, Model, Component, Prop) |>
  mutate(
    Trait = recode(Trait,
                   "tr_len"      = "TL",
                   "tr_coni"     = "TC",
                   "tr_brz_prop" = "TBZP",
                   "in_len"      = "INL",
                   "lsh_ang_1"   = "LSA"),
    Trait = factor(Trait, levels = c("INL","LSA","TBZP","TC","TL")),
    Year     = factor(Year, levels = sort(unique(Year))),
    MutState = factor(MutState, levels = c("pop","Co"), labels = c("pop","pop+Co")),
    Model    = factor(Model,
                      levels = MODEL_KEEP,
                      labels = c("A","A+D","A+AA","A+AA+D")),
    Component = factor(Component,
                       levels = c("Additive","A","D","AA","AD","DD","Field spatial","Residual"))
  )

fill_vals <- c(
  "Additive"      = "#1f78b4",
  "A"             = "#1f78b4",
  "D"             = "#d95f02",
  "AA"            = "#7570b3",
  "AD"            = "#e7298a",
  "DD"            = "#66a61e",
  "Field spatial" = "#cc79a7",
  "Residual"      = "#e69f00"
)

fill_labels <- c(
  "Additive"      = "Additive",
  "A"             = "Additive",
  "D"             = "Dominance",
  "AA"            = "Epistasis (AA)",
  "AD"            = "AD",
  "DD"            = "DD",
  "Field spatial" = "Field spatial",
  "Residual"      = "Residual"
)

# ─────────────────────────────────────────────────────────────
# Final publication-ready plot (no title, minimal axes)
# ─────────────────────────────────────────────────────────────

p <- ggplot(vc, aes(x = Model, y = Prop, fill = Component)) +
  geom_col(width = 0.7) +   # no black outline → cleaner
  ggh4x::facet_nested(Year ~ Trait + MutState, switch = "y") +
  scale_fill_manual(values = fill_vals, labels = fill_labels, name = "Variance component") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0,0)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y      = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(0.4, "lines"),
    panel.spacing.y    = unit(0.6, "lines"),
    strip.background   = element_rect(fill = "gray95", color = "black", linewidth = 0.4),
    strip.text.x       = element_text(size = 10, face = "bold"),
    strip.text.y.left  = element_text(angle = 0, size = 11, face = "bold"),
    legend.position    = "right",
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 11),
    plot.margin        = margin(t = 6, r = 18, b = 8, l = 8)
  )

print(p)

# ── Save as vector graphic (SVG) ─────────────────────────────
out_svg <- file.path(DATA_DIR, "vc_traits_year_model.svg")
ggsave(out_svg, plot = p, width = 16, height = 9, device = svglite)
