# ─────────────────────────────────────────────────────────────
# 0) Reset session
# ─────────────────────────────────────────────────────────────
rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

  library(dplyr)
  library(readr)
  library(tidyr)
  library(devEMF)
  library(corrplot)
  library(RColorBrewer)
  library(svglite)

# ─────────────────────────────────────────────────────────────
# 1) Paths
# ─────────────────────────────────────────────────────────────
DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"
DIR_STW  <- file.path(DATA_DIR, "Single_trait_walkthrough")
OUT_DIR  <- file.path(DIR_STW,   "MV5_Rg_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────
# 2) Load Rg matrices + across-year table
# ─────────────────────────────────────────────────────────────
Rg2017 <- as.matrix(read.csv(file.path(OUT_DIR, "Rg_5trait_2017.csv"),row.names   = 1,check.names = FALSE))
Rg2018 <- as.matrix(read.csv(file.path(OUT_DIR, "Rg_5trait_2018.csv"),row.names   = 1,check.names = FALSE))
rg_year_table <- read.csv(file.path(OUT_DIR, "rg_2017_2018_by_trait.csv"),check.names = FALSE)

# Alphabetical trait order
TRAITS <- sort(rownames(Rg2017))
Rg2017 <- Rg2017[TRAITS, TRAITS]
Rg2018 <- Rg2018[TRAITS, TRAITS]

trait_labels <- c(
  "in_len"      = "INL",
  "lsh_ang_1"   = "LSA",
  "tr_brz_prop" = "TBZP",
  "tr_coni"     = "TC",
  "tr_len"      = "TL")

# ─────────────────────────────────────────────────────────────
# 3) Helper to plot a single-year Rg heatmap (corrplot-style)
# ─────────────────────────────────────────────────────────────
plot_Rg_year <- function(M, traits, title_txt,
                         out_file_emf = NULL,
                         out_file_png = NULL,
                         out_file_svg = NULL) {
  
  pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  
  M_plot <- M[traits, traits]
  rownames(M_plot) <- trait_labels[rownames(M_plot)]
  colnames(M_plot) <- trait_labels[colnames(M_plot)]
  
  # SVG
  if (!is.null(out_file_svg)) {
    svglite::svglite(out_file_svg, width = 8, height = 8)
    corrplot::corrplot(
      M_plot,
      method      = "circle",
      type        = "lower",
      diag        = FALSE,
      col         = pal,
      tl.col      = "black",
      tl.srt      = 45,
      tl.cex      = 1.4,
      addCoef.col = "black",
      number.cex  = 0.8,
      addgrid.col = "grey85",
      mar         = c(0, 0, 3, 0),
      title       = title_txt
    )
    dev.off()
  }
  
  invisible(M_plot)
}

# ─────────────────────────────────────────────────────────────
# 4) Make the two Rg heatmaps (2017 & 2018)
# ─────────────────────────────────────────────────────────────
Cor2017_Rg <- plot_Rg_year(
  M             = Rg2017,
  traits        = TRAITS,
  title_txt     = "Genetic correlations (Rg) — 2017",
  out_file_svg  = file.path(OUT_DIR, "Rg_trait_corr_2017.svg")
)

Cor2018_Rg <- plot_Rg_year(
  M             = Rg2018,
  traits        = TRAITS,
  title_txt     = "Genetic correlations (Rg) — 2018",
  out_file_svg  = file.path(OUT_DIR, "Rg_trait_corr_2018.svg")
)

cat("\nSaved Rg heatmaps (SVG) in:\n", OUT_DIR, "\n")
