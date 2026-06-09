rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
  library(RColorBrewer)
  library(svglite)
})

FQ_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
IN_DIR  <- file.path(FQ_DIR, "GEBV_results")
OUT_DIR <- file.path(FQ_DIR, "GEBV_corr_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TRAIT_ORDER  <- c("brix", "acids", "phenols", "firmness", "area_cm2", "CIE_a", "h_w_ratio")
TRAIT_LABELS <- c("Brix", "Acids", "Phenols", "Firmness", "Area", "CIE_a", "H/W")

# ── Read all GEBV files — trait extracted from filename ───────────────────────
# filenames: gebv_JA_brix.csv, gebv_FRGB_firmness.csv etc.
files <- list.files(IN_DIR, pattern = "^gebv_.*\\.csv$", full.names = TRUE)

df <- purrr::map_dfr(files, function(f) {
  d <- read.csv(f, check.names = FALSE)
  d$trait <- sub("^gebv_(.+)\\.csv$", "\\1", basename(f))
  d
}) %>%
  transmute(
    Genotype = as.character(Genotype),
    trait    = as.character(trait),
    dGEBV    = as.numeric(dGEBV),
    pop_BLUE = as.numeric(pop_BLUE),
    value    = dGEBV + pop_BLUE
  ) %>%
  mutate(trait = factor(trait, levels = TRAIT_ORDER, labels = TRAIT_LABELS)) %>%
  filter(!is.na(trait), is.finite(value))

traits <- levels(df$trait)
cat(sprintf("Traits loaded: %s\n", paste(traits, collapse = ", ")))

# ── Trait-trait correlation heatmap ───────────────────────────────────────────
wide <- df %>%
  select(Genotype, trait, value) %>%
  distinct() %>%
  pivot_wider(names_from = trait, values_from = value) %>%
  select(Genotype, any_of(traits))

X <- wide %>%
  column_to_rownames("Genotype") %>%
  select(where(~ sum(!is.na(.)) > 10))

M   <- cor(X, use = "pairwise.complete.obs")
pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

svglite(file.path(OUT_DIR, "dGEBV_trait_corr_FQ.svg"), width = 8, height = 8)
corrplot(
  M,
  method      = "circle",
  type        = "lower",
  diag        = FALSE,
  col         = pal,
  tl.col      = "black",
  tl.srt      = 45,
  tl.cex      = 1.2,
  addCoef.col = "black",
  number.cex  = 0.85,
  addgrid.col = "grey85",
  mar         = c(0, 0, 1, 0)
)
dev.off()
cat("Saved: dGEBV_trait_corr_FQ.svg\n")

write.csv(as.data.frame(M),
          file.path(OUT_DIR, "dGEBV_trait_corr_matrix.csv"))
print(round(M, 3))