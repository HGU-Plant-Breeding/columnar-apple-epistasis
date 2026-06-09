# ─────────────────────────────────────────────────────────────
# 0) Reset session
# ─────────────────────────────────────────────────────────────
rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
  library(RColorBrewer)
  library(devEMF)
  library(svglite)
})

# ─────────────────────────────────────────────────────────────
# 1) Paths
# ─────────────────────────────────────────────────────────────
BASE_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
IN_DIR   <- file.path(BASE_DIR, "deregression_output")
OUT_DIR  <- file.path(BASE_DIR, "dGEBV_corr_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────
# 2) Read all dGEBV files (NO pop_BLUE used here)
# ─────────────────────────────────────────────────────────────
files <- list.files(IN_DIR, pattern = "^dGEBV_.*_pop\\.csv$", full.names = TRUE)

df <- purrr::map_dfr(files, ~ read.csv(.x, check.names = FALSE)) %>%
  transmute(
    Year     = as.integer(Year),
    Trait    = as.character(Trait),
    Genotype = as.character(Genotype),
    dGEBV    = as.numeric(dGEBV)
  ) %>%
  mutate(
    Trait = recode(Trait,
                   "in_len"       = "INL",
                   "lsh_ang_1"    = "LSA",
                   "tr_brz_prop"  = "TBZP",
                   "tr_coni"      = "TC",
                   "tr_len"       = "TL"),
    Trait = factor(Trait,
                   levels = c("INL","LSA","TBZP","TC","TL"))
  )

traits <- levels(df$Trait)
years  <- sort(unique(df$Year))

# ─────────────────────────────────────────────────────────────
# 3) Within-year triangular heatmaps (trait–trait cor of dGEBV)
# ─────────────────────────────────────────────────────────────
make_corrplot_year <- function(df, year, traits,
                               out_file_emf = NULL,
                               out_file_svg = NULL,
                               out_file_png = NULL,
                               title_txt) {
  wide <- df %>%
    filter(Year == year, is.finite(dGEBV)) %>%
    select(Genotype, Trait, dGEBV) %>%
    distinct() %>%
    pivot_wider(names_from = Trait, values_from = dGEBV)
  
  keep <- intersect(traits, colnames(wide))
  X    <- as.data.frame(wide[, keep, drop = FALSE])
  
  M   <- cor(X, use = "pairwise.complete.obs")
  pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

  # SVG
  if (!is.null(out_file_svg)) {
    svglite::svglite(out_file_svg, width = 8, height = 8)
    corrplot(
      M,
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
      title       = title_txt,
      mar         = c(0, 0, 3, 0)
    )
    dev.off()
  }
  
  invisible(M)
}

Cor2017 <- make_corrplot_year(
  df            = df,
  year          = 2017,
  traits        = traits,
  out_file_svg  = file.path(OUT_DIR, "dGEBV_trait_corr_2017.svg"),
  title_txt     = ""
)

Cor2018 <- make_corrplot_year(
  df            = df,
  year          = 2018,
  traits        = traits,
  out_file_svg  = file.path(OUT_DIR, "dGEBV_trait_corr_2018.svg"),
  title_txt     = ""
)

# ─────────────────────────────────────────────────────────────
# 4) Cross-year correlations per trait (with SE & CI)
# ─────────────────────────────────────────────────────────────
cross_corr_trait <- function(tr) {
  d17 <- df %>% filter(Trait == tr, Year == 2017, is.finite(dGEBV)) %>%
    select(Genotype, dGEBV)
  d18 <- df %>% filter(Trait == tr, Year == 2018, is.finite(dGEBV)) %>%
    select(Genotype, dGEBV)
  
  m <- inner_join(d17, d18, by = "Genotype", suffix = c("_2017", "_2018"))
  N_eff <- nrow(m)
  
  if (N_eff < 4) {
    return(
      tibble(
        Trait    = tr,
        N_eff    = N_eff,
        r        = NA_real_,
        z        = NA_real_,
        SE_z     = NA_real_,
        z_low    = NA_real_,
        z_high   = NA_real_,
        r_low    = NA_real_,
        r_high   = NA_real_,
        SE_r     = NA_real_,
        label    = NA_character_
      )
    )
  }
  
  r <- cor(m$dGEBV_2017, m$dGEBV_2018, use = "pairwise.complete.obs")
  
  # Fisher z transform
  z    <- 0.5 * log((1 + r) / (1 - r))
  SE_z <- 1 / sqrt(N_eff - 3)
  
  z_low  <- z - 1.96 * SE_z
  z_high <- z + 1.96 * SE_z
  
  r_low  <- tanh(z_low)
  r_high <- tanh(z_high)
  
  # approximate SE of r from CI
  SE_r <- (r_high - r_low) / (2 * 1.96)
  
  tibble(
    Trait  = tr,
    N_eff  = N_eff,
    r      = r,
    z      = z,
    SE_z   = SE_z,
    z_low  = z_low,
    z_high = z_high,
    r_low  = r_low,
    r_high = r_high,
    SE_r   = SE_r,
    label  = sprintf("%.3f ± %.3f", r, SE_r)
  )
}

cross_year_tbl <- purrr::map_dfr(traits, cross_corr_trait) %>% arrange(Trait)

# ─────────────────────────────────────────────────────────────
# 5) Save cross-year correlation table
# ─────────────────────────────────────────────────────────────
out_csv <- file.path(OUT_DIR, "cross_year_dGEBV_correlations.csv")
write.csv(cross_year_tbl, file = out_csv, row.names = FALSE)

print(cross_year_tbl)
