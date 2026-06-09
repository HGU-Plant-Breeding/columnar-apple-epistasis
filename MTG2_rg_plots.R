rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

library(tidyverse)
library(devEMF)

#####################################################################
# 0) PATHS & TRAITS
#####################################################################

MTG2_OUT_DIR <- "C:/Users/nguevenc/mtg2_project/mtg2_out"
PLOT_DIR     <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"

traits <- c("in_len", "lsh_ang_1", "tr_brz_prop", "tr_coni", "tr_len")

#####################################################################
# 1) CROSS-YEAR rg BARPLOT (FROM rg_* FILES)
#####################################################################

# Read MTG2 bivariate output files
rg_in_len_lines      <- readLines(file.path(MTG2_OUT_DIR, "rg_in_len"))
rg_lsh_ang_1_lines   <- readLines(file.path(MTG2_OUT_DIR, "rg_lsh_ang_1"))
rg_tr_brz_prop_lines <- readLines(file.path(MTG2_OUT_DIR, "rg_tr_brz_prop"))
rg_tr_coni_lines     <- readLines(file.path(MTG2_OUT_DIR, "rg_tr_coni"))
rg_tr_len_lines      <- readLines(file.path(MTG2_OUT_DIR, "rg_tr_len"))

# Extract rg + SE from one MTG2 rg file
extract_rg_from_lines <- function(lines) {
  idx <- grep("^\\s*cor", lines)
  line  <- lines[idx[1]]
  parts <- strsplit(trimws(line), "\\s+")[[1]]  # "cor  rg  se"
  list(rg = as.numeric(parts[2]),
       se = as.numeric(parts[3]))
}

# Apply to all traits
res_in_len      <- extract_rg_from_lines(rg_in_len_lines)
res_lsh_ang_1   <- extract_rg_from_lines(rg_lsh_ang_1_lines)
res_tr_brz_prop <- extract_rg_from_lines(rg_tr_brz_prop_lines)
res_tr_coni     <- extract_rg_from_lines(rg_tr_coni_lines)
res_tr_len      <- extract_rg_from_lines(rg_tr_len_lines)

cross_rg <- tibble(
  trait = factor(traits, levels = traits),
  rg    = c(
    res_in_len$rg,
    res_lsh_ang_1$rg,
    res_tr_brz_prop$rg,
    res_tr_coni$rg,
    res_tr_len$rg
  ),
  se    = c(
    res_in_len$se,
    res_lsh_ang_1$se,
    res_tr_brz_prop$se,
    res_tr_coni$se,
    res_tr_len$se
  )
) %>%
  mutate(
    ci_lower = pmax(0, rg - 1.96 * se),
    ci_upper = pmin(1, rg + 1.96 * se)
  )

print(cross_rg)

gg_cross_rg <- ggplot(cross_rg, aes(x = trait, y = rg)) +
  geom_col(width = 0.7, fill = "steelblue") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "Cross-year genetic correlation (rg)",
    title = "Genetic correlation between years"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x        = element_text(angle = 30, hjust = 1),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

#####################################################################
# 2) EXTRACT 5×5 GENETIC rg MATRICES FROM mv5_2017 / mv5_2018
#####################################################################

mv5_2017_lines <- readLines(file.path(MTG2_OUT_DIR, "mv5_2017"))
mv5_2018_lines <- readLines(file.path(MTG2_OUT_DIR, "mv5_2018"))

# Helper: build G matrix (additive covariance) from MTG2 mv5_* output
extract_G_from_mv5 <- function(lines, n_trait = 5) {
  # take the last n_trait Va lines and last choose(n_trait,2) cova lines
  idx_va   <- grep("^\\s*Va",   lines)
  idx_cova <- grep("^\\s*cova", lines)
  
  if (length(idx_va)   < n_trait)            stop("Too few Va lines.")
  if (length(idx_cova) < choose(n_trait, 2)) stop("Too few cova lines.")
  
  va_lines   <- lines[tail(idx_va,   n_trait)]
  cova_lines <- lines[tail(idx_cova, choose(n_trait, 2))]
  
  parse_est <- function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    as.numeric(parts[2])  # estimate is the 2nd token
  }
  
  Va   <- vapply(va_lines,   parse_est, numeric(1))
  cova <- vapply(cova_lines, parse_est, numeric(1))
  
  # Build covariance matrix G
  G <- matrix(0, n_trait, n_trait)
  diag(G) <- Va
  
  pairs <- combn(n_trait, 2)  # (1,2), (1,3), ..., (4,5)
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]
    j <- pairs[2, k]
    G[i, j] <- cova[k]
    G[j, i] <- cova[k]
  }
  
  rownames(G) <- traits[seq_len(n_trait)]
  colnames(G) <- traits[seq_len(n_trait)]
  
  G
}

G2017 <- extract_G_from_mv5(mv5_2017_lines, n_trait = 5)
G2018 <- extract_G_from_mv5(mv5_2018_lines, n_trait = 5)

Rg2017 <- cov2cor(G2017)
Rg2018 <- cov2cor(G2018)
#####################################################################
# 3) BUBBLE-STYLE rg PLOT FUNCTION (like your example figure)
#####################################################################

plot_rg_bubbles <- function(Rg_matrix, title) {
  stopifnot(is.matrix(Rg_matrix))
  
  df_long <- as.data.frame(Rg_matrix) %>%
    mutate(trait_row = rownames(Rg_matrix)) %>%
    pivot_longer(
      cols      = all_of(traits),
      names_to  = "trait_col",
      values_to = "rg"
    ) %>%
    mutate(
      trait_row = factor(trait_row, levels = traits),
      trait_col = factor(trait_col, levels = traits)
    ) %>%
    # keep only lower triangle (row > col)
    dplyr::filter(as.integer(trait_row) > as.integer(trait_col))
  
  ggplot(df_long, aes(x = trait_col, y = trait_row)) +
    # circles
    geom_point(aes(size = abs(rg), fill = rg),
               shape = 21, colour = NA, alpha = 0.9) +
    # rg labels
    geom_text(aes(label = sprintf("%.2f", rg)), size = 3) +
    scale_size(range = c(3, 18), guide = "none") +
    scale_fill_gradient2(
      limits   = c(-1, 1),
      midpoint = 0,
      name     = "rg"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = title) +
    coord_fixed() +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

gg_rg_2017 <- plot_rg_bubbles(Rg2017, "Genetic correlations among traits (2017)")
gg_rg_2018 <- plot_rg_bubbles(Rg2018, "Genetic correlations among traits (2018)")


#####################################################################
# 4) SHOW PLOTS
#####################################################################

print(gg_cross_rg)
print(gg_rg_2017)
print(gg_rg_2018)

#####################################################################
# 5) SAVE AS EMF
#####################################################################

# Barplot
emf(file = file.path(PLOT_DIR, "rg_crossyear_barplot.emf"),
    width = 6, height = 5, emfPlus = TRUE)
print(gg_cross_rg)
dev.off()

# Heatmap 2017
emf(file = file.path(PLOT_DIR, "rg_heatmap_2017.emf"),
    width = 6, height = 5, emfPlus = TRUE)
print(gg_rg_2017)
dev.off()

# Heatmap 2018
emf(file = file.path(PLOT_DIR, "rg_heatmap_2018.emf"),
    width = 6, height = 5, emfPlus = TRUE)
print(gg_rg_2018)
dev.off()
