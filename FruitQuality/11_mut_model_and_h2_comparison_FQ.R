rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(dplyr)
  library(ggplot2)
  library(svglite)
})

FQ_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
SNP_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"
GRM_DIR <- file.path(SNP_DIR, "MASTER_Ginv_flt_imp")
OUT_DIR <- file.path(FQ_DIR, "h2_mut_fixed")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Toggles ───────────────────────────────────────────────────────────────────
TRAIT <- "CIE_a"   # JA:   anthocyanin, brix, acids, a, phenols
# FRGB: firmness, area_cm2, diameter_cm, h_w_ratio,
#        avg_R, avg_G, avg_B, circularity, CIE_a
MODEL <- "pop"    # "pop" — year + pop, full genome         (Ginv_FQ.rds)
# "co"  — year + pop + co_loc             (Ginv_FQ_noCo.rds)
#          Co window Chr10: 20–36 Mb removed
# "myb" — year + pop + myb_dose + myb_dom (Ginv_FQ_noMYB.rds)
#          MYB window Chr9: 24–40 Mb removed

# ── Derive dataset from trait ─────────────────────────────────────────────────
JA_TRAITS <- c("anthocyanin", "brix", "acids", "a", "phenols")
DATASET   <- if (TRAIT %in% JA_TRAITS) "JA" else "FRGB"
cat(sprintf("Trait: %s | Dataset: %s | Model: %s\n", TRAIT, DATASET, MODEL))

# ── Load phenotype ────────────────────────────────────────────────────────────
in_file <- file.path(FQ_DIR, if (DATASET == "JA") "juice_antho_traits_OutRmv.csv"
                     else                  "firmness_rgb_OutRmv.csv")

df <- read.csv(in_file, sep = ";") %>%
  mutate(
    Genotype = factor(Genotype),
    pop      = factor(pop),
    location = factor(location),
    year     = factor(year),
    row      = as.integer(row),
    column   = as.integer(column),
    row_f    = factor(row),
    col_f    = factor(column),
    row_ar1  = factor(row),
    col_ar1  = factor(column),
    loc_yr   = factor(paste(location, year, sep = "_")),
    myb_dose = as.numeric(myb),
    myb_dom  = as.integer(myb == 1),
    co_loc   = factor(co_loc)
  )
if (DATASET == "FRGB") df <- df %>% mutate(rep = as.integer(rep))

if (TRAIT == "anthocyanin") {
  df$pop    <- factor(df$pop, levels = c("bibon", "pxw"))
  df        <- df[df$location != "fuchsberg", ]
  df$loc_yr <- droplevels(df$loc_yr)
}

# ── Load pre-built Ginv ───────────────────────────────────────────────────────
ginv_file <- switch(MODEL,
                    "pop" = file.path(GRM_DIR, "Ginv_FQ.rds"),
                    "co"  = file.path(GRM_DIR, "Ginv_FQ_noCo.rds"),
                    "myb" = file.path(GRM_DIR, "Ginv_FQ_noMYB.rds")
)
cat(sprintf("Loading: %s\n", basename(ginv_file)))
Ginv <- readRDS(ginv_file)

# ── NA-out genotypes missing from GRM ────────────────────────────────────────
missing_from_ginv <- setdiff(levels(df$Genotype), attr(Ginv, "rowNames"))
df <- df %>%
  mutate(Genotype = if_else(Genotype %in% missing_from_ginv,
                            factor(NA, levels = levels(Genotype)),
                            Genotype))

# ── Sort & fixed formula ──────────────────────────────────────────────────────
if (DATASET == "JA") {
  df <- df[order(df$loc_yr, df$row, df$column), ]
} else {
  df <- df[order(df$location, df$year, df$row, df$column, df$rep), ]
}

fixed_extra <- switch(MODEL,
                      "pop" = NULL,
                      "co"  = "co_loc",
                      "myb" = c("myb_dose")
)
fixed_frm <- as.formula(paste(TRAIT, "~",
                              paste(c("year", "pop", fixed_extra), collapse = " + ")))

# ── Fit model ─────────────────────────────────────────────────────────────────
if (DATASET == "JA") {
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = ~ vm(Genotype, Ginv),
    residual  = ~ dsum(~ ar1(row_ar1):ar1(col_ar1) | loc_yr),
    data      = df,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE, maxit = 30
  )
} else {
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = ~ vm(Genotype, Ginv)
    + diag(loc_yr):row_f
    + diag(loc_yr):col_f,
    residual  = ~ dsum(~ idv(units) | loc_yr),
    data      = df,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE, maxit = 30
  )
}
while (!m_fit$converge) m_fit <- update(m_fit)

# ── Variance components & Cullis h² ──────────────────────────────────────────
vc         <- summary(m_fit)$varcomp
v_g        <- vc["vm(Genotype, Ginv)", "component"]
resid_pat  <- if (DATASET == "JA") "!R$" else "!units$"
resid_rows <- vc[grep(resid_pat, rownames(vc)), , drop = FALSE]
resid_rows <- resid_rows[resid_rows$bound != "B", , drop = FALSE]
v_e        <- mean(resid_rows$component)

pv    <- predict(m_fit, classify = "Genotype", sed = TRUE)
avsed <- as.numeric(pv$avsed["mean"])
h2_cu <- 1 - (avsed^2) / (2 * v_g)

cat(sprintf("\nh² Cullis (%s | %s) = %.3f  [v_g=%.4f  v_e=%.4f]\n",
            TRAIT, MODEL, h2_cu, v_g, v_e))
if (!m_fit$converge) cat("WARNING: model did not converge.\n")

# ── Save ──────────────────────────────────────────────────────────────────────
window_label <- switch(MODEL,
                       "pop" = "none",
                       "co"  = "Co_Chr10_20-36Mb",
                       "myb" = "MYB_Chr9_24-40Mb"
)

result <- data.frame(
  trait      = TRAIT,
  dataset    = DATASET,
  model      = MODEL,
  window_rmv = window_label,
  h2_Cullis  = round(h2_cu, 3),
  v_g        = round(v_g, 4),
  v_e        = round(v_e, 4),
  loglik     = round(summary(m_fit)$loglik, 3),
  aic        = round(summary(m_fit)$aic, 3)
)

out_file <- file.path(OUT_DIR, sprintf("h2_%s_%s_%s.csv", DATASET, TRAIT, MODEL))
write.csv(result, out_file, row.names = FALSE)
cat(sprintf("Saved: %s\n", basename(out_file)))

# ══════════════════════════════════════════════════════════════════════════════
# PLOT — reads all saved CSVs, shows all traits x models side by side
# Run after fitting all traits and all three MODEL settings
# ══════════════════════════════════════════════════════════════════════════════
all_files <- list.files(OUT_DIR, pattern = "^h2_.*\\.csv$", full.names = TRUE)

if (length(all_files) > 0) {
  
  all_res <- purrr::map_dfr(all_files, read.csv, check.names = FALSE) %>%
    mutate(
      model = factor(model,
                     levels = c("pop", "co", "myb"),
                     labels = c("pop only",
                                "pop + Co\n(Co window Chr10: 20-36 Mb)",
                                "pop + MYB\n(MYB window Chr9: 24-40 Mb)")),
      trait = factor(trait)
    )
  
  model_colors <- c(
    "pop only"                                  = "#4e79a7",
    "pop + Co\n(Co window Chr10: 20-36 Mb)"     = "#f28e2b",
    "pop + MYB\n(MYB window Chr9: 24-40 Mb)"    = "#59a14f"
  )
  
  p <- ggplot(all_res, aes(x = trait, y = h2_Cullis, fill = model)) +
    geom_col(position = position_dodge(0.65), width = 0.6) +
    scale_fill_manual(values = model_colors, name = "Model") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = expression(h^2~"(Cullis)")) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),
      legend.position    = "right",
      legend.title       = element_text(face = "bold"),
      plot.margin        = margin(t = 6, r = 14, b = 8, l = 8)
    )
  
  print(p)
  
  svg_path <- file.path(OUT_DIR, "h2_mut_comparison_all_traits.svg")
  ggsave(svg_path, plot = p,
         width  = max(8, length(unique(all_res$trait)) * 1.8),
         height = 6, device = svglite)
  cat(sprintf("Saved plot: %s\n", svg_path))
}