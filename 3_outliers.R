
# ──────────────────────────────────────────────────────────────────────────────
# 0) Reset session
rm(list=ls(all=TRUE)); graphics.off(); closeAllConnections()

# ──────────────────────────────────────────────────────────────────────────────
# 1) Libraries (minimal)
suppressPackageStartupMessages({
  library(ASRgenomics)
  library(AGHmatrix)
  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# ──────────────────────────────────────────────────────────────────────────────
# 2) Paths / switches — EDIT HERE
DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
SNP_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"
OUT_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough/outlier_removed_raw"

YEAR     <- 2018
TRAIT    <- "tr_len"       ### in_len, lsh_ang_1, tr_brz_prop, tr_coni, tr_len

FILE_17 <- file.path(DATA_DIR, "df_17_final_w_outliers.csv")
FILE_18 <- file.path(DATA_DIR, "df_18_final_w_outliers.csv")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Load this year (keep only needed columns + TRAIT) — NO Co terms anywhere
infile <- if (YEAR == 2017) FILE_17 else FILE_18
stopifnot(file.exists(infile))

meta_cols <- c("Genotype","pop","Loc","Row","Column","rep","Plot_RC_ID","Year", "co_loc")
df <- read.csv(infile, sep=";") %>%
  mutate(
    Genotype   = factor(Genotype),
    pop        = factor(pop),
    Loc   = factor(Loc),
    Row        = factor(Row),
    Column     = factor(Column),
    rep        = factor(rep),
    Plot_RC_ID = factor(Plot_RC_ID),
    Year       = YEAR
  ) %>%
  select(any_of(c(meta_cols, TRAIT)))

if (!TRAIT %in% names(df)) stop("Trait '", TRAIT, "' not found in ", basename(infile))

# ──────────────────────────────────────────────────────────────────────────────
# 4) Load GRMs (VanRaden)
Ginv <- readRDS("C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip/MASTER_Ginv_flt_imp/Ginv.rds")
# Ginv_noCo <- readRDS("C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip/MASTER_Ginv_flt_imp/Ginv_Co_window_removed.rds")




#### RERUN CODE BELOW : ####
# ──────────────────────────────────────────────────────────────────────────────
# 5) Fit single-trait ASReml model (ar1(rep))
fixed_term <- as.formula(paste(TRAIT, "~ pop"))

m_fit <- asreml(
  fixed    = fixed_term,
  random   = ~ vm(Genotype, Ginv) + at(Loc):Row + at(Loc):Column,
  residual = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | Loc),
  data     = df,
  na.action= na.method(y="include", x="include"),
  ai.sing  = FALSE
)
m_fit <- update(m_fit, aom = TRUE)

# Metrics
s   <- summary(m_fit)
LL  <- s$loglik
AIC <- s$aic
cat(sprintf("• LogLik = %.3f | AIC = %.3f\n", LL, AIC))

# ──────────────────────────────────────────────────────────────────────────────
# 5.2) Outliers (|z| > 3)
z <- m_fit$aom$R[, "stdCondRes"]
f <- as.numeric(fitted(m_fit))
r <- as.numeric(residuals(m_fit, type = "response"))

plot_df <- df |>
  mutate(Fitted = f, Resid = r, z = z,
         Outlier = ifelse(!is.na(z) & abs(z) > 3, "Outlier (|z|>3)", "Inlier"))

n_out <- sum(abs(z) > 3, na.rm=TRUE)
if (n_out == 0) {
  cat("\n✅ No outliers detected (|z|>3). Model appears stable.\n")
} else {
  cat("\n⚠️  Detected ", n_out, " outliers (|z|>3) — review plot below.\n", sep = "")
}

print(
  ggplot(plot_df, aes(Fitted, z, color = Outlier)) +
    geom_hline(yintercept = c(-3, 3), linetype = 2, linewidth = 0.5) +
    geom_point(size = 1.6, alpha = 0.9) +
    facet_wrap(~ Loc, scales = "free_x") +
    scale_color_manual(values = c("Outlier (|z|>3)" = "#d73027", "Inlier" = "#000000")) +
    labs(
      title = paste0("Standardized conditional residuals — ", TRAIT, " × ", YEAR),
      x = "Fitted values",
      y = "z = stdCondRes",
      color = ""
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = .5))
)

# Blank outliers and overwrite df in-memory for potential re-run
df[[TRAIT]][abs(z) > 3] <- NA_real_
cat(sprintf("\n• Blanked %d outliers at |z|>3 in df.\n", n_out))



##### RERUN CODE ABOVE :  ##### 
# ──────────────────────────────────────────────────────────────────────────────



# 6) Optional: Save cleaned CSV
SAVE_CSV <- TRUE
if (isTRUE(SAVE_CSV)) {
  suf     <- "pop"
  outfile <- file.path(OUT_DIR, paste0("clean_", TRAIT, "_", YEAR, "_", suf, ".csv"))
  write.csv(df, outfile, row.names = FALSE)
  cat(sprintf("\n💾 Saved: %s\n", outfile))
} else {
  cat("\n(ℹ️  Not saving CSV this run.)\n")
}

