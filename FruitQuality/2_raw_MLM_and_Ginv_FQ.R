rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

  library(ASRgenomics)
  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(dplyr)
  library(data.table)

FQ_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
SNP_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"

DATASET <- "FRGB"   # "JA" or "FRGB"
TRAIT   <- "h_w_ratio"

# JA traits:   anthocyanin, brix, acids, a, phenols
# FRGB traits: firmness, size_cm2, diameter_cm, h_w_ratio, avg_R, avg_G, avg_B, circularity, CIE_a

df_raw <- read.csv(
  file.path(FQ_DIR, if (DATASET == "JA") "juice_antho_traits_FULL.csv" else "firmness_rgb_FULL.csv"),
  sep = ";"
)

df <- df_raw %>%
  mutate(
    Genotype = factor(Genotype),
    pop      = factor(pop),
    rootstock= factor(rootstock),
    location = factor(location),
    year     = factor(year),
    row      = factor(row),
    column   = factor(column),
    loc_yr   = factor(paste(location, year, sep = "_"))
  )
if (DATASET == "FRGB") df <- df %>% mutate(rep = as.integer(rep))

# ── GRM ───────────────────────────────────────────────────────────────────────
# Ginv is cached. Downstream scripts only need to readRDS() from grm_dir.
grm_dir    <- file.path(SNP_DIR, "MASTER_Ginv_flt_imp")
ginv_cache <- file.path(grm_dir, "Ginv_FQ.rds")

if (file.exists(ginv_cache)) {
  Ginv <- readRDS(ginv_cache)
} else {
  c130 <- fread(file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv"),
                sep = ",", data.table = FALSE)
  rownames(c130) <- c130$V1
  c130 <- c130[, -1, drop = FALSE]
  colnames(c130) <- gsub("^AX\\.", "AX-", colnames(c130))
  M     <- as.matrix(c130)
  p_alt <- colMeans(M, na.rm = TRUE) / 2
  M     <- M[, pmin(p_alt, 1 - p_alt) >= 0.05, drop = FALSE]
  for (j in which(colSums(is.na(M)) > 0))
    M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
  M <- round(M)
  G       <- ASRgenomics::G.matrix(M = M, method = "VanRaden", sparseform = FALSE)$G
  diag(G) <- diag(G) + 1e-2
  Ginv    <- ASRgenomics::G.inverse(G = G, sparseform = TRUE)$Ginv
  dir.create(grm_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(Ginv, ginv_cache)
}

# NA-out phenotyped trees absent from the SNP file so asreml doesn't warn
missing_from_ginv <- setdiff(levels(df$Genotype), attr(Ginv, "rowNames"))
df <- df %>%
  mutate(Genotype = if_else(Genotype %in% missing_from_ginv,
                            factor(NA, levels = levels(Genotype)),
                            Genotype))

# ── Model ─────────────────────────────────────────────────────────────────────
if (TRAIT == "anthocyanin") {
  df$pop <- factor(df$pop, levels = c("bibon", "pxw"))
  empty_locyr <- levels(df$loc_yr)[tapply(df[[TRAIT]], df$loc_yr, function(x) all(is.na(x)))]
  df <- df[!df$location %in% "fuchsberg", ]
  df$loc_yr <- droplevels(df$loc_yr)
}


fixed_frm <- as.formula(paste(TRAIT, "~ pop + rootstock"))

if (DATASET == "JA") {
  
  df_fit <- df %>% arrange(loc_yr, row, column)
  
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = ~ corh(loc_yr):vm(Genotype, Ginv) + at(loc_yr):row + at(loc_yr):column,
    residual  = ~ dsum(~ ar1(row):ar1(column)|loc_yr),
    data      = df_fit,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE,
    maxit     = 30,
    workspace = "2gb"
  )
  
} else {
  
  df_fit <- df %>% arrange(location, year, row, column, rep)
  
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = ~ corh(loc_yr):vm(Genotype, Ginv) + at(loc_yr):row + at(loc_yr):column,
    residual  = ~ dsum(~ idv(units)| loc_yr),
    data      = df_fit,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE,
    maxit     = 30,
    workspace = "2gb"
  )
}

m_fit <- update(m_fit, aom = TRUE)
m_fit <- update(m_fit, aom = TRUE)

# ── Diagnostics ───────────────────────────────────────────────────────────────
s  <- summary(m_fit)
s$loglik; s$aic; s$bic


vc <- s$varcomp

print(vc)

if (!m_fit$converge) cat("WARNING: model did not converge.\n")

# ── h² ────────────────────────────────────────────────────────────────────────
rn <- rownames(vc)

v_g_rows   <- vc[grep("loc_yr:vm\\(Genotype, Ginv\\)", rn), , drop = FALSE]
v_g_rows   <- v_g_rows[v_g_rows$bound != "B", , drop = FALSE]

resid_pat  <- if (DATASET == "JA") "!R$" else "!units$"
resid_rows <- vc[grep(resid_pat, rn), , drop = FALSE]
resid_rows <- resid_rows[resid_rows$bound != "B", , drop = FALSE]

# per-environment h² and pooled mean
v_g_vec <- v_g_rows$component
v_e_vec <- resid_rows$component

h2_per_env <- v_g_vec / (v_g_vec + v_e_vec)
names(h2_per_env) <- sub(".*!loc_yr_", "", rownames(v_g_rows))

cat("\nh² per environment:\n")
print(round(h2_per_env, 3))
cat(sprintf("\nh² (mean across environments) = %.3f\n", mean(h2_per_env)))


wald.asreml(m_fit)
plot(m_fit)