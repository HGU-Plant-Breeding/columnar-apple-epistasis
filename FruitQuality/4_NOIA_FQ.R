rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(ASRgenomics)
  library(dplyr)
  library(data.table)
})

FQ_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
SNP_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"
VC_DIR  <- file.path(FQ_DIR, "VC_results")
dir.create(VC_DIR, showWarnings = FALSE, recursive = TRUE)

DATASET     <- "JA"    # "JA" or "FRGB"
TRAIT       <- "phenols"
KERNEL_MODE <- "A"       # "A", "A_D", "A_AA", "A_AA_D"

# JA traits:   anthocyanin, brix, acids, a, phenols
# FRGB traits: firmness, area_cm2, diameter_cm, h_w_ratio,
#              avg_R, avg_G, avg_B, circularity, CIE_a

# ── Load phenotype ────────────────────────────────────────────────────────────
in_file <- file.path(FQ_DIR, if (DATASET == "JA") "juice_antho_traits_OutRmv.csv"
                     else                  "firmness_rgb_OutRmv.csv")

df <- read.csv(in_file, sep = ";") %>%
  mutate(
    Genotype  = factor(Genotype),
    pop       = factor(pop),
    rootstock = factor(rootstock),
    location  = factor(location),
    year      = factor(year),
    row       = factor(as.integer(row)),
    column    = factor(as.integer(column)),
    row_ar1   = factor(as.integer(row)),
    col_ar1   = factor(as.integer(column)),
    loc_yr    = factor(paste(location, year, sep = "_"))
  )
if (DATASET == "FRGB") df <- df %>% mutate(rep = as.integer(rep))

cat(sprintf("Rows: %d | Non-NA %s: %d\n", nrow(df), TRAIT,
            sum(!is.na(df[[TRAIT]]))))

# ── Drop bibon for JA (too sparse for spatial terms) ─────────────────────────
if (DATASET == "JA") {
  df        <- df[df$location != "bibon", ]
  df$loc_yr <- droplevels(df$loc_yr)
  df$pop    <- droplevels(df$pop)
  cat(sprintf("Bibon dropped. Non-NA %s remaining: %d\n", TRAIT,
              sum(!is.na(df[[TRAIT]]))))
}

# ── anthocyanin restriction ───────────────────────────────────────────────────
if (TRAIT == "anthocyanin") {
  df$pop    <- factor(df$pop, levels = c("bibon", "pxw"))
  df        <- df[df$location != "fuchsberg", ]
  df$loc_yr <- droplevels(df$loc_yr)
}

# ── Build genotype matrix M ───────────────────────────────────────────────────
c130 <- fread(file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv"),
              sep = ",", data.table = FALSE)
rownames(c130) <- c130$V1

keep_ids    <- intersect(rownames(c130), as.character(df$Genotype))
c130        <- c130[keep_ids, , drop = FALSE]
M           <- as.matrix(c130[, -1, drop = FALSE])
ids_NOIA    <- rownames(M)
colnames(M) <- gsub("^AX\\.", "AX-", colnames(M))

p_alt <- colMeans(M, na.rm = TRUE) / 2
M     <- M[, pmin(p_alt, 1 - p_alt) >= 0.05, drop = FALSE]

M <- apply(M, 2, function(v) {
  if (all(is.na(v))) return(v)
  vt <- table(v[!is.na(v) & v %in% c(0, 1, 2)])
  if (!length(vt)) return(v)
  mode_val <- as.numeric(names(vt)[which.max(vt)])
  v[is.na(v)] <- mode_val
  pmin(2, pmax(0, round(as.numeric(v))))
})
M <- as.matrix(M)
storage.mode(M) <- "numeric"
rownames(M)     <- ids_NOIA

# ── Kernel construction ───────────────────────────────────────────────────────
scale_by_trace <- function(G) G / (sum(diag(G)) / nrow(G))

build_NOIA_A <- function(M) {
  n <- nrow(M)
  Z <- matrix(0, n, ncol(M))
  for (j in seq_len(ncol(M))) {
    gj <- M[, j]
    f1 <- mean(gj == 1); f2 <- mean(gj == 2)
    pA <- f2 + 0.5 * f1
    a  <- numeric(n)
    a[gj == 0] <- -2 * pA
    a[gj == 1] <-  1 - 2 * pA
    a[gj == 2] <-  2 - 2 * pA
    Z[, j] <- a
  }
  GA <- scale_by_trace(Z %*% t(Z))
  dimnames(GA) <- list(rownames(M), rownames(M))
  GA
}

GA  <- build_NOIA_A(M)
GD  <- scale_by_trace(ASRgenomics::G.matrix(M = M, method = "Vitezica")$G)
GAA <- scale_by_trace(GA * GA)

diag(GA)  <- diag(GA)  + 1e-3
diag(GD)  <- diag(GD)  + 1e-3
diag(GAA) <- diag(GAA) + 1e-3

Ginv_A  <- G.inverse(GA,  sparseform = TRUE)$Ginv
Ginv_D  <- G.inverse(GD,  sparseform = TRUE)$Ginv
Ginv_AA <- G.inverse(GAA, sparseform = TRUE)$Ginv

df$GenotypeK <- factor(
  ifelse(df$Genotype %in% rownames(GA), as.character(df$Genotype), NA),
  levels = rownames(GA)
)

# ── Fixed formula ─────────────────────────────────────────────────────────────
fixed_frm <- as.formula(paste(TRAIT, "~ pop + rootstock"))

# ── Kernel random term ────────────────────────────────────────────────────────
kernel_rhs <- switch(
  KERNEL_MODE,
  "A"      = "vm(GenotypeK, Ginv_A)",
  "A_D"    = "vm(GenotypeK, Ginv_A) + vm(GenotypeK, Ginv_D)",
  "A_AA"   = "vm(GenotypeK, Ginv_A) + vm(GenotypeK, Ginv_AA)",
  "A_AA_D" = "vm(GenotypeK, Ginv_A) + vm(GenotypeK, Ginv_AA) + vm(GenotypeK, Ginv_D)"
)

# ── Fit model ─────────────────────────────────────────────────────────────────
if (DATASET == "JA") {
  # JA: kernel only in random; spatial gradient captured by AR1 residual
  df_fit     <- df[order(df$loc_yr, df$row, df$column), ]
  random_frm <- as.formula(paste("~", kernel_rhs,
                                 "+ at(loc_yr):row + at(loc_yr):column"))
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = random_frm,
    residual  = ~ dsum(~ ar1(row):ar1(column) | loc_yr),
    data      = df_fit,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE,
    maxit     = 30,
    workspace = "1gb"
  )
} else {
  # FRGB: kernel + location-specific row/column spatial in random; IID fruit residual
  df_fit     <- df[order(df$location, df$year, df$row, df$column, df$rep), ]
  random_frm <- as.formula(paste("~", kernel_rhs,
                                 "+ at(loc_yr):row + at(loc_yr):column"))
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = random_frm,
    residual  = ~ dsum(~ idv(units) | loc_yr),
    data      = df_fit,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE,
    maxit     = 30,
    workspace = "1gb"
  )
}

while (!m_fit$converge) m_fit <- update(m_fit)

# ── Diagnostics ───────────────────────────────────────────────────────────────
s  <- summary(m_fit)
vc <- s$varcomp
cat(sprintf("\nTrait: %s | Dataset: %s | Kernel: %s | LogLik = %.3f | AIC = %.3f\n",
            TRAIT, DATASET, KERNEL_MODE, s$loglik, s$aic))
print(vc)
if (!m_fit$converge) cat("WARNING: model did not converge.\n")
wald.asreml(m_fit)

# ── Variance partitioning ─────────────────────────────────────────────────────
rn <- rownames(vc)
v  <- vc$component; names(v) <- rn

var_A  <- sum(v[grepl("Ginv_A\\b", rn) & !grepl("Ginv_AA", rn)], na.rm = TRUE)
var_D  <- sum(v[grepl("Ginv_D\\b", rn)], na.rm = TRUE)
var_AA <- sum(v[grepl("Ginv_AA\\b", rn)], na.rm = TRUE)

# spatial: at(loc_yr):row and at(loc_yr):column terms (FRGB only)
var_sp <- if (any(grepl(":row$|:column$", rn))) {
  sum(v[grepl(":row$|:column$", rn) & vc$bound != "B"], na.rm = TRUE)
} else { 0 }

resid_pat  <- if (DATASET == "JA") "!R$" else "!units$"
resid_rows <- vc[grep(resid_pat, rn), , drop = FALSE]
resid_rows <- resid_rows[resid_rows$bound != "B", , drop = FALSE]
var_e      <- mean(resid_rows$component)

comp_vals <- c(
  A        = var_A,
  D        = if (var_D  > 0) var_D  else NULL,
  AA       = if (var_AA > 0) var_AA else NULL,
  Spatial  = if (var_sp > 0) var_sp else NULL,
  Residual = var_e
)
var_total <- sum(comp_vals)
props     <- comp_vals / var_total

cat(sprintf("\nVariance partitioning (%s):\n", KERNEL_MODE))
for (nm in names(comp_vals))
  cat(sprintf("  %-10s = %.4f  (prop = %.3f)\n", nm, comp_vals[nm], props[nm]))

# ── Save VC table ─────────────────────────────────────────────────────────────
out_df <- data.frame(
  trait     = TRAIT,
  dataset   = DATASET,
  component = names(comp_vals),
  var       = as.numeric(comp_vals),
  prop      = as.numeric(props),
  model     = KERNEL_MODE,
  loglik    = rep(s$loglik, length(comp_vals)),
  aic       = rep(s$aic,    length(comp_vals)),
  stringsAsFactors = FALSE
)

out_file <- file.path(VC_DIR, sprintf("vc_%s_%s.csv", TRAIT, KERNEL_MODE))
write.csv(out_df, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))