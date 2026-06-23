# ─────────────────────────────────────────────────────────────
# 0) Reset session
# ─────────────────────────────────────────────────────────────
rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

# ─────────────────────────────────────────────────────────────
# 1) Libraries
# ─────────────────────────────────────────────────────────────
library(asreml)
vm <- get("asr_vm", envir = asNamespace("asreml"))
library(ASRgenomics)
library(AGHmatrix)
library(dplyr)
library(data.table)

# ─────────────────────────────────────────────────────────────
# 2) Paths / switches — EDIT HERE (trait/year/Co_PCR only)
# ─────────────────────────────────────────────────────────────
DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"
SNP_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"

TRAIT   <- "in_len"      # "in_len", "lsh_ang_1", "tr_brz_prop", "tr_coni", "tr_len"
YEAR    <- 2017          # 2017 or 2018
Co_PCR  <- FALSE       # include co_loc as fixed AND drop LD window around Co locus

KERNEL_MODE <- "A"   # "A", "A_D", "A_AA", "A_AA_D"

# ─────────────────────────────────────────────────────────────
# 3) Load the single cleaned phenotype file
# ─────────────────────────────────────────────────────────────
dir_stw   <- file.path(DATA_DIR, "Single_trait_walkthrough", "outlier_removed_raw")
phenopath <- list.files(
  dir_stw,
  pattern    = sprintf("^clean_%s_%d_.*\\.csv$", TRAIT, YEAR),
  full.names = TRUE
)
stopifnot(length(phenopath) == 1)
phenopath <- phenopath[1]

meta_cols <- c("Genotype","pop","Loc","Row","Column","rep","Plot_RC_ID","co_loc")

df <- read.csv(phenopath, sep = ",") |>
  mutate(
    Genotype   = factor(Genotype),
    pop        = factor(pop),
    Loc   = factor(Loc),
    Row        = factor(Row),
    Column     = factor(Column),
    rep        = factor(rep),
    Plot_RC_ID = factor(Plot_RC_ID),
    co_loc = factor(co_loc)
  ) |>
  dplyr::select(any_of(c(meta_cols, TRAIT)))

cat("Loaded phenotype file: ", basename(phenopath), "\n", sep = "")

# ─────────────────────────────────────────────────────────────
# 4) Genotype matrix M: filter, Co window, MAF, impute 0/1/2
# ─────────────────────────────────────────────────────────────
c130 <- fread(
  file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv"),
  sep = ",", data.table = FALSE
)
rownames(c130) <- c130$V1

keep_ids <- intersect(rownames(c130), as.character(df$Genotype))
stopifnot(length(keep_ids) > 1)
c130 <- c130[keep_ids, , drop = FALSE]

M <- as.matrix(c130[, -1, drop = FALSE])
ids_NOIA <- rownames(M)
colnames(M) <- gsub("^AX\\.", "AX-", colnames(M))

# Optional LD window removal around Columnar locus (Chr10: 20–36 Mb)
if (Co_PCR) {
  GM  <- read.csv(file.path(SNP_DIR, "GM_v11.csv"), sep = ";")
  win <- GM$chrom == 10 & GM$pos >= 2e7 & GM$pos <= 3.6e7
  to_drop <- unique(GM$name[win])
  M <- M[, !colnames(M) %in% to_drop, drop = FALSE]
}

# MAF filter using allele freq from 0/1/2 coding
p_alt <- colMeans(M, na.rm = TRUE) / 2
maf   <- pmin(p_alt, 1 - p_alt)
M     <- M[, maf >= 0.05, drop = FALSE]

# Mode-impute to 0/1/2 per marker (NOIA / Vitezica need discrete genotypes)
M <- apply(M, 2, function(v) {
  if (all(is.na(v))) return(v)
  vt <- table(v[!is.na(v) & v %in% c(0, 1, 2)])
  if (!length(vt)) return(v)
  m <- as.numeric(names(vt)[which.max(vt)])
  v[is.na(v)] <- m
  pmin(2, pmax(0, round(as.numeric(v))))
})
M <- as.matrix(M)
storage.mode(M) <- "numeric"
rownames(M) <- ids_NOIA
stopifnot(!is.null(rownames(M)))

# ─────────────────────────────────────────────────────────────
# 5) Additive (NOIA), dominance, and AA epistasis kernels
# ─────────────────────────────────────────────────────────────
# scale, so trace(G)/n = 1 → average diagonal = 1
scale_by_trace <- function(G) {
  n <- nrow(G)
  s <- sum(diag(G)) / n
  G / s
}

# Vitezica-style additive coding (NOIA)
build_NOIA_A <- function(M) {
  n <- nrow(M)
  p <- ncol(M)
  Z <- matrix(0, n, p)
  
  for (j in seq_len(p)) {
    gj <- M[, j]          # genotypes at marker j as 0/1/2
    
    # observed genotype frequencies in the sample
    f0 <- mean(gj == 0)   # freq(aa)
    f1 <- mean(gj == 1)   # freq(Aa)
    f2 <- mean(gj == 2)   # freq(AA)
    
    # allele frequencies from genotype frequencies
    pA <- f2 + 0.5 * f1   # freq(A)
    qA <- 1 - pA          # freq(a)
    
    # NOIA additive scores (Vitezica idea)
    a <- numeric(n)
    a[gj == 0] <- -2 * pA        # aa
    a[gj == 1] <-  (qA - pA)     # Aa
    a[gj == 2] <-  2 * qA        # AA
    
    Z[, j] <- a
  }
  
  # additive relationship: GA_ij = sum_k a_ik * a_jk
  GA <- Z %*% t(Z)
  dimnames(GA) <- list(rownames(M), rownames(M))
  
  GA <- scale_by_trace(GA)   # normalise so mean diag ≈ 1
  GA
}

# Additive kernel (NOIA / Vitezica coding)
GA <- build_NOIA_A(M)

# Dominance kernel from ASRgenomics (Vitezica dominance matrix)
GD <- ASRgenomics::G.matrix(M = M, method = "Vitezica")$G
GD <- scale_by_trace(GD)


GAA <- scale_by_trace(GA * GA)                                # AA epistasis

# Small ridge on diagonals to stabilise inversion
diag(GA)  <- diag(GA)  + 1e-3
diag(GD)  <- diag(GD)  + 1e-3
diag(GAA) <- diag(GAA) + 1e-3

Ginv_A  <- G.inverse(GA,  sparseform = FALSE)$Ginv
Ginv_D  <- G.inverse(GD,  sparseform = FALSE)$Ginv
Ginv_AA <- G.inverse(GAA, sparseform = FALSE)$Ginv

ids_in_K <- rownames(Ginv_A)

df$GenotypeK <- ifelse(df$Genotype %in% ids_in_K, as.character(df$Genotype), NA)
df$GenotypeK <- factor(df$GenotypeK, levels = ids_in_K)

# ─────────────────────────────────────────────────────────────
# 7) Fixed structure
# ─────────────────────────────────────────────────────────────
fixed_term <- as.formula(paste(TRAIT, ifelse(Co_PCR, "~ co_loc + pop", "~ pop")))

cat("\nFitting single-year NOIA for: ", TRAIT, 
    "  | YEAR: ", YEAR, "\n", sep = "")

# ─────────────────────────────────────────────────────────────
# 8) Kernel
# ─────────────────────────────────────────────────────────────

cat("KERNEL_MODE = ", KERNEL_MODE, "\n", sep = "")

kernel_term <- switch(
  KERNEL_MODE,
  "A"      = "vm(GenotypeK, Ginv_A)",
  "A_D"    = paste("vm(GenotypeK, Ginv_A)",
                   "vm(GenotypeK, Ginv_D)",  sep = " + "),
  "A_AA"   = paste("vm(GenotypeK, Ginv_A)",
                   "vm(GenotypeK, Ginv_AA)", sep = " + "),
  "A_AA_D" = paste("vm(GenotypeK, Ginv_A)",
                   "vm(GenotypeK, Ginv_AA)",
                   "vm(GenotypeK, Ginv_D)",  sep = " + "))

random_formula <- as.formula(paste("~", kernel_term, "+ at(Loc):Row + at(Loc):Column"))

# ─────────────────────────────────────────────────────────────
# 9) Fit single-year NOIA model
# ─────────────────────────────────────────────────────────────
m_fit <- asreml(
  fixed     = fixed_term,
  random    = random_formula,
  residual  = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | Loc),
  data      = df,
  na.action = na.method(y = "include", x = "include"),
  ai.sing   = FALSE
)

m_fit <- update(m_fit, aom = TRUE)

# ─────────────────────────────────────────────────────────────
# 10) Model metrics, var comps → save per kernel for VC plots
# ─────────────────────────────────────────────────────────────
s  <- summary(m_fit)
vc <- s$varcomp
LL <- s$loglik
AIC <- s$aic

cat(sprintf("\nLogLik = %.3f | AIC = %.3f\n\n", LL, AIC))
print(vc)

vc_df <- as.data.frame(vc)
rn    <- rownames(vc_df)
v     <- vc_df$component; names(v) <- rn

idx_A  <- grepl("GenotypeK", rn) & grepl("Ginv_A\\b",  rn)
idx_D  <- grepl("GenotypeK", rn) & grepl("Ginv_D\\b",  rn)
idx_AA <- grepl("GenotypeK", rn) & grepl("Ginv_AA\\b", rn)

var_A  <- sum(v[idx_A],  na.rm = TRUE)
var_D  <- sum(v[idx_D],  na.rm = TRUE)
var_AA <- sum(v[idx_AA], na.rm = TRUE)

var_spatial <- mean(v[grepl("Row|Column", rn)], na.rm = TRUE)
var_resid   <- mean(v[grepl("!R$", rn)],        na.rm = TRUE)

comp_vals <- c(
  A              = var_A,
  D              = var_D,
  AA             = var_AA,
  "Field spatial" = var_spatial,
  Residual       = var_resid
)

total_var <- sum(comp_vals)
props     <- comp_vals / total_var

MODEL_TAG <- paste0("noia_", KERNEL_MODE)

STW_DIR <- file.path(DATA_DIR, "Single_trait_walkthrough")
VC_DIR  <- file.path(STW_DIR, "h2_and_deregression_input")

out_df <- data.frame(
  Trait     = TRAIT,
  Year      = YEAR,
  Component = names(comp_vals),
  Var       = as.numeric(comp_vals),
  Prop      = as.numeric(props),
  Model     = MODEL_TAG,
  stringsAsFactors = FALSE)

out_file <- file.path(VC_DIR, sprintf("vc_%s_%d_%s_%s.csv", TRAIT, YEAR, MODEL_TAG, ifelse(Co_PCR, "Co", "pop")))


write.csv(out_df, out_file, row.names = FALSE)
cat("Saved VC CSV to: ", out_file, "\n")

wald.asreml(m_fit)
