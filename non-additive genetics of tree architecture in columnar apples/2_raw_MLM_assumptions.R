# ──────────────────────────────────────────────────────────────────────────────
# 0) Reset session
rm(list=ls(all=TRUE));
graphics.off(); closeAllConnections()

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
  library(MASS)
})

# ──────────────────────────────────────────────────────────────────────────────
# 2) Paths / switches — EDIT HERE
DATA_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
SNP_DIR   <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"

# Switches
YEAR      <- 2017                 # 2017 or 2018
TRAIT     <- "in_len"        # tr_len, tr_coni, tr_brz_dens, tr_brz_prop, ssh_prop, lsh_ang_1, tsh_count, in_len
Co_PCR    <- FALSE                 # include co_loc in fixed effects?

FILE_17 <- file.path(DATA_DIR, "df_17_final_w_outliers.csv")
FILE_18 <- file.path(DATA_DIR, "df_18_final_w_outliers.csv")

# ──────────────────────────────────────────────────────────────────────────────
# 3) Load selected year (with outliers) & prepare
infile <- if (YEAR == 2017) FILE_17 else FILE_18
if (!file.exists(infile)) stop("Missing file: ", infile)

df <- read.csv(infile, sep=";") %>%
  mutate(
    Genotype   = factor(Genotype),
    pop        = factor(pop),
    Loc   = factor(Loc),
    Row        = factor(Row),
    Column     = factor(Column),
    rep        = factor(rep),
    Plot_RC_ID = factor(Plot_RC_ID),
    co_loc = factor(co_loc),
    Year       = YEAR
  )

# ──────────────────────────────────────────────────────────────────────────────
# 4) Build RMs + maf-filter for genotypes present this year   (+remove LD window)
c130 <- read.csv(file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv"), sep=",")
rownames(c130) <- c130$X
c130 <- c130[rownames(c130) %in% df$Genotype, , drop=FALSE]
M <- as.matrix(c130[,-1, drop=FALSE])
colnames(M) <- gsub("^AX\\.", "AX-", colnames(M))

# Remove LD window if Co_PCR = TRUE
if (Co_PCR) {
  GM  <- read.csv(file.path(SNP_DIR, "GM_v11.csv"), sep=";")
  win <- GM$chrom == 10 & GM$pos >= 2e7 & GM$pos <= 3.6e7
  hap <- unique(GM$name[win])
  M   <- M[, !colnames(M) %in% hap, drop=FALSE]
}

p_alt <- colMeans(M, na.rm = TRUE) / 2
maf   <- pmin(p_alt, 1 - p_alt)
keep <- which(maf >= 0.05)
M <- M[, keep, drop = FALSE]

for (j in which(colSums(is.na(M))>0)) M[is.na(M[,j]), j] <- mean(M[,j], na.rm=TRUE)

###  van Raden GRM  ###
grm_res <- snpReady::G.matrix(M=M, method="VanRaden", format="sparse")
A       <- grm_res$Ga
diag(A) <- diag(A) + 1e-2
Ginv    <- G.inverse(A, sparseform=TRUE)$Ginv

# ──────────────────────────────────────────────────────────────────────────────
# 5) Fit single-trait model
fixed_term <- as.formula(paste(TRAIT, "~ pop", if (Co_PCR) " + co_loc" else ""))

cat("\nFitting model for:", TRAIT, "×", YEAR,
    if (Co_PCR) "(+ Co fixed)" else "(no Co)")
m_fit <- asreml(
  fixed    = fixed_term,
  random   = ~ vm(Genotype, Ginv) + at(Loc):Row + at(Loc):Column,
  residual = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | Loc),
  data     = df,
  na.action= na.method(y="include", x="include"),
  ai.sing  = FALSE
)

# update to stabilize aom/residuals
m_fit <- update(m_fit, aom = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 6) check model metrics and varcomps
s <- summary(m_fit)
LogLik <- s$loglik; AIC <- s$aic
s$sigma

cat("\nModel metrics — ", TRAIT, "×", YEAR,
    if (Co_PCR) " (+ Co fixed)" else " (no Co)",
    "\nLogLik = ", LogLik, ", AIC = ", AIC, "\n", sep = "")

s$varcomp
wald.asreml(m_fit)

# 7) asreml default diagnostic plots
plot(m_fit)


