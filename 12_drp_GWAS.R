
rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

# ── 0) SWITCHES ──────────────────────────────────────────────
YEAR <- 2018   # 2017 or 2018

# ── 1) PATHS ─────────────────────────────────────────────────
BASE_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory"
SNP_DIR <- file.path(BASE_DIR, "SNP Chip")
DER_DIR <- file.path(BASE_DIR, "Lidar", "Single_trait_walkthrough", "deregression_output")
GM_path <- file.path(SNP_DIR, "GM_v11.csv")
GT_path <- file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv")

OUT_DIR <- file.path(BASE_DIR, "Lidar", "Single_trait_walkthrough", sprintf("BLINK_dGEBV_%d", YEAR))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

# ── 2) PACKAGES + GAPIT ─────────────────────────────────────
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(data.table)
  library(readr)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

# ─────────────────────────────────────────────────────────────
# 3) GENOTYPE PREP (GT -> M, Co-window removal, basic align)
# ─────────────────────────────────────────────────────────────
c130 <- read.csv(GT_path, sep = ",", check.names = FALSE)
rownames(c130) <- c130[, 1]
c130[, 1] <- NULL

M <- as.matrix(c130)
storage.mode(M) <- "numeric"
colnames(M) <- gsub("^AX\\.", "AX-", colnames(M))

# GM: marker annotation
GM <- read.csv(GM_path, sep = ";", stringsAsFactors = FALSE)
stopifnot(all(c("name", "chrom", "pos") %in% names(GM)))
GM$name <- gsub("^AX\\.", "AX-", GM$name)

# remove chr10 Co window (20–36 Mb) in both M and GM
win <- GM$chrom == 10 & GM$pos >= 2e7 & GM$pos <= 3.6e7
hap <- unique(GM$name[win])
M  <- M[, !colnames(M) %in% hap, drop = FALSE]
GM <- GM[!GM$name %in% hap, , drop = FALSE]

# basic intersection & alignment of SNPs (no MAF yet)
common_snps0 <- intersect(colnames(M), GM$name)
M  <- M[, common_snps0, drop = FALSE]
GM <- GM[match(common_snps0, GM$name), , drop = FALSE]
stopifnot(identical(colnames(M), GM$name))

# ─────────────────────────────────────────────────────────────
# 4) PHENOTYPE PREP (dGEBVs + pop_BLUE)
# ─────────────────────────────────────────────────────────────

files <- list.files(
  DER_DIR,
  pattern = "^dGEBV_.*_pop\\.csv$",
  full.names = TRUE
)

pheno_list <- lapply(files, function(f) {
  df <- read.csv(f, check.names = FALSE)
  
  # dGEBV := dGEBV + pop_BLUE (no new column)
  df$dGEBV <- as.numeric(df$dGEBV) + as.numeric(df$pop_BLUE)
  
  df %>%
    dplyr::transmute(
      Trait      = as.character(Trait),
      Year       = as.integer(Year),
      Genotype   = as.character(Genotype),
      dGEBV      = as.numeric(dGEBV),
      co_loc = co_loc
    )
})
names(pheno_list) <- basename(files)
d_all <- dplyr::bind_rows(pheno_list)

make_year_df <- function(df, yr) {
  d_year <- df %>% filter(Year == yr)
  wide <- d_year %>%
    pivot_wider(
      id_cols = c(Genotype, co_loc),
      names_from = Trait,
      values_from = dGEBV
    )
  
  # Rename traits safely
  rename_map <- c(
    "in_len"       = "INL",
    "tr_brz_prop"  = "TBZP",
    "tr_len"       = "TL",
    "tr_coni"      = "TC",
    "lsh_ang_1"    = "LSA"
  )
  
  names(wide) <- ifelse(names(wide) %in% names(rename_map),
                        rename_map[names(wide)],
                        names(wide))
  
  wide$co_loc <- factor(wide$co_loc)
  names(wide)[names(wide) == "Genotype"] <- "Taxa"
  wide
}

Y_2017 <- make_year_df(d_all, 2017)
Y_2018 <- make_year_df(d_all, 2018)

Y <- if (YEAR == 2017) Y_2017 else Y_2018
Y$Taxa <- as.character(Y$Taxa)

# ─────────────────────────────────────────────────────────────
# 5) Subset genotypes to 509 phenotyped IDs, then MAF + impute
# ─────────────────────────────────────────────────────────────
M <- M[Y$Taxa, , drop = FALSE]
rownames(M) <- Y$Taxa


p_alt <- colMeans(M, na.rm = TRUE) / 2
maf   <- pmin(p_alt, 1 - p_alt)
keep <- which(maf >= 0.05)
M <- M[, keep, drop = FALSE]

for (j in which(colSums(is.na(M))>0)) M[is.na(M[,j]), j] <- mean(M[,j], na.rm=TRUE)

GM <- GM[keep, , drop = FALSE]
stopifnot(identical(colnames(M), GM$name))

# ─────────────────────────────────────────────────────────────
# 6) Build GD
# ─────────────────────────────────────────────────────────────
names(GM) <- c("SNP", "Chromosome", "Position")

GD <- data.frame(Taxa = rownames(M), M,
                 check.names = FALSE, stringsAsFactors = FALSE)
GD$Taxa <- as.character(GD$Taxa)

stopifnot(identical(colnames(GD)[-1], GM$SNP))
stopifnot(identical(Y$Taxa, GD$Taxa))

CV <- data.frame(
  Taxa = Y$Taxa,
  co   = as.integer(as.character(Y$co_loc)),
  stringsAsFactors = FALSE
)

Y <- as.data.frame(Y, stringsAsFactors = FALSE)[,-2]

# ─────────────────────────────────────────────────────────────
# 7) Run GAPIT BLINK
# ─────────────────────────────────────────────────────────────
out <- GAPIT(
  Y           = Y,
  GD          = GD,
  GM          = GM,
  CV          = CV,
  model       = "BLINK",
  PCA.total   = 3,
  file.output = TRUE
)

cat("GAPIT BLINK completed.\n")


write.csv(CV, file.path(DER_DIR, "Co_CV.csv"), row.names = FALSE)
