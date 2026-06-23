rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(dplyr)
  library(data.table)
  library(stringr)
  library(snpReady)
  library(ASRgenomics)

# ── 0) SWITCHES ───────────────────────────────────────────────────────────────
  TRAIT    <- "INL"    # "INL", "LSA", "TBZP", "TC", "TL"
  YEAR     <- 2017      # 2017 or 2018
  threshold <- 0.05     # FDR threshold for GWAS hits , modify under step 3

  # ── Trait legend ─────────
  trait_file_map <- c(
    "INL"  = "in_len",
    "TBZP" = "tr_brz_prop",
    "TL"   = "tr_len",
    "TC"   = "tr_coni",
    "LSA"  = "lsh_ang_1"
  )
  
  TRAIT_FILE <- trait_file_map[TRAIT]
  
  
# ── 1) PATHS ──────────────────────────────────────────────────────────────────
BASE_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory"

LIDAR_DIR  <- file.path(BASE_DIR, "Lidar")
STW_DIR    <- file.path(LIDAR_DIR, "Single_trait_walkthrough")
SNP_DIR    <- file.path(BASE_DIR, "SNP Chip")
GT_path    <- file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv")
GM_path    <- file.path(SNP_DIR, "GM_v11.csv")
LD_path    <- file.path(SNP_DIR, "LD.txt")
BLINK_DIR <- file.path(STW_DIR, sprintf("BLINK_dGEBV_%d", YEAR))
# Trait-specific GWAS NYC file (BLINK)
gwas_file <- file.path(BLINK_DIR,sprintf("GAPIT.Association.GWAS_Results.BLINK.%s(NYC).csv", TRAIT))
pheno_dir  <- file.path(STW_DIR, "outlier_removed_raw")
pheno_file <- file.path(pheno_dir, sprintf("clean_%s_%d_pop.csv", TRAIT_FILE, YEAR))

# ── 2) READ PHENOTYPE ────────────────────────────────────────────────────────
df <- read.csv(pheno_file, check.names = FALSE)
names(df)[names(df) == TRAIT_FILE] <- TRAIT


# ── 3) READ TRAIT-SPECIFIC GWAS FILE & APPLY FDR  ─────────────────────────────
gwas <- read.csv(gwas_file, check.names = FALSE)
gwas <- gwas %>% dplyr::select(SNP, Chr, Pos, P.value)
gwas$sgnf <- p.adjust(gwas$P.value, method = "BH")
gwas_hits <- gwas %>% filter(sgnf < threshold)
n_hits    <- nrow(gwas_hits)
gwas_hits; n_hits

# ── 4) READ GT, GM, LD (GT + LD via fread) ───────────────────────────────────
GT <- fread(GT_path)
GM <- fread(GM_path)
LD <- fread(LD_path)

# GT: first column = genotype IDs (e.g. "V1")
colnames(GT)[1] <- "Genotype"
M      <- as.matrix(GT[, -1, drop = FALSE])
rownames(M) <- GT$Genotype
colnames(M) <- gsub("^AX\\.", "AX-", colnames(M))

# ── 5) FILTER GT GENOTYPES, MAF ≥ 0.05, MEAN-IMPUTE ─────────────────
ids_nonNA <- unique(df$Genotype[!is.na(df$Genotype)])
M <- M[ids_nonNA, , drop = FALSE]

# MAF filter (≥ 0.05)
p_alt    <- colMeans(M, na.rm = TRUE) / 2
maf      <- pmin(p_alt, 1 - p_alt)
keep_maf <- which(maf >= 0.05)
M <- M[, keep_maf, drop = FALSE]

# Mean-impute GT
for (j in which(colSums(is.na(M))>0)) M[is.na(M[,j]), j] <- mean(M[,j], na.rm=TRUE)

# Filter GM to markers present in M, sort by chrom, pos
GM <- GM[name %in% colnames(M)]
GM <- GM[order(chrom, pos)]

sign_SNPs       <- gwas_hits$SNP[gwas_hits$SNP %in% colnames(M)]

# ── 6) CO WINDOW (chr10: 20–36 Mb) ───────────────────────────────────
co_window_snps <- GM[chrom == 10L & pos >= 2e7 & pos <= 3.6e7, unique(name)]

# ── 8) LD ≥ 0.8 WINDOWS AROUND ALL sign_SNPs ────────────────────
LD <- LD[LD >= 0.8 & Chrom %in% gwas_hits$Chr 
            & (Name1 %in% gwas_hits$SNP | Name2 %in% gwas_hits$SNP) 
            & Name1 %in% GM$name & Name2 %in% GM$name]
  
  # 3) Long table: (Query = sign SNP, Partner = LD neighbour)
  ld_hits <- rbind(
    LD[Name1 %in% gwas_hits$SNP,.(Query = Name1, Partner = Name2, Chrom, LD)],
    LD[Name2 %in% gwas_hits$SNP,.(Query = Name2, Partner = Name1, Chrom, LD)])

  ld_hits$Pos_query <- GM$pos[match(ld_hits$Query,   GM$name)]
  ld_hits$Pos_partner <- GM$pos[match(ld_hits$Partner, GM$name)]

  # collapse to one LD block per SNP
  pad_bp <- 2e4L          # ±20 kb extra
  ld_windows <- ld_hits[, .(Pos_query   = unique(Pos_query),
                            min_partner = min(Pos_partner), max_partner = max(Pos_partner), 
                            Chrom = unique(Chrom)), by = Query ][, .(Query,Chrom, 
                            win_start = pmax(0L, pmin(Pos_query, min_partner) - pad_bp),
                            win_end = pmax(Pos_query, max_partner) + pad_bp)]
    
# ── 9) ADD sign_SNP DOSAGES TO PHENOTYPE ─────────────
  snp_cov <- as.data.frame(M[, sign_SNPs, drop = FALSE])
  snp_cov$Genotype <- rownames(M)
  df <- df %>% left_join(snp_cov, by = "Genotype")
  
  # make SNP column names syntactically safe (AX-... -> AX.)
  snp_model_names <- make.names(sign_SNPs)
  names(df)[match(sign_SNPs, names(df))] <- snp_model_names
  
  nrow(gwas_hits)
  
# ── 10) REMOVE sign_SNPs & LD WINDOWS (and Co WINDOW) FROM GRM MARKERS ──────
snps_to_drop <- unique(c(co_window_snps))

  snps_to_drop_from_ld <- GM %>%
    inner_join(ld_windows, by = c("chrom" = "Chrom"),
               relationship = "many-to-many") %>%
    filter(pos >= win_start & pos <= win_end) %>%
    pull(name)
  

# append these to Co-window's snps_to_drop
snps_to_drop <- unique(c(snps_to_drop, snps_to_drop_from_ld))

# add sign_SNPs to drop from GRM
snps_to_drop <- unique(c(snps_to_drop, sign_SNPs))
snps_to_drop <- snps_to_drop[snps_to_drop %in% colnames(M)]

M <- M[, setdiff(colnames(M), snps_to_drop), drop = FALSE]

cat("Markers remaining for GRM: ", ncol(M), "\n\n", sep = "")

# ── 11) BUILD VAN RADEN GRM & INVERSE ────────────────────────────────────────
grm_res <- snpReady::G.matrix(M = M, method = "VanRaden", format = "wide")
A       <- grm_res$Ga
diag(A) <- diag(A) + 1e-2
Ginv    <- ASRgenomics::G.inverse(A, sparseform = FALSE)$Ginv
g_levels <- rownames(Ginv)

# align df$Genotype with Ginv but keep full grid design
df$Genotype <- factor(df$Genotype, levels = g_levels)

# ── 12) SET FACTORS FOR ASREML ────────────────────────────────────────────────
df$pop         <- factor(df$pop)
df$co_loc  <- factor(df$co_loc)
df$Loc    <- factor(df$Loc)
df$Row         <- factor(df$Row)
df$Column      <- factor(df$Column)
df$Plot_RC_ID  <- factor(df$Plot_RC_ID)
df$rep         <- factor(df$rep)

# ── 13) FIT ASREML GBLUP MODEL ───────────────────────────────────────────────
fixed_terms <- c("pop", "co_loc")

  snp_terms    <- snp_model_names
  snp_co_terms <- paste0("co_loc:", snp_model_names)
  
    snp_pairs <- t(combn(snp_model_names, 2))
    snp_pair_terms <- apply(
      snp_pairs, 1,
      function(x) paste0(x[1], ":", x[2]))
    
  fixed_terms <- c(fixed_terms, snp_terms, snp_co_terms, snp_pair_terms)

fixed_rhs     <- paste(fixed_terms, collapse = " + ")
fixed_formula <- as.formula(paste0(TRAIT, " ~ ", fixed_rhs))

m_fit <- asreml(
  fixed     = fixed_formula,
  random    = ~ vm(Genotype, Ginv) +
    at(Loc):Row +
    at(Loc):Column,
  residual  = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | Loc),
  data      = df,
  na.action = na.method(y = "include", x = "include"),
  ai.sing   = FALSE,
  workspace = "512mb",
  pworkspace = "512mb"
)

m_fit <- update(m_fit, aom = TRUE)

cat("  LogLik: ", summary(m_fit)$loglik, "\n", sep = "")
cat("  AIC   : ", summary(m_fit)$aic,    "\n\n", sep = "")

# ── 15) WALD TEST, ADD FDR, WRITE CSV ────────────────────────────────────────
W_df  <- as.data.frame(wald.asreml(m_fit))
W_df$Effect <- rownames(W_df)
rownames(W_df) <- NULL
pcol <- "Pr(Chisq)"

# mark SNP-related rows (main SNPs, Co:SNP, SNP:SNP)
snps_idx <- grepl("AX\\.", W_df$Effect)

# BH FDR adjustment on SNP-related rows only
W_df$FDR_BH[snps_idx] <- p.adjust(W_df[[pcol]][snps_idx], method = "BH")
W_df

wald_outfile <- file.path(
  BLINK_DIR,
  sprintf("wald_%s_%d_SNPfixed_GBLUP.csv", TRAIT, YEAR)
)
write.csv(W_df, wald_outfile, row.names = FALSE)

cat("Wald table written to:\n  ", wald_outfile, "\n", sep = "")

