rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

library(asreml); vm <- get("asr_vm", envir = asNamespace("asreml"))
library(data.table)
library(dplyr)
library(snpReady)
library(ASRgenomics)

# ── SWITCHES ───────────────────────────────────────────────────────────────
TRAIT   <- "in_len"   # in_len, lsh_ang_1, tr_brz_prop, tr_coni, tr_len
YEAR    <- 2017    # 2017, 2018
thr_fdr <- 0.05
Co_PCR  <- FALSE   # TRUE, FALSE where TRUE = fix co_loc + remove Co window

trait_map <- c(
  in_len      = "INL",
  tr_len      = "TL",
  tr_coni     = "TC",
  lsh_ang_1   = "LSA",
  tr_brz_prop = "TBZP")
TRAIT_GWAS <- trait_map[[TRAIT]]

# ── PATHS ───────────────────────────────────────────────────────────────────
BASE   <- "C:/Users/nguevenc/Desktop/R_Working_Directory"
STW    <- file.path(BASE, "Lidar", "Single_trait_walkthrough")
SNPDIR <- file.path(BASE, "SNP Chip")

pheno_file <- file.path(STW, "outlier_removed_raw", sprintf("clean_%s_%d_pop.csv", TRAIT, YEAR))
gwas_file <- file.path(STW,sprintf("BLINK_dGEBV_%d", YEAR),sprintf("GAPIT.Association.GWAS_Results.BLINK.%s(NYC).csv", TRAIT_GWAS))
GT_path <- file.path(SNPDIR, "GT_filtered_numeric_transposed.csv")
GM_path <- file.path(SNPDIR, "GM_v11.csv")
LD_path <- file.path(SNPDIR, "LD.txt")
VC_DIR <- file.path(STW, "sign_Kernel_VCs")

# ── PHENO + GWAS HITS ───────────────────────────────────────────────────────
df <- read.csv(pheno_file, check.names = FALSE)
gwas <- read.csv(gwas_file, check.names = FALSE) |> dplyr::select(SNP, Chr, Pos, P.value)
gwas$sgnf <- p.adjust(gwas$P.value, method = "BH")
hits <- gwas[gwas$sgnf < thr_fdr, ]
sign_SNPs <- hits$SNP

# ── GT/GM/LD ────────────────────────────────────────────────────────────────
GT <- fread(GT_path)
GM <- fread(GM_path)
LD <- fread(LD_path)
setnames(GT, 1, "Genotype")
M <- as.matrix(GT[, -1, with = FALSE])
rownames(M) <- GT$Genotype
colnames(M) <- gsub("^AX\\.", "AX-", colnames(M))

ids <- unique(df$Genotype[!is.na(df$Genotype)])
M <- M[ids, , drop = FALSE]
p_alt <- colMeans(M, na.rm = TRUE) / 2
maf   <- pmin(p_alt, 1 - p_alt)
M <- M[, which(maf >= 0.05), drop = FALSE]
for (j in which(colSums(is.na(M)) > 0))
M[is.na(M[, j]), j] <- mean(M[, j], na.rm = TRUE)
GM <- GM[name %in% colnames(M)][order(chrom, pos)]
sign_SNPs <- sign_SNPs[sign_SNPs %in% colnames(M)]

# ── CO WINDOW + LD EXPANSION FOR SIGNAL SET ──────────────────────────────────
co_window <- GM[chrom == 10L & pos >= 2e7 & pos <= 3.6e7, unique(name)]

LD <- LD[LD >= 0.8 & Chrom %in% hits$Chr & (Name1 %in% sign_SNPs | Name2 %in% sign_SNPs) & Name1 %in% GM$name & Name2 %in% GM$name]

ld_hits <- rbind(
  LD[Name1 %in% sign_SNPs, .(Query = Name1, Partner = Name2, Chrom, LD)],
  LD[Name2 %in% sign_SNPs, .(Query = Name2, Partner = Name1, Chrom, LD)])

pad_bp <- 2e4L
snps_ld <- character(0)

if (nrow(ld_hits)) {
  ld_hits[, Pos_query   := GM$pos[match(Query,   GM$name)]]
  ld_hits[, Pos_partner := GM$pos[match(Partner, GM$name)]]
  
  ld_windows <- ld_hits[, .(
    Chrom = unique(Chrom),
    win_start = pmax(0L, pmin(unique(Pos_query), min(Pos_partner)) - pad_bp),
    win_end   = pmax(unique(Pos_query), max(Pos_partner)) + pad_bp
  ), by = Query]
  
snps_ld <- GM |> inner_join(ld_windows, by = c("chrom" = "Chrom"), relationship = "many-to-many") |> filter(pos >= win_start & pos <= win_end) |> pull(name) |> unique()}
S_sig <- unique(c(sign_SNPs, snps_ld))

co_drop <- if (Co_PCR) co_window else character(0)
S_sig   <- setdiff(S_sig, co_drop)

M_sig <- M[, S_sig, drop = FALSE]
M_bg  <- M[, setdiff(colnames(M), c(S_sig, co_drop)), drop = FALSE]

# ── GRMs + INVERSES ──────────────────────────────────────────────────────────
A_sig <- snpReady::G.matrix(M = M_sig, method = "VanRaden", format = "wide")$Ga
A_bg  <- snpReady::G.matrix(M = M_bg,  method = "VanRaden", format = "wide")$Ga

diag(A_sig) <- diag(A_sig) + 1e-2
diag(A_bg)  <- diag(A_bg)  + 1e-2

A_sig <- A_sig / mean(diag(A_sig))
A_bg  <- A_bg  / mean(diag(A_bg))

Ginv_sig <- ASRgenomics::G.inverse(A_sig, sparseform = FALSE)$Ginv
Ginv_bg  <- ASRgenomics::G.inverse(A_bg,  sparseform = FALSE)$Ginv

df$Genotype   <- factor(df$Genotype, levels = rownames(Ginv_bg))
df$pop        <- factor(df$pop)
df$co_loc <- factor(df$co_loc)
df$Loc   <- factor(df$Loc)
df$Row        <- factor(df$Row)
df$Column     <- factor(df$Column)
df$Plot_RC_ID <- factor(df$Plot_RC_ID)
df$rep        <- factor(df$rep)

# ── FIT DOUBLE-KERNEL ─────────────────────────────────────────────────────────
fixed_formula <- if (Co_PCR) as.formula(paste0(TRAIT, " ~ pop + co_loc")) else as.formula(paste0(TRAIT, " ~ pop"))

m_fit <- asreml(
  fixed     = fixed_formula,
  random    = ~ vm(Genotype, Ginv_sig) + vm(Genotype, Ginv_bg) + at(Loc):Row + at(Loc):Column,
  residual  = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | Loc),
  data      = df,
  na.action = na.method(y = "include", x = "include"),
  ai.sing   = FALSE,
  workspace = "512mb",
  pworkspace = "512mb"
)
m_fit <- update(m_fit, aom = TRUE)

vc <- summary(m_fit)$varcomp
print(vc)

# ── SAVE VC PROPORTIONS (AVERAGE across Locs for spatial/residual) ───────
V_sig <- as.numeric(vc["vm(Genotype, Ginv_sig)", "component"])
V_bg  <- as.numeric(vc["vm(Genotype, Ginv_bg)",  "component"])
V_sp <- mean(vc[grep("^at\\(Loc, '.+'\\):(Row|Column)$", rownames(vc)), "component"],na.rm = TRUE)
V_res <- mean(vc[grep("^Loc_.+!R$", rownames(vc)), "component"],na.rm = TRUE)
V_tot <- V_sig + V_bg + V_sp + V_res

out <- data.table(
  Trait     = TRAIT,
  Year      = YEAR,
  Model     = if (Co_PCR) "pop + Co" else "pop",
  Component = c("A_sig","A_bg","Field spatial","Residual"),
  Var       = c(V_sig, V_bg, V_sp, V_res),
  Prop  = 100 * c(V_sig, V_bg, V_sp, V_res) / V_tot
)

outfile <- file.path(
  VC_DIR,
  sprintf("vc_signKernel_%s_%d_%s.csv", TRAIT, YEAR, if (Co_PCR) "pop + Co" else "pop")
)
fwrite(out, outfile)

