rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

library(dplyr)

#####################################################################
# 0) SWITCHES & PATHS
#####################################################################

BASE_DIR   <- "C:/Users/nguevenc/Desktop/R_Working_Directory"

SNP_DIR    <- file.path(BASE_DIR, "SNP Chip")  # no longer used directly, but kept for reference
MTG2_PROJ  <- "C:/Users/nguevenc/mtg2_project"
INPUT_DIR  <- file.path(MTG2_PROJ, "input")

PLINK_BIN  <- "C:/Users/nguevenc/plink_win64_20250819/plink.exe"
GCTA_BIN   <- "C:/Users/nguevenc/gcta-1.95.0-Win-x86_64/gcta-1.95.0-Win-x86_64/bin/gcta64.exe"

#####################################################################
# 1) READ FILTERED SNP GENOTYPES (COMMON ID SET)
#####################################################################

gt_file <- file.path(MTG2_PROJ, "GT_filtered_ids_common.csv")
gt <- read.csv(gt_file, sep = ",", check.names = FALSE)
colnames(gt)[1] <- "Genotype"
ids_keep <- as.character(gt$Genotype)

#####################################################################
# 2) SNP MATRIX
#####################################################################

M <- as.matrix(gt[, -1])   # SNPs only
n_ind <- nrow(M)
n_snp <- ncol(M)

#####################################################################
# 4) CREATE FAM FILE (apple.fam) FOR MTG2 / PLINK
#####################################################################

fam <- data.frame(
  FID    = ids_keep,
  IID    = ids_keep,
  FATHER = 0,
  MOTHER = 0,
  SEX    = 0,
  PHENO  = -9
)

fam_path <- file.path(INPUT_DIR, "apple.fam")

write.table(fam, fam_path,
            quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)

#####################################################################
# 5) CREATE PED / MAP FOR PLINK
#####################################################################

ped_base <- file.path(INPUT_DIR, "apple")  # apple.ped / apple.map

# Convert 0/1/2 -> allele pairs (A/T), 0/NA -> missing "0"
A1 <- ifelse(M == 0, "A",
             ifelse(M == 1, "A",
                    ifelse(M == 2, "T", "0")))

A2 <- ifelse(M == 0, "A",
             ifelse(M == 1, "T",
                    ifelse(M == 2, "T", "0")))

geno_plink <- matrix(NA_character_, nrow = n_ind, ncol = 2 * n_snp)
geno_plink[, seq(1, 2 * n_snp, by = 2)] <- A1
geno_plink[, seq(2, 2 * n_snp, by = 2)] <- A2

ped <- data.frame(
  FID   = ids_keep,
  IID   = ids_keep,
  PID   = 0,
  MID   = 0,
  SEX   = 0,
  PHENO = -9,
  geno_plink,
  check.names = FALSE
)

ped_path <- paste0(ped_base, ".ped")
map_path <- paste0(ped_base, ".map")

write.table(ped, ped_path,
            quote = FALSE, sep = " ",
            row.names = FALSE, col.names = FALSE)

# Dummy MAP file (needed for PLINK but not used)
snp_names <- colnames(gt)[-1]
map <- data.frame(
  CHR = 1,
  SNP = snp_names,
  CM  = 0,
  BP  = seq_len(n_snp),
  stringsAsFactors = FALSE
)

write.table(map, map_path,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

#####################################################################
# 6) RUN PLINK: PED/MAP -> BED/BIM/FAM
#####################################################################

plink_cmd <- sprintf('"%s" --file "%s" --make-bed --out "%s"',
                     PLINK_BIN, ped_base, ped_base)

plink_status <- system(plink_cmd)

#####################################################################
# 7) RUN GCTA: BED/BIM/FAM -> GRM
#####################################################################

grm_base <- file.path(INPUT_DIR, "apple_grm")

gcta_cmd <- sprintf('"%s" --bfile "%s" --make-grm --out "%s"',
                    GCTA_BIN, ped_base, grm_base)

gcta_status <- system(gcta_cmd)
