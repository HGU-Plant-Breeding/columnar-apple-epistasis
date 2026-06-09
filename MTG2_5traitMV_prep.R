rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

library(dplyr)

#####################################################################
# 0) SWITCHES & PATHS
#####################################################################

YEAR <- 2017   # <<< EDIT: 2017 or 2018

BASE_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory"
SNP_DIR   <- file.path(BASE_DIR, "SNP Chip")
DER_DIR   <- file.path(BASE_DIR, "Lidar", "Single_trait_walkthrough", "deregression_output")
MTG2_PROJ <- "C:/Users/nguevenc/mtg2_project"
dir.create(MTG2_PROJ, showWarnings = FALSE, recursive = TRUE)

#####################################################################
# 1) READ SNP FILE AND GET SNP IDS
#####################################################################

gt_file <- file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv")

gt <- read.csv(gt_file, sep = ",", check.names = FALSE)
colnames(gt)[1] <- "Genotype"
ids_snp <- as.character(gt$Genotype)

#####################################################################
# 2) READ EACH TRAIT FILE FOR THIS YEAR (SEPARATELY, ALPHABETICAL)
#####################################################################

# 1) in_len
f_in_len <- file.path(DER_DIR, sprintf("dGEBV_in_len_%d_pop.csv", YEAR))
in_len_raw <- read.csv(f_in_len, sep = ",", check.names = FALSE)

# 2) lsh_ang_1
f_lsh_ang_1 <- file.path(DER_DIR, sprintf("dGEBV_lsh_ang_1_%d_pop.csv", YEAR))
lsh_ang_1_raw <- read.csv(f_lsh_ang_1, sep = ",", check.names = FALSE)

# 3) tr_brz_prop
f_tr_brz_prop <- file.path(DER_DIR, sprintf("dGEBV_tr_brz_prop_%d_pop.csv", YEAR))
tr_brz_prop_raw <- read.csv(f_tr_brz_prop, sep = ",", check.names = FALSE)

# 4) tr_coni
f_tr_coni <- file.path(DER_DIR, sprintf("dGEBV_tr_coni_%d_pop.csv", YEAR))
tr_coni_raw <- read.csv(f_tr_coni, sep = ",", check.names = FALSE)

# 5) tr_len
f_tr_len <- file.path(DER_DIR, sprintf("dGEBV_tr_len_%d_pop.csv", YEAR))
tr_len_raw <- read.csv(f_tr_len, sep = ",", check.names = FALSE)

#####################################################################
# 3) CHECK REQUIRED COLUMNS
#####################################################################

needed_cols <- c("Genotype", "dGEBV", "weight")

#####################################################################
# 4) FIND COMMON GENOTYPE IDS ACROSS SNP + ALL 5 TRAITS
#####################################################################

ids_common <- Reduce(
  intersect,
  list(ids_snp,
       in_len_raw$Genotype,
       lsh_ang_1_raw$Genotype,
       tr_brz_prop_raw$Genotype,
       tr_coni_raw$Genotype,
       tr_len_raw$Genotype)
)
ids_common <- sort(ids_common)

# Sanity check: ordering must match
stopifnot(all(in_len_raw$Genotype      == ids_common),
          all(lsh_ang_1_raw$Genotype   == ids_common),
          all(tr_brz_prop_raw$Genotype == ids_common),
          all(tr_coni_raw$Genotype     == ids_common),
          all(tr_len_raw$Genotype      == ids_common))

### save ids_common to a text file
ids_file <- file.path(MTG2_PROJ, "ids_common.txt")
writeLines(ids_common, ids_file)

### save SNP matrix filtered to these IDs (for GRM prep later)
gt_common <- gt %>%
  dplyr::filter(Genotype %in% ids_common) %>%
  dplyr::arrange(match(Genotype, ids_common))

stopifnot(all(as.character(gt_common$Genotype) == ids_common))

gt_common_path <- file.path(MTG2_PROJ, "GT_filtered_ids_common.csv")
write.csv(gt_common, gt_common_path, row.names = FALSE, quote = FALSE)

#####################################################################
# 5) Z-STANDARDISE dGEBVs AND EXTRACT WEIGHTS (PER TRAIT)
#####################################################################

# dGEBV z-scores
in_len_z      <- scale(in_len_raw$dGEBV)[, 1]
lsh_ang_1_z   <- scale(lsh_ang_1_raw$dGEBV)[, 1]
tr_brz_prop_z <- scale(tr_brz_prop_raw$dGEBV)[, 1]
tr_coni_z     <- scale(tr_coni_raw$dGEBV)[, 1]
tr_len_z      <- scale(tr_len_raw$dGEBV)[, 1]

# weights
in_len_w      <- in_len_raw$weight
lsh_ang_1_w   <- lsh_ang_1_raw$weight
tr_brz_prop_w <- tr_brz_prop_raw$weight
tr_coni_w     <- tr_coni_raw$weight
tr_len_w      <- tr_len_raw$weight

#####################################################################
# 6) BUILD PHENOTYPE AND WEIGHT MATRICES (ALPHABETICAL ORDER)
#####################################################################

pheno_out <- data.frame(
  FID        = ids_common,
  IID        = ids_common,
  in_len     = in_len_z,
  lsh_ang_1  = lsh_ang_1_z,
  tr_brz_prop= tr_brz_prop_z,
  tr_coni    = tr_coni_z,
  tr_len     = tr_len_z
)

weights_out <- data.frame(
  FID        = ids_common,
  IID        = ids_common,
  in_len     = in_len_w,
  lsh_ang_1  = lsh_ang_1_w,
  tr_brz_prop= tr_brz_prop_w,
  tr_coni    = tr_coni_w,
  tr_len     = tr_len_w
)

#####################################################################
# 7) WRITE allpheno_YEAR_std.dat AND weights_YEAR.dat
#####################################################################

allpheno_path <- file.path(MTG2_PROJ, sprintf("allpheno_%d_std.dat", YEAR))
weights_path  <- file.path(MTG2_PROJ, sprintf("weights_%d.dat", YEAR))

write.table(pheno_out, allpheno_path,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE, na = "")

write.table(weights_out, weights_path,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE, na = "")
