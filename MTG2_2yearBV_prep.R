rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

library(dplyr)
library(readr)

#####################################################################
# 0) SWITCHES
#####################################################################

TRAIT    <- "tr_len"   # tr_len, tr_coni, tr_brz_prop, in_len, lsh_ang_1

BASE_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory"
DER_DIR   <- file.path(BASE_DIR, "Lidar", "Single_trait_walkthrough", "deregression_output")
MTG2_PROJ <- "C:/Users/nguevenc/mtg2_project"

#####################################################################
# 1) READ dGEBV FILES FOR THIS TRAIT (BOTH YEARS)
#####################################################################

y17 <- file.path(DER_DIR, sprintf("dGEBV_%s_2017_pop.csv", TRAIT))
y18 <- file.path(DER_DIR, sprintf("dGEBV_%s_2018_pop.csv", TRAIT))

d1 <- read.csv(y17, sep = ",", check.names = FALSE)
d2 <- read.csv(y18, sep = ",", check.names = FALSE)

#####################################################################
# 2) USE PRE-SAVED COMMON ID LIST FROM 5-TRAIT PREP
#####################################################################

ids_common_file <- file.path(MTG2_PROJ, "ids_common.txt")
ids_common_all <- readLines(ids_common_file)

ids_y1 <- as.character(d1$Genotype)
ids_y2 <- as.character(d2$Genotype)

# Intersection of: IDs used in multivariate 5-trait + this trait in both years
ids_both <- Reduce(intersect, list(ids_common_all, ids_y1, ids_y2))
ids_both <- sort(ids_both)

#####################################################################
# 3) ORDER BOTH YEARS BY THE COMMON ID LIST
#####################################################################

d1 <- d1 %>%
  filter(Genotype %in% ids_both) %>%
  arrange(match(Genotype, ids_both))

d2 <- d2 %>%
  filter(Genotype %in% ids_both) %>%
  arrange(match(Genotype, ids_both))

stopifnot(all(as.character(d1$Genotype) == ids_both),
          all(as.character(d2$Genotype) == ids_both))

#####################################################################
# 4) Z-STANDARDISE dGEBVs PER YEAR (AFTER ID FILTERING)
#####################################################################

z1 <- scale(d1$dGEBV)[, 1]
z2 <- scale(d2$dGEBV)[, 1]

#####################################################################
# 5) BUILD rg_TRAIT.dat and rg_TRAIT.wt
#####################################################################

rg_dat <- data.frame(
  FID = ids_both,
  ID  = ids_both,
  y1  = z1,
  y2  = z2
)

rg_wt <- data.frame(
  FID = ids_both,
  ID  = ids_both,
  w1  = d1$weight,
  w2  = d2$weight
)

out_dat <- file.path(MTG2_PROJ, sprintf("rg_%s.dat", TRAIT))
out_wt  <- file.path(MTG2_PROJ, sprintf("rg_%s.wt",  TRAIT))

write.table(rg_dat, out_dat,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(rg_wt,  out_wt,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)