# 0) Reset session
rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

# 1) Libraries
library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(dplyr)
  library(readr)
  library(tidyr)
  library(rlang) 

# 2) Paths / switches
DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"
SNP_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"

DIR_STW   <- file.path(DATA_DIR, "Single_trait_walkthrough")
IN_CLEAN  <- file.path(DIR_STW,   "outlier_removed_raw")

TRAITS <- c("tr_len", "in_len", "tr_coni", "tr_brz_prop", "lsh_ang_1")

## <<< YEAR SWITCH HERE >>>
YEAR <- 2017   # 2017, 2018
## <<< YEAR SWITCH HERE >>>

meta_cols <- c(
  "Genotype","pop","Loc","Row","Column","rep",
  "Plot_RC_ID","co_loc"
)

# 3) Helper: load single-trait file
load_clean_trait_year <- function(trait, year) {
  pattern <- sprintf("^clean_%s_%d_pop\\.csv$", trait, year)
  phenopath <- list.files(IN_CLEAN, pattern = pattern, full.names = TRUE)
  
  df <- read.csv(phenopath, sep = ",") %>%
    mutate(
      Genotype   = as.character(Genotype),
      pop        = factor(pop),
      Loc   = factor(Loc),
      Row        = factor(Row),
      Column     = factor(Column),
      rep        = factor(rep),
      Plot_RC_ID = factor(Plot_RC_ID),
      co_loc = factor(co_loc),
      !!sym(trait) := as.numeric(!!sym(trait)) 
    ) %>%
    dplyr::select(any_of(c(meta_cols, trait)))
  df
}

# 4) Build 5-trait data frame for the chosen YEAR
message(sprintf("\n--- Loading and Merging Data Files for YEAR = %d ---", YEAR))

df_year <- load_clean_trait_year(TRAITS[1], YEAR)

for (tr in TRAITS[-1]) {
  tmp <- load_clean_trait_year(tr, YEAR)
  df_year[[tr]] <- tmp[[tr]]
}

# 5) Load precomputed Ginv
Ginv <- readRDS(file.path(SNP_DIR, "MASTER_Ginv_flt_imp", "Ginv.rds"))
df_year$Genotype <- factor((as.character(df_year$Genotype)) )
  
# 6) Long format + scaling
unstack_and_scale <- function(df) {
  df_long <- df %>%
    as.data.frame() %>% 
    pivot_longer(
      cols = all_of(TRAITS),
      names_to = "TraitFactor",
      values_to = "value"
    ) %>%
    mutate(
      value = ave(value, TraitFactor, FUN = function(x) as.numeric(scale(x))),
      TraitFactor = factor(TraitFactor, levels = TRAITS)
    ) %>%
    as.data.frame()
  df_long
}

df_long <- unstack_and_scale(df_year)

# explicit Trait × Loc factor for residual nesting
df_long$TFLoc <- interaction(df_long$TraitFactor,
                             df_long$Loc,
                             drop = TRUE)

df_long <- df_long %>%
  arrange(TFLoc, Plot_RC_ID, rep)

# 7) Multivariate 5-trait model
message(sprintf(
  "\n--- 5-trait MV GBLUP for YEAR = %d ---",
  YEAR
))

m_mv5 <- asreml(
  fixed    = value ~ TraitFactor - 1 + TraitFactor:pop, 
  random   = ~ us(TraitFactor):vm(Genotype, Ginv) +
    TraitFactor:at(Loc):Row +
    TraitFactor:at(Loc):Column,
  residual = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | TFLoc),
  data     = df_long,
  na.action = na.method(y = "include", x = "include"),
  maxit = 200,
  workspace = "512mb"
)
m_mv5 <- update(m_mv5, aom = TRUE, data = df_long)


# 8) Extract G and Rg from 5-trait MV model
extract_G_Rg <- function(model, year_label) {
  message(sprintf("\n--- Extracting G and Rg for %s ---", year_label))
  
  vc <- summary(model)$varcomp
  rn <- rownames(vc)
  comp <- vc[, "component"]
  names(comp) <- rn
  
  n_traits <- length(TRAITS)
  Gmat <- matrix(NA_real_, n_traits, n_traits, dimnames = list(TRAITS, TRAITS))
  
  for (i in seq_along(TRAITS)) {
    for (j in seq_along(TRAITS)) {
      t1 <- TRAITS[i]
      t2 <- TRAITS[j]

      key <- sprintf("TraitFactor:vm(Genotype, Ginv)!TraitFactor_%s:%s", t1, t2)
      
      if (key %in% rn) {
        val <- as.numeric(comp[key])
        Gmat[t1, t2] <- val
      } else {
      }
    }
  }

  # enforce symmetry
  for (i in seq_along(TRAITS)) {
    for (j in seq_along(TRAITS)) {
      if (i != j) {
        if (is.na(Gmat[i, j]) && !is.na(Gmat[j, i])) {
          Gmat[i, j] <- Gmat[j, i]
        } else if (!is.na(Gmat[i, j]) && is.na(Gmat[j, i])) {
          Gmat[j, i] <- Gmat[i, j]
        }
      }
    }
  }

  # genetic correlation matrix
  sd_g <- sqrt(diag(Gmat))
  Rg   <- Gmat / (sd_g %o% sd_g)
  diag(Rg) <- 1
  
  list(
    Year = year_label,
    G    = Gmat,
    Rg   = Rg
  )
}

mv_res <- extract_G_Rg(m_mv5, YEAR)

cat(sprintf("\n## 📊 %d Genetic Correlation Matrix (Rg) ##\n", YEAR))
print(round(mv_res$Rg, 3))

# 9) Save Rg matrix for heatmaps
OUT_DIR  <- file.path(DIR_STW, "MV5_Rg_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(OUT_DIR, sprintf("Rg_5trait_%d.csv", YEAR))
write.csv(
  mv_res$Rg,
  file = out_file,
  row.names = TRUE
)