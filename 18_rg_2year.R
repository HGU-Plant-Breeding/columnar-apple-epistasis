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
library(dplyr)
library(readr)

# ─────────────────────────────────────────────────────────────
# 2) Paths / switches — EDIT HERE
# ─────────────────────────────────────────────────────────────
DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"
SNP_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"

DIR_STW  <- file.path(DATA_DIR, "Single_trait_walkthrough")
IN_CLEAN <- file.path(DIR_STW, "outlier_removed_raw")

# Choose ONE trait here, then rerun whole script
TRAIT <- "in_len"   # "in_len", "lsh_ang_1", "tr_brz_prop", "tr_coni", "tr_len"

meta_cols <- c("Genotype","pop","Loc","Row","Column","rep","Plot_RC_ID","co_loc")

# ─────────────────────────────────────────────────────────────
# 3) Helper: load ONE clean single-trait file (trait + year)
# ─────────────────────────────────────────────────────────────
load_clean_trait_year <- function(trait, year) {
  pattern   <- sprintf("^clean_%s_%d_pop\\.csv$", trait, year)
  phenopath <- list.files(IN_CLEAN, pattern = pattern, full.names = TRUE)
  stopifnot(length(phenopath) == 1, file.exists(phenopath))
  
  df <- read.csv(phenopath, sep = ",") %>%
    mutate(
      Genotype   = as.character(Genotype),
      pop        = factor(pop),
      Loc   = factor(Loc),
      Row        = factor(Row),
      Column     = factor(Column),
      rep        = factor(rep),
      Plot_RC_ID = factor(Plot_RC_ID),
      co_loc = factor(co_loc)
    )
  
  df[[trait]] <- as.numeric(df[[trait]])
  
  df <- df %>%
    dplyr::select(any_of(c(meta_cols, trait)))
  
  return(df)
}

# ─────────────────────────────────────────────────────────────
# 4) Load 2017 + 2018 data for THIS TRAIT
# ─────────────────────────────────────────────────────────────

df2017 <- load_clean_trait_year(TRAIT, 2017)
df2018 <- load_clean_trait_year(TRAIT, 2018)

all_genotypes <- sort(unique(c(df2017$Genotype, df2018$Genotype)))
df2017$Genotype <- factor(df2017$Genotype, levels = all_genotypes)
df2018$Genotype <- factor(df2018$Genotype, levels = all_genotypes)

# ─────────────────────────────────────────────────────────────
# 5) Load precomputed Ginv
# ─────────────────────────────────────────────────────────────
Ginv <- readRDS(file.path(SNP_DIR, "MASTER_Ginv_flt_imp", "Ginv.rds"))

# ─────────────────────────────────────────────────────────────
# 6) Build bivariate (2-year) data set MANUALLY
# ─────────────────────────────────────────────────────────────
d17 <- df2017 %>% mutate(
  Year       = 2017,
  Genotype   = factor(Genotype, levels = all_genotypes),
  TY         = "Y2017",
  TYLoc      = interaction("Y2017", Loc, drop = TRUE))

d18 <- df2018 %>% mutate(
  Year       = 2018,
  Genotype   = factor(Genotype, levels = all_genotypes),
  TY         = "Y2018",
  TYLoc      = interaction("Y2018", Loc, drop = TRUE))

# stack
dat <- bind_rows(d17, d18)

# rename response column to 'value'
dat$value <- dat[[TRAIT]]

# standardise within year
dat <- dat %>%
  group_by(Year) %>%
  mutate(value = as.numeric(scale(value))) %>%
  ungroup()

# make factors used in model
dat$YearFactor <- factor(dat$TY, levels = c("Y2017", "Y2018"))
dat$Loc   <- factor(dat$Loc)
dat$Row        <- factor(dat$Row)
dat$Column     <- factor(dat$Column)
dat$rep        <- factor(dat$rep)
dat$Plot_RC_ID <- factor(dat$Plot_RC_ID)
dat$TYLoc      <- factor(dat$TYLoc)
dat <- dat %>% arrange(TYLoc, Plot_RC_ID, rep)
# ─────────────────────────────────────────────────────────────
# 7) Fit 2-year bivariate GBLUP for THIS TRAIT
# ─────────────────────────────────────────────────────────────
m_biv <- asreml(
  fixed    = value ~ YearFactor - 1 + YearFactor:pop,
  random   = ~ us(YearFactor):vm(Genotype, Ginv) +
    YearFactor:at(Loc):Row +
    YearFactor:at(Loc):Column,
  residual = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | TYLoc),
  data     = dat,
  na.action = na.method(y = "include", x = "include"),
  maxit    = 200,
  workspace = "512mb"
)

m_biv <- update(m_biv, aom = TRUE)

# ─────────────────────────────────────────────────────────────
# 8) Extract 2×2 G matrix and Rg
# ─────────────────────────────────────────────────────────────
vc   <- summary(m_biv)$varcomp
rn   <- rownames(vc)
comp <- vc[, "component"]
names(comp) <- rn

# Expected ASReml rownames
name_g11  <- "YearFactor:vm(Genotype, Ginv)!YearFactor_Y2017:Y2017"
name_g22  <- "YearFactor:vm(Genotype, Ginv)!YearFactor_Y2018:Y2018"
name_g12a <- "YearFactor:vm(Genotype, Ginv)!YearFactor_Y2017:Y2018"
name_g12b <- "YearFactor:vm(Genotype, Ginv)!YearFactor_Y2018:Y2017"

g11 <- as.numeric(comp[name_g11])
g22 <- as.numeric(comp[name_g22])

if (!is.na(comp[name_g12a])) {
  g12 <- as.numeric(comp[name_g12a])
} else {
  g12 <- as.numeric(comp[name_g12b])
}

Gmat <- matrix(
  c(g11, g12,
    g12, g22),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Y2017", "Y2018"), c("Y2017", "Y2018"))
)

sd_g <- sqrt(diag(Gmat))
Rg   <- Gmat / (sd_g %o% sd_g)
diag(Rg) <- 1

rg_2017_2018 <- Rg["Y2017", "Y2018"]

cat("\n## Across-year genetic correlation (2017 vs 2018) ##\n")
print(
  data.frame(
    Trait        = TRAIT,
    rg_2017_2018 = round(rg_2017_2018, 3),
    row.names    = NULL
  )
)

# ─────────────────────────────────────────────────────────────
# 9) Effective N = genotypes with data in BOTH years (for 95% CI)
# ─────────────────────────────────────────────────────────────
ids17 <- df2017 %>%
  filter(!is.na(.data[[TRAIT]])) %>%
  pull(Genotype) %>%
  unique()

ids18 <- df2018 %>%
  filter(!is.na(.data[[TRAIT]])) %>%
  pull(Genotype) %>%
  unique()

N_eff_current <- length(intersect(ids17, ids18))

# one-row result for this trait
res_row <- data.frame(
  Trait         = TRAIT,
  rg_2017_2018  = rg_2017_2018,
  N_eff         = N_eff_current,
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────
# 10) Cumulative CSV saving
# ─────────────────────────────────────────────────────────────
OUT_DIR  <- file.path(DIR_STW, "MV5_Rg_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(OUT_DIR, "rg_2017_2018_by_trait.csv")

if (file.exists(out_file)) {
  old <- read.csv(out_file, check.names = FALSE, stringsAsFactors = FALSE)
  old <- old[old$Trait != TRAIT, , drop = FALSE]
  rg_year_table <- rbind(old, res_row)
} else {
  rg_year_table <- res_row
}

rg_year_table$Trait <- factor(
  rg_year_table$Trait,
  levels = sort(unique(as.character(rg_year_table$Trait)))
)
rg_year_table <- rg_year_table[order(rg_year_table$Trait), ]
rg_year_table$Trait <- as.character(rg_year_table$Trait)

write.csv(rg_year_table, file = out_file, row.names = FALSE)

cat("\nBivariate rg CSV updated. File written to:\n", out_file, "\n")
