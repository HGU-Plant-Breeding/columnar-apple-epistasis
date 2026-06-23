# ──────────────────────────────────────────────────────────────────────────────
# 0) Reset session
rm(list=ls(all=TRUE)); graphics.off(); closeAllConnections()

# ──────────────────────────────────────────────────────────────────────────────
# 1) Libraries
  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(ASRgenomics)
  library(AGHmatrix)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)

# ──────────────────────────────────────────────────────────────────────────────
# 2) Paths / switches — EDIT HERE
DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"
SNP_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"

TRAIT <- "in_len"   # in_len, lsh_ang_1, tr_brz_prop, tr_coni, tr_len
YEAR  <- 2017          # 2017 or 2018

# Folders
dir_stw   <- file.path(DATA_DIR, "Single_trait_walkthrough")
in_clean  <- file.path(dir_stw, "outlier_removed_raw")                   # ← input folder
out_dir   <- file.path(dir_stw, "h2_and_deregression_input")         # ← output folder
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Load the single cleaned file: clean_{TRAIT}_{YEAR}_pop.csv (assume exactly one)
phenopath <- list.files(in_clean, pattern = sprintf("^clean_%s_%d_pop\\.csv$", TRAIT, YEAR), full.names = TRUE)[1]
cat("Loading ", basename(phenopath), "\n", sep = "")

meta_cols <- c("Genotype","pop","Loc","Row","Column","rep","Plot_RC_ID", "co_loc")
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

# ──────────────────────────────────────────────────────────────────────────────
# 4) Load precomputed Ginv
Ginv <- readRDS("C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip/MASTER_Ginv_flt_imp/Ginv.rds")

# ──────────────────────────────────────────────────────────────────────────────
# 5) Fit single-year GBLUP (pop fixed)
cat("\nFitting single-year GBLUP for ", TRAIT, " (", YEAR, ") — Fixed: pop", sep = "")
m_fit <- asreml(
  fixed     = as.formula(sprintf("%s ~ pop", TRAIT)),
  random    = ~ vm(Genotype, Ginv) + at(Loc):Row + at(Loc):Column,
  residual  = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | Loc),
  data      = df,
  na.action = na.method(y = "include", x = "include"),
  ai.sing   = FALSE
)
m_fit <- update(m_fit, aom = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 6) Metrics, varcomps, Wald
s   <- summary(m_fit)
vc  <- s$varcomp
LL  <- s$loglik
AIC <- s$aic
cat(sprintf("\nLogLik = %.3f | AIC = %.3f\n", LL, AIC))
print(vc)
print(wald.asreml(m_fit))

# ──────────────────────────────────────────────────────────────────────────────
# 7) Cullis h2 (single year)
Gvar  <- as.numeric(vc["vm(Genotype, Ginv)", "component"])   # σ²_A
predG <- predict(m_fit, classify = "Genotype", sed = TRUE)

# AVSED and h2 (use the 'mean' entry of avsed)
AVSED     <- as.numeric(predG$avsed["mean"])
h2_cullis <- (1 - (AVSED^2) / (2 * Gvar))

cat(sprintf("Gvar(vm[Genotype]): %.6f\nAVSED: %.6f\nh2_Cullis: %.6f\n", Gvar, AVSED, h2_cullis))

# ──────────────────────────────────────────────────────────────────────────────
# 9) BLUPs+SE with AVSED and h2 columns, merged with population BLUEs
# ──────────────────────────────────────────────────────────────────────────────
# (a) Get genotype-level predictions (BLUPs)
BLUPs_and_SEs <- predG$pvals |> mutate(AVSED = AVSED, h2_Cullis = h2_cullis, Gvar = Gvar)

# b) Blank BLUPs for genotypes with no phenotypic records for this trait×year
no_obs_genos <- df |>
  group_by(Genotype) |>
  summarise(all_na = all(is.na(.data[[TRAIT]])), .groups = "drop") |>
  filter(all_na) |>
  pull(Genotype) |>
  as.character()

if (length(no_obs_genos) > 0) {BLUPs_and_SEs <- BLUPs_and_SEs |> mutate(
  predicted.value = ifelse(Genotype %in% no_obs_genos, NA_real_, predicted.value),
  std.error = ifelse(Genotype %in% no_obs_genos & "std.error" %in% names(BLUPs_and_SEs), NA_real_, std.error))}

# (c) Add co_loc and pop from the original cleaned data
BLUPs_and_SEs <- BLUPs_and_SEs |> left_join(df |> dplyr::select(Genotype, pop, co_loc) |> distinct(), by = "Genotype")

# (d) Compute population-level BLUEs
pop_pred <- predict(m_fit, classify = "pop", sed = TRUE)
pop_blues <- pop_pred$pvals |> dplyr::transmute(pop, pop_BLUE = predicted.value)

# (d) Add pop_BLUE (merge by population)
BLUPs_and_SEs <- BLUPs_and_SEs |> left_join(pop_blues, by = "pop")

# (e) Save combined table
outfile_blup <- file.path(out_dir, sprintf("BLUPs_SEs_h2_%s_%d_pop.csv", TRAIT, YEAR))
write.csv(BLUPs_and_SEs, outfile_blup, row.names = FALSE, quote = TRUE)

