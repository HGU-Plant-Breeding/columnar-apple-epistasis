rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(dplyr)
})

FQ_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
SNP_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"
OUT_DIR <- file.path(FQ_DIR, "GEBV_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATASET <- "JA"     # "JA" or "FRGB"
TRAIT   <- "anthocyanin"

# JA traits:   anthocyanin, brix, acids, a, phenols
# FRGB traits: firmness, area_cm2, diameter_cm, h_w_ratio,
#              avg_R, avg_G, avg_B, circularity, CIE_a

# ── Load phenotype ────────────────────────────────────────────────────────────
in_file <- file.path(FQ_DIR, if (DATASET == "JA") "juice_antho_traits_OutRmv.csv"
                     else                  "firmness_rgb_OutRmv.csv")

df <- read.csv(in_file, sep = ";") %>%
  mutate(
    Genotype  = factor(Genotype),
    pop       = factor(pop),
    rootstock = factor(rootstock),
    location  = factor(location),
    year      = factor(year),
    row       = factor(row),
    column    = factor(column),
    loc_yr    = factor(paste(location, year, sep = "_"))
  )
if (DATASET == "FRGB") df <- df %>% mutate(rep = as.integer(rep))

if (TRAIT == "anthocyanin") {
  df$pop    <- factor(df$pop, levels = c("bibon", "pxw"))
  df        <- df[df$location != "fuchsberg", ]
  df$loc_yr <- droplevels(df$loc_yr)
}

# metadata before subsetting — geno_pop VOR bibon-Drop (bibon-Genotypen behalten pop-Info)
geno_pop <- df %>%
  filter(!is.na(Genotype), !is.na(pop)) %>%
  transmute(Genotype = as.character(Genotype),
            pop      = as.character(pop)) %>%
  distinct()

# drop bibon for JA BEFORE extracting geno_with_pheno
if (DATASET == "JA") {
  df        <- df[df$location != "bibon", ]
  df$loc_yr <- droplevels(df$loc_yr)
  df$pop    <- droplevels(df$pop)
}

# geno_with_pheno NACH bibon-Drop — damit die Genotyp-Namen mit blup_df übereinstimmen
geno_with_pheno <- df %>%
  filter(!is.na(.data[[TRAIT]]), !is.na(Genotype)) %>%
  pull(Genotype) %>% as.character() %>% unique()

# ── Ginv ──────────────────────────────────────────────────────────────────────
Ginv <- readRDS(file.path(SNP_DIR, "MASTER_Ginv_flt_imp", "Ginv_FQ.rds"))
missing_from_ginv <- setdiff(levels(df$Genotype), attr(Ginv, "rowNames"))
df <- df %>%
  mutate(Genotype = if_else(Genotype %in% missing_from_ginv,
                            factor(NA, levels = levels(Genotype)),
                            Genotype))

# ── Sort ──────────────────────────────────────────────────────────────────────
if (DATASET == "JA") {
  df <- df[order(df$loc_yr, df$row, df$column), ]
} else {
  df <- df[order(df$location, df$year, df$row, df$column, df$rep), ]
}

# ── Fit model ─────────────────────────────────────────────────────────────────
fixed_frm <- as.formula(paste(TRAIT, "~ pop"))

if (DATASET == "JA") {
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = ~ corh(loc_yr):vm(Genotype, Ginv)
    + at(loc_yr):row + at(loc_yr):column,
    residual  = ~ dsum(~ ar1(row):ar1(column) | loc_yr),
    data      = df,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE, maxit = 30, workspace = "2gb"
  )
} else {
  m_fit <- asreml(
    fixed     = fixed_frm,
    random    = ~ corh(loc_yr):vm(Genotype, Ginv)
    + at(loc_yr):row + at(loc_yr):column,
    residual  = ~ dsum(~ idv(units) | loc_yr),
    data      = df,
    na.action = na.method(y = "include", x = "include"),
    ai.sing   = FALSE, maxit = 30, workspace = "2gb"
  )
}
while (!m_fit$converge) m_fit <- update(m_fit)

# ── Variance components ───────────────────────────────────────────────────────
vc <- summary(m_fit)$varcomp
rn <- rownames(vc)

vg_idx  <- grep("loc_yr:vm\\(Genotype, Ginv\\)!loc_yr_", rn)
vg_idx  <- vg_idx[!grepl("!loc_yr!cor", rn[vg_idx])]
vg_rows <- vc[vg_idx, , drop = FALSE]
vg_rows <- vg_rows[vg_rows$bound != "B", , drop = FALSE]
v_g     <- mean(vg_rows$component)

resid_pat  <- if (DATASET == "JA") "!R$" else "!units$"
resid_rows <- vc[grep(resid_pat, rn), , drop = FALSE]
resid_rows <- resid_rows[resid_rows$bound != "B", , drop = FALSE]
v_e <- mean(resid_rows$component)
h2  <- v_g / (v_g + v_e)

cat(sprintf("Trait: %s | h²=%.3f | v_g=%.4f | v_e=%.4f | LogLik=%.2f\n",
            TRAIT, h2, v_g, v_e, summary(m_fit)$loglik))

# ── Population BLUEs ─────────────────────────────────────────────────────────
pop_blues <- predict(m_fit, classify = "pop")$pvals %>%
  transmute(pop = as.character(pop), pop_BLUE = predicted.value)

# ── Extract BLUPs from random effect solutions ────────────────────────────────
blup_raw <- summary(m_fit, coef = TRUE)$coef.random
geno_idx <- grep("vm\\(Genotype, Ginv\\)", rownames(blup_raw))
blup_df  <- data.frame(
  Genotype = sub(".*_", "", rownames(blup_raw)[geno_idx]),
  GEBV     = blup_raw[geno_idx, "solution"],
  GEBV_se  = blup_raw[geno_idx, "std.error"],
  stringsAsFactors = FALSE
)

gebv <- blup_df %>%
  mutate(
    PEV         = GEBV_se^2,
    reliability = 1 - PEV / v_g
  ) %>%
  # NA setzen für Genotypen ohne Phänotypdaten
  mutate(across(c(GEBV, GEBV_se, PEV, reliability),
                ~ if_else(Genotype %in% geno_with_pheno, .x, NA_real_))) %>%
  left_join(geno_pop,  by = "Genotype") %>%
  left_join(pop_blues, by = "pop") %>%
  mutate(
    dGEBV  = GEBV / reliability,
    weight = (1 - h2 / reliability) / (1 - reliability)^2
  )

# ── Save ──────────────────────────────────────────────────────────────────────
out_file <- file.path(OUT_DIR, sprintf("gebv_%s.csv", TRAIT))
write.csv(gebv, out_file, row.names = FALSE)
cat(sprintf("Saved: %s\n", out_file))