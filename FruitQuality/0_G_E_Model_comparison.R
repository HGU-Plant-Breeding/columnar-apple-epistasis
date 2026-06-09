rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(asreml)
  vm <- get("asr_vm", envir = asNamespace("asreml"))
  library(dplyr)
  library(data.table)
})

FQ_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
SNP_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip"
GRM_DIR <- file.path(SNP_DIR, "MASTER_Ginv_flt_imp")
OUT_DIR <- file.path(FQ_DIR, "GxE_model_comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Toggles ───────────────────────────────────────────────────────────────────
DATASET <- "JA"   # "JA" or "FRGB"
TRAIT   <- "brix"

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

cat(sprintf("Rows: %d | Non-NA %s: %d\n", nrow(df), TRAIT,
            sum(!is.na(df[[TRAIT]]))))

# anthocyanin restriction
if (TRAIT == "anthocyanin") {
  df$pop    <- factor(df$pop, levels = c("bibon", "pxw"))
  df        <- df[df$location != "fuchsberg", ]
  df$loc_yr <- droplevels(df$loc_yr)
}

# ── Load pre-built Ginv ───────────────────────────────────────────────────────
cat("Loading Ginv_FQ.rds...\n")
Ginv <- readRDS(file.path(GRM_DIR, "Ginv_FQ.rds"))

missing_from_ginv <- setdiff(levels(df$Genotype), attr(Ginv, "rowNames"))
df <- df %>%
  mutate(Genotype = if_else(Genotype %in% missing_from_ginv,
                            factor(NA, levels = levels(Genotype)),
                            Genotype))

# ── Drop bibon for JA GxE models (too sparse: 7 and 13 obs) ─────────────────
if (DATASET == "JA") {
  cat(sprintf("Dropping bibon — obs before: %d\n", sum(!is.na(df[[TRAIT]]))))
  df        <- df[df$location != "bibon", ]
  df$loc_yr <- droplevels(df$loc_yr)
  df$pop    <- droplevels(df$pop)
  cat(sprintf("Obs after bibon removal: %d\n", sum(!is.na(df[[TRAIT]]))))
}
# ── Sort ──────────────────────────────────────────────────────────────────────
if (DATASET == "JA") {
  df_fit <- df[order(df$loc_yr, df$row, df$column), ]
} else {
  df_fit <- df[order(df$location, df$year, df$row, df$column, df$rep), ]
}

fixed_frm <- as.formula(paste(TRAIT, "~ pop + rootstock"))

# ── VC extraction function ────────────────────────────────────────────────────
# Handles diag, corh, us, fa — all produce different rowname patterns
extract_vc <- function(m, dataset, gxe_model) {
  
  if (!m$converge) warning("Model did not converge")
  s  <- summary(m)
  vc <- s$varcomp
  rn <- rownames(vc)
  
  # ── Fit metrics ─────────────────────────────────────────────────────────────
  loglik <- s$loglik
  aic    <- s$aic
  bic    <- s$bic
  n_par  <- attr(aic, "parameters")
  
  # ── Residual variance per environment ────────────────────────────────────────
  resid_pat  <- if (dataset == "JA") "!R$" else "!units$"
  resid_rows <- vc[grep(resid_pat, rn), , drop = FALSE]
  resid_rows <- resid_rows[resid_rows$bound != "B", , drop = FALSE]
  v_e_mean   <- mean(resid_rows$component)
  v_e_sd     <- sd(resid_rows$component)
  
  # ── Genetic variance per environment ─────────────────────────────────────────
  # diag / corh: pattern "loc_yr:vm(Genotype, Ginv)!loc_yr_ENV"
  # us:          pattern "loc_yr:vm(Genotype, Ginv)!loc_yr_ENV:ENV" (diagonal only)
  # fa:          pattern "fa(loc_yr):vm(Genotype, Ginv)!ENV!var" + "!ENVfa1"
  
  if (gxe_model %in% c("diag", "corh")) {
    vg_idx  <- grep("loc_yr:vm\\(Genotype, Ginv\\)!loc_yr_", rn)
    # exclude the correlation row from corh (!loc_yr!cor)
    vg_idx  <- vg_idx[!grepl("!loc_yr!cor", rn[vg_idx])]
    vg_rows <- vc[vg_idx, , drop = FALSE]
    vg_rows <- vg_rows[vg_rows$bound != "B", , drop = FALSE]
    env_names <- sub(".*!loc_yr_", "", rownames(vg_rows))
    v_g_vec   <- setNames(vg_rows$component, env_names)
    
    # corh: also extract the common correlation
    cor_val <- if (gxe_model == "corh") {
      cr <- vc[grep("!loc_yr!cor", rn), , drop = FALSE]
      if (nrow(cr)) cr$component[1] else NA_real_
    } else NA_real_
    
  } else if (gxe_model == "us") {
    # diagonal elements only: ENV:ENV
    vg_idx  <- grep("loc_yr:vm\\(Genotype, Ginv\\)!loc_yr_", rn)
    # keep only diagonal (colname pattern ends with :SAME_ENV)
    diag_idx <- vg_idx[sapply(rn[vg_idx], function(x) {
      parts <- strsplit(sub(".*!loc_yr_", "", x), ":")[[1]]
      length(parts) == 2 && parts[1] == parts[2]
    })]
    vg_rows   <- vc[diag_idx, , drop = FALSE]
    vg_rows   <- vg_rows[vg_rows$bound != "B", , drop = FALSE]
    env_names <- sub(":.*", "", sub(".*!loc_yr_", "", rownames(vg_rows)))
    v_g_vec   <- setNames(vg_rows$component, env_names)
    cor_val   <- NA_real_
    
  } else if (gxe_model == "fa") {
    # specific variance (!var) + factor loading (!fa1) per environment
    var_idx   <- grep("!var$", rn)
    var_rows  <- vc[var_idx, , drop = FALSE]
    fa1_idx   <- grep("!fa1$", rn)
    fa1_rows  <- vc[fa1_idx, , drop = FALSE]
    env_names <- sub(".*\\)!", "", sub("!var$", "", rownames(var_rows)))
    # total genetic variance per env = specific var + loading²
    # (under FA1 model: σ²_g_env = ψ_env + λ_env²)
    v_g_vec   <- setNames(
      var_rows$component + fa1_rows$component^2,
      env_names
    )
    cor_val   <- NA_real_
  }
  
  # ── Per-environment h² ───────────────────────────────────────────────────────
  # align v_e by environment name
  env_order  <- names(v_g_vec)
  v_e_named  <- setNames(resid_rows$component,
                         sub(".*_([^_]+_[0-9]{4}).*", "\\1", rownames(resid_rows)))
  # match order
  v_e_vec    <- v_e_named[env_order]
  h2_per_env <- v_g_vec / (v_g_vec + v_e_vec)
  
  # ── Spatial variance summary ──────────────────────────────────────────────────
  sp_idx  <- grep(":(row|column)$", rn)
  sp_rows <- vc[sp_idx, , drop = FALSE]
  sp_rows <- sp_rows[sp_rows$bound != "B", , drop = FALSE]
  v_sp_mean <- if (nrow(sp_rows)) mean(sp_rows$component) else NA_real_
  
  # ── Assemble output rows ─────────────────────────────────────────────────────
  env_df <- data.frame(
    trait      = TRAIT,
    dataset    = dataset,
    gxe_model  = gxe_model,
    converged  = m$converge,
    loglik     = loglik,
    aic        = as.numeric(aic),
    bic        = as.numeric(bic),
    n_par      = n_par,
    common_cor = cor_val,
    v_e_mean   = v_e_mean,
    v_e_sd     = v_e_sd,
    v_sp_mean  = v_sp_mean,
    loc_yr     = names(h2_per_env),
    v_g        = as.numeric(v_g_vec),
    v_e        = as.numeric(v_e_vec),
    h2         = as.numeric(h2_per_env),
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("\n── %s | loglik=%.2f AIC=%.2f BIC=%.2f ──\n",
              gxe_model, loglik, as.numeric(aic), as.numeric(bic)))
  cat("h² per environment:\n")
  print(round(h2_per_env, 3))
  cat(sprintf("Mean h² = %.3f\n", mean(h2_per_env, na.rm = TRUE)))
  
  env_df
}

# ── Fit all four models ───────────────────────────────────────────────────────
results <- list()

## 1. diag ────────────────────────────────────────────────────────────────────
cat("\n====== Fitting: diag ======\n")
tryCatch({
  if (DATASET == "JA") {
    m <- asreml(fixed = fixed_frm,
                random   = ~ diag(loc_yr):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ ar1(row):ar1(column) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = FALSE, maxit = 30, workspace = "2gb")
  } else {
    m <- asreml(fixed = fixed_frm,
                random   = ~ diag(loc_yr):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ idv(units) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = FALSE, maxit = 30, workspace = "2gb")
  }
  while (!m$converge) m <- update(m)
  results[["diag"]] <- extract_vc(m, DATASET, "diag")
}, error = function(e) cat(sprintf("diag FAILED: %s\n", e$message)))

## 2. corh ────────────────────────────────────────────────────────────────────
cat("\n====== Fitting: corh ======\n")
tryCatch({
  if (DATASET == "JA") {
    m <- asreml(fixed = fixed_frm,
                random   = ~ corh(loc_yr):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ ar1(row):ar1(column) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = FALSE, maxit = 30, workspace = "2gb")
  } else {
    m <- asreml(fixed = fixed_frm,
                random   = ~ corh(loc_yr):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ idv(units) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = FALSE, maxit = 30, workspace = "2gb")
  }
  while (!m$converge) m <- update(m)
  results[["corh"]] <- extract_vc(m, DATASET, "corh")
}, error = function(e) cat(sprintf("corh FAILED: %s\n", e$message)))

## 3. us ──────────────────────────────────────────────────────────────────────
cat("\n====== Fitting: us ======\n")
tryCatch({
  if (DATASET == "JA") {
    m <- asreml(fixed = fixed_frm,
                random   = ~ us(loc_yr):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ ar1(row):ar1(column) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = TRUE, maxit = 30, workspace = "2gb")
  } else {
    m <- asreml(fixed = fixed_frm,
                random   = ~ us(loc_yr):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ idv(units) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = TRUE, maxit = 30, workspace = "2gb")
  }
  while (!m$converge) m <- update(m)
  results[["us"]] <- extract_vc(m, DATASET, "us")
}, error = function(e) cat(sprintf("us FAILED: %s\n", e$message)))

## 4. fa(1) ───────────────────────────────────────────────────────────────────
cat("\n====== Fitting: fa(1) ======\n")
tryCatch({
  if (DATASET == "JA") {
    m <- asreml(fixed = fixed_frm,
                random   = ~ fa(loc_yr, 1):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ ar1(row):ar1(column) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = FALSE, maxit = 30, workspace = "2gb")
  } else {
    m <- asreml(fixed = fixed_frm,
                random   = ~ fa(loc_yr, 1):vm(Genotype, Ginv)
                + at(loc_yr):row + at(loc_yr):column,
                residual = ~ dsum(~ idv(units) | loc_yr),
                data = df_fit, na.action = na.method(y="include", x="include"),
                ai.sing = FALSE, maxit = 30, workspace = "2gb")
  }
  while (!m$converge) m <- update(m)
  results[["fa"]] <- extract_vc(m, DATASET, "fa")
}, error = function(e) cat(sprintf("fa FAILED: %s\n", e$message)))

# ── Combine and save ──────────────────────────────────────────────────────────
if (length(results) == 0) stop("All models failed.")

new_rows <- bind_rows(results)

out_csv <- file.path(OUT_DIR, sprintf("GxE_comparison_%s.csv", DATASET))

if (file.exists(out_csv)) {
  existing <- read.csv(out_csv, check.names = FALSE)
  existing <- existing[!(existing$trait == TRAIT & 
                           existing$gxe_model %in% names(results)), ]
  combined <- bind_rows(existing, new_rows)
} else {
  combined <- new_rows
}
write.csv(combined, out_csv, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_csv))

# ── Summary table ─────────────────────────────────────────────────────────────
summary_tbl <- combined %>%
  filter(trait == TRAIT) %>%
  group_by(gxe_model) %>%
  summarise(
    loglik     = first(loglik),
    AIC        = first(aic),
    BIC        = first(bic),
    n_par      = first(n_par),
    common_cor = first(common_cor),
    mean_h2    = mean(h2, na.rm = TRUE),
    mean_v_g   = mean(v_g, na.rm = TRUE),
    mean_v_e   = mean(v_e, na.rm = TRUE),
    converged  = first(converged),
    .groups    = "drop"
  ) %>%
  arrange(AIC)

cat(sprintf("\n── Model comparison: %s (%s) ──\n", TRAIT, DATASET))
print(summary_tbl, digits = 3)