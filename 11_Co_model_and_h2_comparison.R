# ──────────────────────────────────────────────────────────────────────────────
# Part A — pop + Co fixed Cullis h2 (toggle TRAIT/YEAR manually)
#          writes: .../h2_comparison/h2_<trait>_<year>_popCo.csv
# ──────────────────────────────────────────────────────────────────────────────

rm(list=ls(all=TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors=FALSE, scipen=999)

library(dplyr)
library(readr)
library(ggplot2)
library(forcats)
library(scales)
library(svglite)   
library(asreml)
library(data.table)
library(snpReady)

# ── SWITCHES ────────────────────────────────────────────────────────────────
TRAIT <- "in_len"  # in_len, lsh_ang_1, tr_brz_prop, tr_coni, tr_len
YEAR  <- 2017      # 2017 , 2018

# ── PATHS ────────────────────────────────────────────────────────────────────
BASE_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory"
STW_DIR  <- file.path(BASE_DIR, "Lidar", "Single_trait_walkthrough")
SNP_DIR  <- file.path(BASE_DIR, "SNP Chip")

pheno_file <- file.path(STW_DIR, "outlier_removed_raw", sprintf("clean_%s_%d_pop.csv", TRAIT, YEAR))
GT_path    <- file.path(SNP_DIR, "GT_filtered_numeric_transposed.csv")
GM_path    <- file.path(SNP_DIR, "GM_v11.csv")

stopifnot(file.exists(pheno_file), file.exists(GT_path), file.exists(GM_path))

OUT_DIR  <- file.path(STW_DIR, "h2_comparison")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
out_file <- file.path(OUT_DIR, sprintf("h2_%s_%d_popCo.csv", TRAIT, YEAR))

# ── PHENO ────────────────────────────────────────────────────────────────────
df <- read.csv(pheno_file, check.names=FALSE)
stopifnot(all(c("Genotype","pop","co_loc","Loc","Row","Column","Plot_RC_ID","rep",TRAIT) %in% names(df)))

# ── GT ───────────────────────────────────────────────────────────────────────
GT <- fread(GT_path, data.table=FALSE)
colnames(GT)[1] <- "Genotype"
M <- as.matrix(GT[, -1, drop=FALSE])
rownames(M) <- GT$Genotype
colnames(M) <- gsub("^AX\\.", "AX-", colnames(M))

ids <- intersect(unique(df$Genotype[!is.na(df$Genotype)]), rownames(M))
M <- M[ids, , drop=FALSE]

# MAF ≥ 0.05
p_alt <- colMeans(M, na.rm=TRUE) / 2
maf   <- pmin(p_alt, 1 - p_alt)
M     <- M[, maf >= 0.05, drop=FALSE]

# Mean-impute
na_cols <- which(colSums(is.na(M)) > 0)
if (length(na_cols)) for (j in na_cols) M[is.na(M[, j]), j] <- mean(M[, j], na.rm=TRUE)

# ── DROP CO WINDOW (Chr10: 20–36 Mb) ─────────────────────────────────────────
GM <- fread(GM_path, data.table=FALSE)
co_snps <- unique(GM$name[GM$chrom == 10L & GM$pos >= 2e7 & GM$pos <= 3.6e7])
if (length(co_snps)) M <- M[, !colnames(M) %in% co_snps, drop=FALSE]

# ── GRM + INV ────────────────────────────────────────────────────────────────
A <- snpReady::G.matrix(M=M, method="VanRaden", format="wide")$Ga
diag(A) <- diag(A) + 1e-2
Ginv <- ASRgenomics::G.inverse(A, sparseform = FALSE)$Ginv
glev <- rownames(Ginv)

df$Genotype   <- factor(df$Genotype, levels=glev)
df$pop        <- factor(df$pop)
df$co_loc <- factor(df$co_loc)
df$Loc   <- factor(df$Loc)
df$Row        <- factor(df$Row)
df$Column     <- factor(df$Column)
df$Plot_RC_ID <- factor(df$Plot_RC_ID)
df$rep        <- factor(df$rep)

# ── FIT ──────────────────────────────────────────────────────────────────────
m_fit <- asreml(
  fixed     = as.formula(paste0(TRAIT, " ~ pop + co_loc")),
  random    = ~ vm(Genotype, Ginv) + at(Loc):Row + at(Loc):Column,
  residual  = ~ dsum(~ id(Plot_RC_ID):ar1(rep) | Loc),
  data      = df,
  na.action = na.method(y="include", x="include"),
  ai.sing   = FALSE
)
m_fit <- update(m_fit, aom=TRUE)

# ── CULLIS h2 ────────────────────────────────────────────────────────────────
vc <- summary(m_fit)$varcomp
vc_df <- as.data.frame(vc)
rn <- rownames(vc_df)
v  <- vc_df$component; names(v) <- rn

idx_vm <- grepl("^vm\\(Genotype,\\s*Ginv\\)", rn)
stopifnot(any(idx_vm))
Gvar <- sum(v[idx_vm], na.rm=TRUE)

predG <- predict(m_fit, classify="Genotype", sed=TRUE)
AVSED <- as.numeric(predG$avsed["mean"])
h2_Cullis <- 1 - (AVSED^2)/(2*Gvar)
write.csv(data.frame(Trait=TRAIT, Year=YEAR, Method="pop+Co fixed", Gvar=Gvar, AVSED=AVSED, h2_Cullis=h2_Cullis),
          out_file, row.names=FALSE)

# ──────────────────────────────────────────────────────────────────────────────
# Part B — Load ALL trait×year h2 (pop fixed & pop+Co fixed) and plot one figure
# ──────────────────────────────────────────────────────────────────────────────
rm(list=ls(all=TRUE)); graphics.off(); closeAllConnections()

DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"
TRAITS   <- c("in_len","lsh_ang_1","tr_brz_prop","tr_coni","tr_len")
YEARS    <- c(2017, 2018)
# =====================

dir_stw     <- file.path(DATA_DIR, "Single_trait_walkthrough")
dir_popfix  <- file.path(dir_stw, "h2_and_deregression_input") # pop-fixed h2 source
dir_compare <- file.path(dir_stw, "h2_comparison")             # pop+Co h2 source/plot target
dir.create(dir_compare, showWarnings = FALSE, recursive = TRUE)

read_pop_fixed <- function(tr, yr){
  f <- file.path(dir_popfix, sprintf("BLUPs_SEs_h2_%s_%d_pop.csv", tr, yr))
  stopifnot(file.exists(f))
  tab <- read.csv(f, check.names = FALSE)
  stopifnot(all(c("AVSED","h2_Cullis") %in% names(tab)))
  data.frame(
    Trait    = tr,
    Year     = yr,
    Method   = "pop fixed",
    Gvar     = NA_real_,
    AVSED    = tab$AVSED[1],
    h2_Cullis = tab$h2_Cullis[1],
    check.names = FALSE
  )
}

read_popCo <- function(tr, yr){
  f <- file.path(dir_compare, sprintf("h2_%s_%d_popCo.csv", tr, yr))
  stopifnot(file.exists(f))
  read.csv(f, check.names = FALSE) |>
    mutate(Method = "pop+Co fixed") |>
    select(Trait, Year, Method, Gvar, AVSED, h2_Cullis)
}

pop_fixed <- do.call(
  rbind,
  lapply(TRAITS, function(tr)
    do.call(rbind, lapply(YEARS, function(yr) read_pop_fixed(tr, yr))))
)

popCo <- do.call(
  rbind,
  lapply(TRAITS, function(tr)
    do.call(rbind, lapply(YEARS, function(yr) read_popCo(tr, yr))))
)
plot_df <- bind_rows(pop_fixed, popCo) |>
  mutate(
    Trait = recode(Trait,
                   "in_len"      = "INL",
                   "lsh_ang_1"   = "LSA",
                   "tr_brz_prop" = "TBZP",
                   "tr_coni"     = "TC",
                   "tr_len"      = "TL"),
    Trait  = factor(Trait,
                    levels = c("INL","LSA","TBZP","TC","TL")),
    Year   = factor(Year, levels = sort(unique(Year))),
    Method = factor(Method, levels = c("pop fixed","pop+Co fixed"))
  )

method_cols <- c("pop fixed" = "#1f77b4", "pop+Co fixed" = "#ffa500")  

p <- ggplot(plot_df, aes(Trait, h2_Cullis, fill = Method)) +
  geom_col(
    position = position_dodge(width = 0.58),
    width    = 0.52
  ) +   # no border color → cleaner
  facet_grid(Year ~ ., switch = "y") +
  scale_fill_manual(values = method_cols, name = "Model") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "h\u00B2"   # keep y-label for clarity
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),
    panel.spacing      = unit(0.8, "lines"),
    
    axis.text.x        = element_text(angle = 30, hjust = 1, vjust = 1, size = 11),
    axis.title.y       = element_text(margin = margin(r = 8)),
    
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "white", colour = "grey60", linewidth = 0.4),
    strip.text.y.left  = element_text(angle = 0, face = "bold", size = 12),
    
    legend.position    = "right",
    legend.direction   = "vertical",
    legend.title       = element_text(size = 13),
    legend.text        = element_text(size = 11),
    legend.key.width   = unit(14, "pt"),
    legend.key.height  = unit(10, "pt"),
    
    plot.title         = element_blank(),  # no title in the figure
    plot.margin        = margin(t = 6, r = 14, b = 8, l = 8)
  )

print(p)

# Save as SVG (vector, editable)
svg_path <- file.path(dir_compare, "h2_pop_vs_popCo.svg")
ggsave(svg_path, plot = p, width = 12, height = 8, device = svglite)
cat("✔ Saved SVG to: ", svg_path, "\n", sep = "")
