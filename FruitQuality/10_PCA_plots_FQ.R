rm(list = ls()); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(car)
  library(tidyverse)
  library(svglite)
})

FQ_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
IN_DIR <- file.path(FQ_DIR, "GEBV_results")

TRAIT_ORDER  <- c("brix", "acids", "phenols", "firmness", "area_cm2", "CIE_a", "h_w_ratio")
TRAIT_LABELS <- c("Brix", "Acids", "Phenols", "Firmness", "Area", "CIE_a", "H/W")

pop_colors <- c(bibon = "#1f78b4", pxw = "#ff7f00", pxa = "#e31a1c", gxr = "#33a02c")
co_shapes  <- c("Wildtype" = 1,  "Mutant" = 16)
myb_shapes <- c("0 (wt)" = 1, "1 (het)" = 16, "2 (hom)" = 17)

# ── Read GEBV files ───────────────────────────────────────────────────────────
files <- list.files(IN_DIR, pattern = "^gebv_.*\\.csv$", full.names = TRUE)

gebv_raw <- purrr::map_dfr(files, function(f) {
  d <- read.csv(f, check.names = FALSE)
  d$trait <- sub("^gebv_(.+)\\.csv$", "\\1", basename(f))
  d
}) %>%
  transmute(
    Genotype = as.character(Genotype),
    trait    = as.character(trait),
    value    = as.numeric(dGEBV) + as.numeric(pop_BLUE)
  ) %>%
  mutate(trait = factor(trait, levels = TRAIT_ORDER, labels = TRAIT_LABELS)) %>%
  filter(!is.na(trait), is.finite(value))

cat(sprintf("Traits loaded: %s\n", paste(levels(gebv_raw$trait), collapse = ", ")))

# ── Metadata ──────────────────────────────────────────────────────────────────
meta_raw <- bind_rows(
  read.csv(file.path(FQ_DIR, "firmness_rgb_OutRmv.csv"),       sep = ";"),
  read.csv(file.path(FQ_DIR, "juice_antho_traits_OutRmv.csv"), sep = ";")
) %>%
  filter(!is.na(Genotype)) %>%
  transmute(
    Genotype = as.character(Genotype),
    pop      = as.character(pop),
    co_lab   = case_when(co_loc == 1 ~ "Mutant",
                         co_loc == 0 ~ "Wildtype",
                         TRUE        ~ NA_character_) %>%
      factor(levels = names(co_shapes)),
    myb_lab  = case_when(myb == 0 ~ "0 (wt)",
                         myb == 1 ~ "1 (het)",
                         myb == 2 ~ "2 (hom)",
                         TRUE     ~ NA_character_) %>%
      factor(levels = names(myb_shapes))
  ) %>%
  distinct(Genotype, .keep_all = TRUE)

# ── Wide matrix: z-standardise within trait, impute mean for missing ──────────
wide <- gebv_raw %>%
  group_by(Genotype, trait) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(trait) %>%
  mutate(z = scale(value)[, 1]) %>%
  ungroup() %>%
  select(Genotype, trait, z) %>%
  pivot_wider(names_from = trait, values_from = z)

cat(sprintf("Genotypes before imputation: %d | with complete data: %d\n",
            nrow(wide), sum(complete.cases(wide))))

X <- wide %>%
  column_to_rownames("Genotype") %>%
  mutate(across(everything(),
                ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

cat(sprintf("Genotypes in PCA: %d\n", nrow(X)))

meta <- wide %>%
  select(Genotype) %>%
  left_join(meta_raw, by = "Genotype") %>%
  mutate(pop = factor(pop, levels = names(pop_colors)))

# ── PCA ───────────────────────────────────────────────────────────────────────
pr      <- prcomp(X, center = TRUE, scale. = FALSE)
var_exp <- 100 * pr$sdev^2 / sum(pr$sdev^2)

scores <- as.data.frame(pr$x[, 1:2]) %>%
  rownames_to_column("Genotype") %>%
  left_join(meta, by = "Genotype")

# remove extreme outliers (>4 SD on either PC) before plotting and arrow scaling
scores_clean <- scores %>%
  filter(abs(PC1) < 4 * sd(PC1) & abs(PC2) < 4 * sd(PC2))

cat(sprintf("Genotypes removed as outliers (>4 SD): %d\n",
            nrow(scores) - nrow(scores_clean)))

# scale arrows to score cloud of cleaned data
loadings    <- as.data.frame(pr$rotation[, 1:2])
scale_arrow <- 0.4 * mean(c(diff(range(scores_clean$PC1)),
                            diff(range(scores_clean$PC2))))
loadings_sc <- loadings * scale_arrow
loadings_sc$Trait <- rownames(loadings_sc)

arrow_layers <- list(
  geom_segment(data = loadings_sc,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.18, "cm")),
               colour = "black", linewidth = 0.6),
  geom_text(data = loadings_sc,
            aes(PC1 * 1.15, PC2 * 1.15, label = Trait),
            inherit.aes = FALSE,
            fontface = "bold", size = 3.8),
  coord_fixed(),
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])),
  theme_minimal(base_size = 14),
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
)

# ── Plot 1: colour = pop, shape = MYB10 ──────────────────────────────────────
p1 <- scores_clean %>%
  filter(!is.na(myb_lab)) %>%
  ggplot(aes(PC1, PC2, colour = pop, shape = myb_lab)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = pop_colors, name = "Population") +
  scale_shape_manual(values  = myb_shapes, name = "MYB10") +
  arrow_layers +
  ggtitle("PCA — population (colour) × MYB10 (shape)")

print(p1)
svglite(file.path(FQ_DIR, "pca_pop_myb_FQ.svg"), width = 9, height = 7)
print(p1); dev.off()
cat("Saved: pca_pop_myb_FQ.svg\n")

# ── Plot 2: colour = pop, shape = Co locus ────────────────────────────────────
p2 <- scores_clean %>%
  filter(!is.na(co_lab)) %>%
  ggplot(aes(PC1, PC2, colour = pop, shape = co_lab)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_colour_manual(values = pop_colors, name = "Population") +
  scale_shape_manual(values  = co_shapes,  name = "Co locus") +
  arrow_layers +
  ggtitle("PCA — population (colour) × Co locus (shape)")

print(p2)
svglite(file.path(FQ_DIR, "pca_pop_co_FQ.svg"), width = 9, height = 7)
print(p2); dev.off()
cat("Saved: pca_pop_co_FQ.svg\n")

# ── Variance partitioning ─────────────────────────────────────────────────────
vpart <- purrr::map_dfr(c("PC1", "PC2"), function(pc) {
  sub <- scores_clean %>% filter(!is.na(pop), !is.na(co_lab))
  mod <- lm(reformulate(c("co_lab", "pop"), response = pc), data = sub)
  a   <- Anova(mod, type = 2)
  ss_total <- sum((sub[[pc]] - mean(sub[[pc]]))^2)
  tibble(PC = pc, Term = rownames(a), Semi_R2 = a$`Sum Sq` / ss_total)
})
cat("\nVariance partitioning of PC1 and PC2:\n")
print(vpart)
write.csv(vpart, file.path(FQ_DIR, "pc_variance_partitioning_FQ.csv"), row.names = FALSE)