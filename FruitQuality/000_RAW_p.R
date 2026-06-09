library(tidyverse)
library(ggplot2)
library(patchwork)
library(svglite)

# ── Data ──────────────────────────────────────────────────────────────────────
FRGB <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality/firmness_rgb_OutRmv.csv", sep = ";")
JA   <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality/juice_antho_traits_OutRmv.csv", sep = ";")

trait_cols   <- c("anthocyanin", "brix", "acids", "phenols",
                  "firmness", "area_cm2", "h_w_ratio", "CIE_a")
trait_labels <- c("Anthocyanin", "Brix", "Acids", "Phenols",
                  "Firmness", "Area (cm²)", "H/W ratio", "CIE a*")

# ── Genotyp-Mittelwerte pro Jahr ──────────────────────────────────────────────
ja_geno <- JA %>%
  filter(!is.na(Genotype)) %>%
  group_by(Genotype, year) %>%
  summarise(across(c(anthocyanin, brix, acids, phenols),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

frgb_geno <- FRGB %>%
  filter(!is.na(Genotype)) %>%
  group_by(Genotype, year) %>%
  summarise(across(c(firmness, area_cm2, h_w_ratio, CIE_a),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

all_geno <- full_join(ja_geno, frgb_geno, by = c("Genotype", "year"))

# ── Panel A: Heatmaps 2024 & 2025 ────────────────────────────────────────────
make_cor_df <- function(data, yr) {
  data %>%
    filter(year == yr) %>%
    select(all_of(trait_cols)) %>%
    cor(use = "pairwise.complete.obs", method = "pearson") %>%
    as.data.frame() %>%
    rownames_to_column("trait1") %>%
    pivot_longer(-trait1, names_to = "trait2", values_to = "r") %>%
    mutate(
      trait1 = factor(trait1, levels = trait_cols),
      trait2 = factor(trait2, levels = trait_cols),
      keep   = as.integer(trait1) > as.integer(trait2)
    ) %>%
    filter(keep) %>%
    mutate(
      label      = sprintf("%.2f", r),
      trait1_lbl = factor(trait_labels[as.integer(trait1)], levels = trait_labels),
      trait2_lbl = factor(trait_labels[as.integer(trait2)], levels = rev(trait_labels)),
      year       = as.character(yr)
    )
}

cor_df <- bind_rows(make_cor_df(all_geno, 2024),
                    make_cor_df(all_geno, 2025))

p_heat <- ggplot(cor_df, aes(x = trait1_lbl, y = trait2_lbl, fill = r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = abs(r) > 0.35), size = 2.9) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"),
                     guide = "none") +
  scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#993C1D",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson r") +
  scale_x_discrete(position = "bottom") +
  facet_wrap(~ year, ncol = 2) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y       = element_text(size = 9),
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(1.2, "cm"),
    strip.text        = element_text(size = 11),
    plot.tag          = element_text(size = 13, face = "bold")
  )

# ── Panel B: Cross-year dotplot ───────────────────────────────────────────────
all_wide <- full_join(
  ja_geno   %>% pivot_wider(names_from = year, values_from = c(anthocyanin, brix, acids, phenols),
                            names_glue = "{.value}_{year}"),
  frgb_geno %>% pivot_wider(names_from = year, values_from = c(firmness, area_cm2, h_w_ratio, CIE_a),
                            names_glue = "{.value}_{year}"),
  by = "Genotype"
)

cy_df <- lapply(seq_along(trait_cols), function(i) {
  tr <- trait_cols[i]
  c24 <- paste0(tr, "_2024"); c25 <- paste0(tr, "_2025")
  if (!c24 %in% names(all_wide) || !c25 %in% names(all_wide)) return(NULL)
  x <- all_wide[[c24]]; y <- all_wide[[c25]]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 10) return(NULL)
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  data.frame(trait = trait_labels[i], r = ct$estimate,
             ci_low = ct$conf.int[1], ci_high = ct$conf.int[2],
             n = sum(ok), stringsAsFactors = FALSE)
}) %>% bind_rows()

p_cy <- ggplot(cy_df, aes(x = r, y = reorder(trait, r))) +
  geom_vline(xintercept = c(0, 0.5, 1), linetype = c("dashed","dotted","dotted"),
             color = c("grey50","grey80","grey80"), linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.25, linewidth = 0.6, color = "grey50") +
  geom_point(aes(fill = r), shape = 21, size = 4, color = "white", stroke = 0.5) +
  geom_text(aes(label = sprintf("r=%.2f (n=%d)", r, n)),
            hjust = -0.15, size = 2.9, color = "grey30") +
  scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#993C1D",
                       midpoint = 0, limits = c(-1, 1), guide = "none") +
  scale_x_continuous(limits = c(-0.2, 1.3), breaks = seq(0, 1, 0.25)) +
  labs(x = "Pearson r (2024 vs. 2025)", y = NULL, tag = "B") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 10),
    plot.tag           = element_text(size = 13, face = "bold")
  )

# ── Kombinieren & speichern ───────────────────────────────────────────────────
p_combined <- p_heat / p_cy +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title    = "Phenotypic correlations",
    subtitle = "Genotype means · pairwise complete observations",
    theme    = theme(
      plot.title    = element_text(size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

out_path <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality/pheno_correlations.svg"
svglite(out_path, width = 11, height = 13)
print(p_combined)
dev.off()
cat(sprintf("Saved: %s\n", out_path))