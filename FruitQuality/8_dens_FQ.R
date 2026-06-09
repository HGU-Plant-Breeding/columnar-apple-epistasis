rm(list = ls()); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(svglite)
})

FQ_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"
IN_DIR <- file.path(FQ_DIR, "GEBV_results")

# ── Trait order and labels ─────────────────────────────────────────────
TRAIT_ORDER  <- c("brix", "acids", "phenols", "firmness", "area_cm2", "CIE_a", "h_w_ratio")
TRAIT_LABELS <- c("Brix", "Acids", "Phenols", "Firmness", "Area", "CIE_a", "H/W")

# ── Palettes ───────────────────────────────────────────────────────────
pop_colors <- c(
  bibon = "#1f78b4",
  pxw   = "#ff7f00",
  pxa   = "#e31a1c",
  gxr   = "#33a02c"
)

myb_colors <- c(
  "0 (wt)"  = "#1b9e77",
  "1 (het)" = "#d95f02",
  "2 (hom)" = "#7570b3"
)

co_colors <- c(
  "Wildtype" = "#1b9e77",
  "Mutant"   = "#d95f02"
)

# ── Read GEBV files + extract trait ────────────────────────────────────
files <- list.files(IN_DIR, pattern = "^gebv_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("No gebv_*.csv files found in ", IN_DIR)

gebv_raw <- map_dfr(files, function(f) {
  read.csv(f, check.names = FALSE) %>%
    mutate(trait = gsub("^gebv_|\\.csv$", "", basename(f)))
})

# ── Read metadata ──────────────────────────────────────────────────────
meta <- read.csv(file.path(FQ_DIR, "firmness_rgb_OutRmv.csv"), sep = ",") %>%
  filter(!is.na(Genotype)) %>%
  transmute(
    Genotype = as.character(Genotype),
    pop      = as.character(pop),
    myb      = as.character(myb),
    co_loc   = as.character(co_loc)
  ) %>%
  distinct()

# ── Combine and prepare ────────────────────────────────────────────────
df <- gebv_raw %>%
  transmute(
    Genotype = as.character(Genotype),
    trait    = as.character(trait),
    GEBV     = as.numeric(GEBV),
    pop_BLUE = as.numeric(pop_BLUE)
  ) %>%
  left_join(meta, by = "Genotype") %>%
  filter(!is.na(GEBV)) %>%
  mutate(
    trait = factor(trait, levels = TRAIT_ORDER, labels = TRAIT_LABELS),
    
    pop = factor(pop, levels = names(pop_colors)),
    
    myb_lab = case_when(
      myb == "0" ~ "0 (wt)",
      myb == "1" ~ "1 (het)",
      myb == "2" ~ "2 (hom)"
    ) %>% factor(levels = names(myb_colors)),
    
    co_lab = case_when(
      co_loc == "1" ~ "Mutant",
      co_loc == "0" ~ "Wildtype"
    ) %>% factor(levels = names(co_colors))
  ) %>%
  group_by(trait) %>%
  mutate(
    sum_val = GEBV + pop_BLUE,
    zsum    = (sum_val - mean(sum_val, na.rm = TRUE)) / sd(sum_val, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(is.finite(zsum), !is.na(trait))

# ── Theme + scale ──────────────────────────────────────────────────────
theme_dens <- theme_minimal(base_size = 13) +
  theme(
    strip.text.x     = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92", colour = NA),
    axis.ticks.y     = element_blank(),
    axis.text.y      = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.spacing    = unit(1.0, "lines")
  )

x_scale <- scale_x_continuous(
  breaks = -3:3,
  limits = c(-3.5, 3.5),
  expand = expansion(mult = 0.03)
)

# ── MYB10 ──────────────────────────────────────────────────────────────
p_myb <- df %>%
  filter(!is.na(myb_lab)) %>%
  ggplot(aes(x = zsum, fill = myb_lab)) +
  geom_density(alpha = 0.55, na.rm = TRUE) +
  geom_vline(
    data = df %>%
      filter(!is.na(myb_lab)) %>%
      group_by(trait, myb_lab) %>%
      summarise(m = mean(zsum), .groups = "drop"),
    aes(xintercept = m, color = myb_lab),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  facet_wrap(~ trait, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = myb_colors, name = "MYB10") +
  scale_color_manual(values = myb_colors) +
  x_scale +
  labs(x = NULL, y = "MYB10") +
  theme_dens

# ── Co locus ───────────────────────────────────────────────────────────
p_co <- df %>%
  filter(!is.na(co_lab)) %>%
  ggplot(aes(x = zsum, fill = co_lab)) +
  geom_density(alpha = 0.55, na.rm = TRUE) +
  geom_vline(
    data = df %>%
      filter(!is.na(co_lab)) %>%
      group_by(trait, co_lab) %>%
      summarise(m = mean(zsum), .groups = "drop"),
    aes(xintercept = m, color = co_lab),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  facet_wrap(~ trait, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = co_colors, name = "Co locus") +
  scale_color_manual(values = co_colors) +
  x_scale +
  labs(x = NULL, y = "Co locus") +
  theme_dens +
  theme(strip.text.x = element_blank(), strip.background = element_blank())

# ── Population ─────────────────────────────────────────────────────────
p_pop <- df %>%
  filter(!is.na(pop)) %>%
  ggplot(aes(x = zsum, fill = pop)) +
  geom_density(alpha = 0.55, na.rm = TRUE) +
  geom_vline(
    data = df %>%
      filter(!is.na(pop)) %>%
      group_by(trait, pop) %>%
      summarise(m = mean(zsum), .groups = "drop"),
    aes(xintercept = m, color = pop),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  facet_wrap(~ trait, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pop_colors, name = "Population") +
  scale_color_manual(values = pop_colors) +
  x_scale +
  labs(x = "standardised (dGEBV + population BLUE)", y = "Population") +
  theme_dens +
  theme(strip.text.x = element_blank(), strip.background = element_blank())

# ── Combine ────────────────────────────────────────────────────────────
combined <- (p_myb / p_co / p_pop) +
  plot_layout(heights = c(1, 1, 1))

print(combined)

# ── Save as SVG ────────────────────────────────────────────────────────
out_svg <- file.path(FQ_DIR, "dGEBV_density_FruitQuality.svg")

ggsave(
  filename = out_svg,
  plot     = combined,
  width    = 20,
  height   = 10,
  device   = svglite
)

cat("Saved:", out_svg, "\n")