library(tidyverse)

# ── Load data ────────────────────────────────────────────────────────────────
ja   <- read_delim("C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality/juice_antho_traits_OutRmv.csv",
                   delim = ";", na = "NA")
frgb <- read_delim("C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality/firmness_rgb_OutRmv.csv",
                   delim = ";", na = "NA")

# ── Trait columns to keep ────────────────────────────────────────────────────
ja_traits   <- c("anthocyanin", "brix", "acids", "a", "phenols")
frgb_traits <- c("firmness", "area_cm2", "h_w_ratio", "CIE_a")

# ── Count replicates per Genotype × pop × year × trait ───────────────────────
count_reps <- function(df, traits) {
  df %>%
    filter(!is.na(Genotype)) %>%
    pivot_longer(all_of(traits), names_to = "trait", values_to = "value") %>%
    group_by(Genotype, pop, year, trait) %>%
    summarise(reps = sum(!is.na(value)), .groups = "drop")
}

reps_all <- bind_rows(
  count_reps(ja,   ja_traits),
  count_reps(frgb, frgb_traits)
) %>%
  mutate(
    reps  = factor(reps, levels = as.character(0:8)),
    trait = factor(trait, levels = c(ja_traits, frgb_traits)),
    year  = factor(year)
  )

# ── Summarise: how many genotypes per pop have each rep-count? ───────────────
plot_df <- reps_all %>%
  count(pop, trait, year, reps, name = "n_geno")

# ── Colour palette ────────────────────────────────────────────────────────────
teal_pal <- c(
  "0" = "#C8BEB2",
  "1" = "#AADDD0",
  "2" = "#7ECEC0",
  "3" = "#52BFB0",
  "4" = "#36B0A0",
  "5" = "#269080",
  "6" = "#1B7060",
  "7" = "#105040",
  "8" = "#053020"
)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(plot_df, aes(x = pop, y = n_geno, fill = reps)) +
  geom_col(width = 0.7) +
  facet_grid(year ~ trait) +
  scale_fill_manual(values = teal_pal, name = "Reps", drop = FALSE) +
  scale_y_continuous(n.breaks = 8) +
  labs(
    title = "Replicate count by population, trait and year",
    x     = "Population",
    y     = "Number of genotypes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold"),
    panel.grid.minor  = element_blank(),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 11),
    legend.position   = "right",
    plot.title        = element_text(hjust = 0.5, face = "bold")
  )