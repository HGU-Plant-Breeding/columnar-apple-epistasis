rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

FQ_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality"

DATASET <- "JA"      # "JA" or "FRGB"
TRAIT   <- "acids"

# JA traits:   anthocyanin, brix, acids, a, phenols
# FRGB traits: firmness, size_cm2, diameter_cm, h_w_ratio,
#              avg_R, avg_G, avg_B, circularity, CIE_a

# ── Load ──────────────────────────────────────────────────────────────────────
in_file <- file.path(FQ_DIR, if (DATASET == "JA") "juice_antho_traits_OutRmv.csv"
                     else                  "firmness_rgb_OutRmv.csv")

df <- read.csv(in_file) %>%
  mutate(
    Genotype = factor(Genotype),
    pop      = factor(pop),
    location = factor(location),
    year     = factor(year),
    myb      = factor(myb)
  )

if (!TRAIT %in% names(df))
  stop(paste0("Trait '", TRAIT, "' not found in dataset '", DATASET, "'"))

# convenience: year panels + "Both years" combined
df_plot <- bind_rows(
  df,
  df %>% mutate(year = factor("Both"))
) %>%
  mutate(year = factor(year, levels = c(sort(unique(as.character(df$year))), "Both")))

# ── 1) Histogram ──────────────────────────────────────────────────────────────
print(
  df_plot %>%
    filter(!is.na(.data[[TRAIT]])) %>%
    ggplot(aes(x = .data[[TRAIT]])) +
    geom_histogram(bins = 30, alpha = 0.85, linewidth = 0.2) +
    facet_wrap(~ year, scales = "free_y") +
    labs(title = paste0("Histogram — ", TRAIT, " (", DATASET, ")"),
         x = TRAIT, y = "Count") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = .5),
          panel.grid.minor = element_blank())
)

# ── 2) Q-Q plot ───────────────────────────────────────────────────────────────
print(
  df_plot %>%
    filter(!is.na(.data[[TRAIT]])) %>%
    ggplot(aes(sample = .data[[TRAIT]])) +
    stat_qq(size = .6, alpha = .6) +
    stat_qq_line(linewidth = 0.6) +
    facet_wrap(~ year, scales = "free") +
    labs(title = paste0("Q-Q plot — ", TRAIT, " (", DATASET, ")"),
         x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = .5),
          panel.grid.minor = element_blank())
)

# ── 3) Boxplot by population ───────────────────────────────────────────────────
print(
  df_plot %>%
    filter(!is.na(.data[[TRAIT]]), !is.na(pop)) %>%
    ggplot(aes(x = pop, y = .data[[TRAIT]], fill = pop)) +
    geom_boxplot(outlier.shape = 16, outlier.size = .8, linewidth = 0.4) +
    facet_wrap(~ year, scales = "free_y") +
    labs(title = paste0("Population distributions — ", TRAIT, " (", DATASET, ")"),
         x = "Population", y = TRAIT) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = .5),
          panel.grid.minor = element_blank(),
          legend.position = "none")
)

# ── 4) Boxplot by location ────────────────────────────────────────────────────
print(
  df_plot %>%
    filter(!is.na(.data[[TRAIT]])) %>%
    ggplot(aes(x = location, y = .data[[TRAIT]], fill = location)) +
    geom_boxplot(outlier.shape = 16, outlier.size = .8, linewidth = 0.4) +
    facet_wrap(~ year, scales = "free_y") +
    labs(title = paste0("Location distributions — ", TRAIT, " (", DATASET, ")"),
         x = "Location", y = TRAIT) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = .5),
          panel.grid.minor = element_blank(),
          legend.position = "none")
)

# ── 5) MYB10 effect (JA anthocyanin only) ────────────────────────────────────
if (DATASET == "JA" && TRAIT == "anthocyanin" && "myb" %in% names(df)) {
  print(
    df_plot %>%
      filter(!is.na(.data[[TRAIT]]), !is.na(myb)) %>%
      ggplot(aes(x = myb, y = .data[[TRAIT]], fill = myb)) +
      geom_boxplot(outlier.shape = 16, outlier.size = .8, linewidth = 0.4) +
      facet_wrap(~ year, scales = "free_y") +
      labs(title = paste0("MYB10 dosage effect — ", TRAIT, " (", DATASET, ")"),
           x = "MYB10 dosage (0/1/2)", y = TRAIT) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = .5),
            panel.grid.minor = element_blank(),
            legend.position = "none")
  )
}

# ── 6) Replication balance ────────────────────────────────────────────────────
# for FRGB: reps = fruit measurements per tree; for JA: 0 or 1 per tree×year
rep_df <- bind_rows(
  df %>% mutate(year_grp = as.character(year)),
  df %>% mutate(year_grp = "Both")
) %>%
  mutate(year_grp = factor(year_grp,
                           levels = c(sort(unique(as.character(df$year))), "Both"))) %>%
  group_by(year_grp, Genotype) %>%
  summarise(n_obs = sum(!is.na(.data[[TRAIT]])), .groups = "drop") %>%
  mutate(n_obs = pmin(n_obs, if (DATASET == "FRGB") 8 else 2)) %>%
  count(year_grp, n_obs, name = "n_geno") %>%
  group_by(year_grp) %>%
  mutate(share = n_geno / sum(n_geno)) %>%
  ungroup() %>%
  mutate(n_obs = factor(n_obs))

print(
  ggplot(rep_df, aes(x = n_obs, y = share, fill = n_obs)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, .02))) +
    facet_wrap(~ year_grp) +
    labs(title = paste0("Observation balance — ", TRAIT, " (", DATASET, ")"),
         x = "Observations per genotype", y = "Share of genotypes") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = .5),
          panel.grid.minor = element_blank(),
          legend.position = "none")
)