# --------------------------------------------------------------
# Minimal diagnostics (one TRAIT × one YEAR) — RAW values, with outliers
# --------------------------------------------------------------
rm(list=ls(all=TRUE)); graphics.off(); closeAllConnections()

library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------
# Global switches
# --------------------
YEAR   <- 2017                 # <- 2017 or 2018
TRAIT  <- "in_len"        # <- one of: tr_len, in_len, tr_brz_prop, lsh_ang_1, tr_coni

# --------------------
# Paths & metadata
# --------------------
DATA_DIR  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
meta_cols <- c("Genotype","pop","Loc","Row","Column","rep","Plot_RC_ID","co_loc","Year")

# --------------------
# Load the selected year (with outliers)
# --------------------
infile <- if (YEAR == 2017) "df_17_final_w_outliers.csv" else "df_18_final_w_outliers.csv"

df <- read.csv(file.path(DATA_DIR, infile), sep = ";") %>%
  mutate(
    Genotype   = factor(Genotype),
    pop        = factor(pop),
    Loc   = factor(Loc),
    Row        = factor(Row),
    Column     = factor(Column),
    rep        = factor(rep),
    Plot_RC_ID = factor(Plot_RC_ID),
    co_loc = factor(co_loc),
    Year       = YEAR
  )

# Guardrails
if (!TRAIT %in% names(df))
  stop(paste0("Trait '", TRAIT, "' not found in file: ", infile))

# Keep just metadata + the chosen trait
df1 <- df %>%
  select(any_of(c(meta_cols, TRAIT)))

# --------------------
# Long format for one trait (RAW)
# --------------------
long1 <- df1 %>%
  pivot_longer(cols = all_of(TRAIT), names_to = "Trait", values_to = "y") %>%
  mutate(
    Year  = as.character(Year),
    Trait = factor(Trait, levels = TRAIT)  # single level
  )

title_suffix <- paste0(TRAIT, " × ", YEAR, " (raw, with outliers)")


# --------------------------------------------------------------
# 1) Replication balance (0/1/2/3 reps) for chosen TRAIT × YEAR
# --------------------------------------------------------------
rep_counts <- long1 %>%
  group_by(Year, Trait, Genotype) %>%
  summarise(n_reps = sum(!is.na(y)), .groups = "drop") %>%
  mutate(n_reps = pmin(n_reps, 3))

rep_share <- rep_counts %>%
  count(Year, Trait, n_reps, name = "n_geno") %>%
  group_by(Year, Trait) %>%
  mutate(share = n_geno / sum(n_geno)) %>%
  ungroup() %>%
  mutate(
    n_reps = factor(n_reps, levels = 0:3,
                    labels = c("0 reps","1 rep","2 reps","3 reps"))
  )

p_rep <- ggplot(rep_share, aes(x = n_reps, y = share, fill = n_reps)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, .02))) +
  scale_fill_brewer(palette = "Blues", direction = 1,
                    name = "Tree reps") +
  labs(title = paste0("Replication balance — ", title_suffix),
       x = "Observed reps per genotype", y = "Share of genotypes") +
  theme_bw(base_size = 13) +
  theme(
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold", hjust = .5),
    panel.grid.minor  = element_blank()
  )
print(p_rep)


# --------------------------------------------------------------
# 2) Histogram (RAW)
# --------------------------------------------------------------
p_hist <- ggplot(long1, aes(x = y)) +
  geom_histogram(bins = 30, alpha = 0.9, linewidth = 0.2) +
  geom_density(na.rm = TRUE, linewidth = 0.7) +
  labs(title = paste0("Histogram — ", title_suffix),
       x = TRAIT, y = "Count") +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = .5),
    panel.grid.minor = element_blank()
  )
print(p_hist)

# --------------------------------------------------------------
# 3) Boxplot by mutation (RAW)
# --------------------------------------------------------------
if (!all(is.na(long1$co_loc))) {
  p_mut <- long1 %>%
    filter(!is.na(co_loc)) %>%
    ggplot(aes(x = co_loc, y = y)) +
    geom_boxplot(outlier.shape = 16, outlier.size = .8, linewidth = 0.4) +
    labs(title = paste0("Mutation effect — ", title_suffix),
         x = "co_loc (mutation status)", y = TRAIT) +
    theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", hjust = .5),
      panel.grid.minor = element_blank()
    )
  print(p_mut)
} else {
  cat("Mutation plot skipped: 'co_loc' not available or all NA.\n")
}

# --------------------------------------------------------------
# 4) Boxplot by population (RAW)
# --------------------------------------------------------------
p_box <- long1 %>%
  filter(!is.na(pop)) %>%
  ggplot(aes(x = pop, y = y)) +
  geom_boxplot(outlier.shape = 16, outlier.size = .8, linewidth = 0.4) +
  labs(title = paste0("Population distributions — ", title_suffix),
       x = "Population", y = TRAIT) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = .5),
    panel.grid.minor = element_blank()
  )
print(p_box)

# --------------------------------------------------------------
# 5) Q–Q plot (RAW)
# --------------------------------------------------------------
qq_df <- long1 %>% filter(!is.na(y))
p_qq <- ggplot(qq_df, aes(sample = y)) +
  stat_qq(size = .6, alpha = .6) +
  stat_qq_line(linewidth = 0.6) +
  labs(title = paste0("Q–Q plot — ", title_suffix),
       x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = .5),
    panel.grid.minor = element_blank()
  )
print(p_qq)

