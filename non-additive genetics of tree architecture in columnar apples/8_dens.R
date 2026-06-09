# ─────────────────────────────────────────────────────────────
#  Setup
# ─────────────────────────────────────────────────────────────
rm(list = ls()); graphics.off(); closeAllConnections()

library(tidyverse)
library(patchwork)   # for combining plots
library(grid)        # for unit(), margin()

# ─────────────────────────────────────────────────────────────
#  Paths
# ─────────────────────────────────────────────────────────────
BASE_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
IN_DIR   <- file.path(BASE_DIR, "deregression_output")
setwd(BASE_DIR)

# ─────────────────────────────────────────────────────────────
#  Palettes
# ─────────────────────────────────────────────────────────────
population_colors <- c(
  pxw = "#ff7f00",  # PxW
  pxa = "#e31a1c",  # P28xA14
  gxr = "#33a02c"   # GxR
)

mutation_colors <- c(
  Wildtype = "#1b9e77",
  Mutant   = "#d95f02"
)

# ─────────────────────────────────────────────────────────────
#  Read & prep
# ─────────────────────────────────────────────────────────────
files <- list.files(IN_DIR, pattern = "^dGEBV_.*_pop\\.csv$", full.names = TRUE)

df <- purrr::map_dfr(files, read.csv, check.names = FALSE) %>%
  transmute(
    Year       = as.integer(Year),
    Trait      = as.character(Trait),
    Genotype,
    pop        = na_if(trimws(as.character(pop)), ""),
    co_locus   = co_loc,
    dGEBV      = as.numeric(dGEBV),
    pop_BLUE   = as.numeric(pop_BLUE)
  ) %>%
  mutate(
    co_lab = case_when(
      co_locus %in% c(1, "1", "Mutant")   ~ "Mutant",
      co_locus %in% c(0, "0", "Wildtype") ~ "Wildtype",
      TRUE ~ NA_character_
    )
  )

trait_levels <- sort(unique(df$Trait))

# Z-standardized (dGEBV + pop_BLUE) within Year × Trait
df <- df %>%
  mutate(
    Trait = recode(Trait,
                   "in_len"       = "INL",
                   "lsh_ang_1"    = "LSA",
                   "tr_brz_prop"  = "TBZP",
                   "tr_coni"      = "TC",
                   "tr_len"       = "TL"),
    Trait = factor(Trait,
                   levels = c("INL","LSA","TBZP","TC","TL"))
  ) %>%
  group_by(Year, Trait) %>%
  mutate(
    sum_val = dGEBV + pop_BLUE,
    zsum    = (sum_val - mean(sum_val, na.rm = TRUE)) /
      sd(sum_val,  na.rm = TRUE)
  ) %>%
  ungroup()
# ─────────────────────────────────────────────────────────────
#  Function: make combined density plot for a given year
# ─────────────────────────────────────────────────────────────
make_density_plot_year <- function(year_val, data, out_file) {
  
  df_year <- data %>%
    filter(Year == year_val, is.finite(zsum))
  
  if (nrow(df_year) == 0L) {
    message("No data for year ", year_val, "; skipping.")
    return(invisible(NULL))
  }
  
  # Top row: by population
  df_pop <- df_year %>%
    filter(!is.na(pop)) %>%
    mutate(pop = factor(pop, levels = c("pxw","pxa","gxr")))
  
  p_pop <- ggplot(df_pop, aes(x = zsum, fill = pop)) +
    geom_density(alpha = 0.55) +
    facet_wrap(~ Trait, nrow = 1, scales = "free_x") +
    scale_fill_manual(
      values = population_colors,
      breaks = c("pxw","pxa","gxr"),
      labels = c("P×W",
                 "P28×A14",
                 "GxR"),
      name   = "Population"
    ) +
    scale_x_continuous(
      breaks = -3:3,
      limits = c(-3, 3),
      expand = expansion(mult = 0.05)
    ) +
    scale_y_continuous(breaks = NULL) +
    labs(
      x = NULL,
      y = "By population"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text.x     = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey90", colour = NA),
      axis.text.x      = element_text(size = 10),
      axis.ticks.y     = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      panel.spacing    = unit(1.2, "lines"),
      plot.margin      = margin(5, 5, 0, 5)
    )
  
  # Bottom row: by mutation (Co locus)
  df_co <- df_year %>%
    filter(!is.na(co_lab)) %>%
    mutate(co_lab = factor(co_lab, levels = c("Wildtype","Mutant")))
  
  p_co <- ggplot(df_co, aes(x = zsum, fill = co_lab)) +
    geom_density(alpha = 0.55) +
    facet_wrap(~ Trait, nrow = 1, scales = "free_x") +
    scale_fill_manual(
      values = mutation_colors,
      name   = "Co locus"
    ) +
    scale_x_continuous(
      breaks = -3:3,
      limits = c(-3, 3),
      expand = expansion(mult = 0.05)
    ) +
    scale_y_continuous(breaks = NULL) +
    labs(
      x = "standardized (dGEBV+population_BLUE)",
      y = "By mutation"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text.x     = element_blank(),         # show only on top row
      strip.background = element_blank(),
      axis.text.x      = element_text(size = 10),
      axis.ticks.y     = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      panel.spacing    = unit(1.2, "lines"),
      plot.margin      = margin(0, 5, 5, 5)
    )
  
  # Combine vertically with a general title
  combined <- (p_pop / p_co) + plot_annotation(
    title = paste0(""),
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face  = "bold",
        size  = 15
      )
    )
  )
  
  # Save as SVG (vector; journal-friendly)
  svg(filename = out_file, width = 9, height = 6)
  print(combined)
  dev.off()
  
  invisible(combined)
}

# ─────────────────────────────────────────────────────────────
#  Create one plot per year (2017 and 2018)
# ─────────────────────────────────────────────────────────────
years_to_plot <- sort(unique(df$Year))

for (yy in years_to_plot) {
  out_svg <- file.path(
    BASE_DIR,
    paste0("Density_", yy, "_combined.svg")
  )
  make_density_plot_year(yy, df, out_svg)
}
