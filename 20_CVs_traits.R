rm(list = ls(all = TRUE)); graphics.off()
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)

# ── File paths ────────────────────────────────────────────────────────────────
DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"
dir_stw <- file.path(DATA_DIR, "Single_trait_walkthrough")
dir_raw <- file.path(dir_stw, "outlier_removed_raw")
dir.create(dir_stw, showWarnings = FALSE, recursive = TRUE)

files <- c(
  "clean_in_len_2017_pop.csv", "clean_in_len_2018_pop.csv",
  "clean_lsh_ang_1_2017_pop.csv","clean_lsh_ang_1_2018_pop.csv",
  "clean_tr_brz_1_2017_pop.csv","clean_tr_brz_1_2018_pop.csv",
  "clean_tr_brz_2_2017_pop.csv","clean_tr_brz_2_2018_pop.csv",
  "clean_tr_brz_prop_2017_pop.csv","clean_tr_brz_prop_2018_pop.csv",
  "clean_tr_coni_2017_pop.csv","clean_tr_coni_2018_pop.csv",
  "clean_tr_dia_1_2017_pop.csv","clean_tr_dia_1_2018_pop.csv",
  "clean_tr_dia_2_2017_pop.csv","clean_tr_dia_2_2018_pop.csv",
  "clean_tr_len_2017_pop.csv","clean_tr_len_2018_pop.csv"
)
files <- file.path(dir_raw, files)

# ── Helper to pick numeric trait column ───────────────────────────────────────
get_value_col <- function(dt) {
  num_cols <- names(dt)[sapply(dt, is.numeric)]
  ignore <- c("Year","Row","Column","rep","co_loc_PCR")
  setdiff(num_cols, ignore)[1]
}

# ── Read data + compute CV per genotype ───────────────────────────────────────
cv_tab <- rbindlist(lapply(files, function(f){
  dt <- fread(f)
  val_col <- get_value_col(dt)
  
  dt_summary <- dt %>%
    group_by(Genotype, Year) %>%
    summarise(
      Mean = mean(.data[[val_col]], na.rm = TRUE),
      SD   = sd(.data[[val_col]], na.rm = TRUE),
      CV   = 100*SD/Mean,
      .groups = "drop"
    ) %>%
    mutate(Trait = gsub("clean_|_pop.csv","",basename(f)))
  
  data.table(dt_summary)
}), fill = TRUE)

# ── Aggregate overall CV per trait × year ─────────────────────────────────────
cv_overall <- cv_tab %>%
  group_by(Trait, Year) %>%
  summarise(CV = mean(CV, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Trait = gsub("_2017|_2018","",Trait),
    Trait = recode(Trait,
                   "in_len"       = "INL",
                   "lsh_ang_1"    = "LSA",
                   "tr_len"       = "TL",
                   "tr_coni"      = "TC",
                   "tr_brz_prop"  = "TBZP",
                   "tr_brz_1"     = "TBZL",
                   "tr_brz_2"     = "TBZU",
                   "tr_dia_1"     = "TDB",
                   "tr_dia_2"     = "TDT"),
    Trait = factor(Trait,
                   levels = c("INL","LSA",
                              "TBZL","TBZU","TBZP",
                              "TDB","TDT",
                              "TC","TL")),
    Year  = factor(Year, levels = c(2017,2018))
  )

# ── Colour scheme like h² plot ───────────────────────────────────────────────
year_cols <- c("2017"="#1f78b4","2018"="#fdae61")

p <- ggplot(cv_overall, aes(x = fct_relevel(Trait, sort(unique(Trait))), y = CV, fill = Year)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.55) +
  scale_fill_manual(values = year_cols, name = "Year") +
  labs(x = NULL, y = "Coefficient of variation (%)") +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),
    
    axis.text.x        = element_text(angle = 30, hjust = 1, vjust = 1, size = 12),
    axis.title.y       = element_text(margin = margin(r = 8)),
    
    legend.position    = "right",
    legend.direction   = "vertical",
    legend.title       = element_text(size = 13),
    legend.text        = element_text(size = 11),
    legend.key.width   = unit(14, "pt"),
    legend.key.height  = unit(10, "pt"),
    
    plot.title         = element_blank(),
    plot.margin        = margin(t = 6, r = 14, b = 8, l = 8)
  )

# ── Save as SVG ─────────────────────────────────────────────────────────────
svg_path <- file.path(dir_stw, "Suppl_Fig_CV_raw_traits_genotypewise.svg")
ggsave(svg_path, plot = p, width = 12, height = 6, device = "svg")
cat("✔ Saved SVG to:", svg_path, "\n")
