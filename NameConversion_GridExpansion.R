library(dplyr)
library(tidyr)

# --- 1. LOAD DATASETS ---
PI17_GDR   <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/2017_GDR.csv", sep = ";")
PI17_LA    <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/2017_LA.csv", sep = ";")
PI18_LA    <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/2018_LA.csv", sep = ";")
PI18_GDR   <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/2018_GDR.csv", sep = ";")

Names      <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/NameConversion.csv", sep = ";")
co_loc_PCR <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/PCR_confirmation_and_GWAS/co_status.csv", sep = ";")

# --- 2. REPLACE NAMES & SET NON-MATCHES TO NA ---
update_genotype_names <- function(df, name_ref) {
  df %>%
    left_join(name_ref, by = c("Genotype" = "KENNAME")) %>%
    mutate(Genotype = MYNAME) %>% 
    select(-MYNAME)
}

PI17_GDR <- update_genotype_names(PI17_GDR, Names)
PI17_LA  <- update_genotype_names(PI17_LA, Names)
PI18_GDR <- update_genotype_names(PI18_GDR, Names)
PI18_LA  <- update_genotype_names(PI18_LA, Names)

# --- 3. EXPAND SPATIAL GRIDS (3 REPS PER PLOT) ---
expand_to_spatial_grid <- function(df) {
  df %>%
    group_by(Loc) %>%
    reframe(
      expand_grid(
        R = 1:max(R, na.rm = TRUE),
        N = 1:max(N, na.rm = TRUE),
        rep = 1:3
      )
    ) %>%
    left_join(df, by = c("Loc", "R", "N")) %>%
    arrange(Loc, R, N, rep)
}

PI17_GDR_grid <- expand_to_spatial_grid(PI17_GDR)
PI17_LA_grid  <- expand_to_spatial_grid(PI17_LA)
PI18_GDR_grid <- expand_to_spatial_grid(PI18_GDR)
PI18_LA_grid  <- expand_to_spatial_grid(PI18_LA)

# --- 4. LOAD OLD DATA & DEFINE PHENOTYPE COLUMNS ---
df_17_old <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/df_17.csv", sep = ",")
df_18_old <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/df_18.csv", sep = ",")

cols_to_move <- c("Genotype", "rep", "tr_len", "tr_brz_1", "tr_brz_2", "tr_dia_1", 
                  "tr_dia_2", "ssh_count", "lsh_count", "in_len", "tsh_count", 
                  "tr_brz_len", "tr_dia_m", "tr_slend", "tr_coni", "tr_vol", 
                  "ssh_prop", "lsh_prop", "tr_brz_prop", "tr_brz_dens", 
                  "tr_len_incr", "lsh_len", "lsh_cord", "lsh_ang_1", "lsh_ang_2", 
                  "lsh_dia_1", "lsh_dia_2", "lsh_ang_cord", "lsh_dia_m", 
                  "lsh_slend", "lsh_coni", "lsh_ang_bend", "lsh_cord_bend", "pop")

# Deduplicate old measurement data
df_17_clean <- df_17_old %>% select(all_of(cols_to_move)) %>% distinct(Genotype, rep, .keep_all = TRUE)
df_18_clean <- df_18_old %>% select(all_of(cols_to_move)) %>% distinct(Genotype, rep, .keep_all = TRUE)

# --- 5. JOIN PHENOTYPES & CLEAN UP DUPLICATES ---
finalize_dataset <- function(grid_df, old_subset) {
  grid_df %>%
    left_join(old_subset, by = c("Genotype", "rep")) %>%
    mutate(Genotype = ifelse(is.na(pop), NA_character_, Genotype)) %>%
    group_by(Loc, R, N, rep) %>%
    slice(1) %>%
    ungroup()
}

PI17_GDR_final <- finalize_dataset(PI17_GDR_grid, df_17_clean)
PI17_LA_final  <- finalize_dataset(PI17_LA_grid, df_17_clean)
PI18_GDR_final <- finalize_dataset(PI18_GDR_grid, df_18_clean)
PI18_LA_final  <- finalize_dataset(PI18_LA_grid, df_18_clean)

# --- 6. RENAME Loc ENTRIES IN GDR DATASETS ---
PI17_GDR_final <- PI17_GDR_final %>% mutate(Loc = gsub("La_18_", "GDR_18_", Loc, ignore.case = TRUE))
PI18_GDR_final <- PI18_GDR_final %>% mutate(Loc = gsub("La_18_", "GDR_18_", Loc, ignore.case = TRUE))

# --- 7. MERGE, RENAME, REORDER AND JOIN STATUS ---

# Deduplicate PCR status
co_loc_PCR_clean <- co_loc_PCR %>% distinct(Genotype, .keep_all = TRUE)

# Processing 2017
df_17_final <- bind_rows(
  PI17_GDR_final %>% mutate(Year = 2017),
  PI17_LA_final  %>% mutate(Year = 2017)
) %>%
  left_join(co_loc_PCR_clean, by = "Genotype") %>%
  mutate(Plot_RC_ID = paste0(Loc, "_R", R, "N", N)) %>%
  rename(Row = R, Column = N) %>%
  # Move requested columns to the front
  select(pop, co_loc, Plot_RC_ID, Year, Loc, Row, Column, rep, Genotype, everything())

# Processing 2018
df_18_final <- bind_rows(
  PI18_GDR_final %>% mutate(Year = 2018),
  PI18_LA_final  %>% mutate(Year = 2018)
) %>%
  left_join(co_loc_PCR_clean, by = "Genotype") %>%
  mutate(Plot_RC_ID = paste0(Loc, "_R", R, "N", N)) %>%
  rename(Row = R, Column = N) %>%
  # Move requested columns to the front
  select(pop, co_loc, Plot_RC_ID, Year, Loc, Row, Column, rep, Genotype, everything())

# --- 8. EXPORT ---

# Set Output Directory
out_path <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough/"

# Save as CSVs
write.csv(df_17_final, paste0(out_path, "df_17_final_w_outliers_30traits.csv"), row.names = FALSE)
write.csv(df_18_final, paste0(out_path, "df_18_final_w_outliers_30traits.csv"), row.names = FALSE)

print("Final dataframes created with Row/Column naming and reordered columns. Files saved.")