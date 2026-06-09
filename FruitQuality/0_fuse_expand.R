#' -----------------------------------------------------------------------------
#' Script: Fruit Quality Data Preparation (Strict Field Boundaries)
#' Purpose: Create rectangular spatial grids for ASReml with hard-coded limits:
#'          Fuchsberg (Max Row 3), Lachacker (Max Row 6), Bibon (Max Row 2).
#' -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

base_path <- "C:/Users/nguevenc/Desktop/R_Working_Directory/FruitQuality/"

# --- 1. Load and Standardize Data ---
load_standardized <- function(filename) {
  path <- paste0(base_path, filename)
  if(!file.exists(path)) return(NULL)
  
  df <- read.csv(path, sep = ";", na.strings = c("NA", ""), check.names = FALSE)
  
  df %>%
    mutate(across(-any_of("Genotype"), ~{
      if(is.character(.)) as.numeric(gsub(",", ".", .)) else .
    }))
}

# Load Phenotypes
df_firmness <- bind_rows(
  load_standardized("firmness2024.csv") %>% mutate(year = 2024),
  load_standardized("firmness2025.csv") %>% mutate(year = 2025) %>% 
    group_by(Genotype) %>% mutate(rep = row_number()) %>% ungroup()
) %>% distinct(Genotype, year, rep, .keep_all = TRUE)

df_rgb <- bind_rows(
  load_standardized("RGB_traits_2024.csv") %>% mutate(year = 2024),
  load_standardized("RGB_traits_2025.csv") %>% mutate(year = 2025)
) %>% rename(rep = apple_number) %>% distinct(Genotype, year, rep, .keep_all = TRUE)

df_antho <- bind_rows(load_standardized("Anthocyanin_2024.csv") %>% mutate(year = 2024),
                      load_standardized("Anthocyanin_2025.csv") %>% mutate(year = 2025)) %>% 
  distinct(Genotype, year, .keep_all = TRUE)

df_juice <- bind_rows(load_standardized("Juice2024.csv") %>% mutate(year = 2024),
                      load_standardized("Juice2025.csv") %>% mutate(year = 2025)) %>% 
  distinct(Genotype, year, .keep_all = TRUE)

# --- 2. Build and Filter Spatial Skeleton ---
all_gts <- unique(c(df_antho$Genotype, df_juice$Genotype, df_firmness$Genotype, df_rgb$Genotype))
all_gts <- all_gts[!is.na(all_gts)]

spatial_master <- data.frame(Genotype = as.character(all_gts)) %>%
  mutate(
    location = case_when(
      str_detect(Genotype, "^B") ~ "bibon",
      str_detect(Genotype, "^F") ~ "fuchsberg",
      str_detect(Genotype, "^L") ~ "lachacker",
      TRUE ~ NA_character_
    ),
    row = as.numeric(str_extract(Genotype, "(?<=R)\\d+(?=N)")),
    column = as.numeric(str_extract(Genotype, "(?<=N)\\d+")),
    pop = case_when(
      location == "fuchsberg" ~ "gxr",
      location == "bibon" ~ "bibon",
      location == "lachacker" & ((row == 1 & column >= 69) | (row %in% 2:3) | 
                                   (row == 4 & column <= 44)) ~ "pxw",
      location == "lachacker" ~ "pxa",
      TRUE ~ NA_character_
    )
  ) %>%
  # PRE-FILTER: Remove malformed IDs or those exceeding known orchard rows
  filter(!is.na(location), !is.na(row), !is.na(column)) %>%
  filter(!(location == "fuchsberg" & row > 3),
         !(location == "lachacker" & row > 6),
         !(location == "bibon" & row > 2)) %>%
  distinct(location, row, column, .keep_all = TRUE)

# --- 3. Isolated Expansion (Strict Limits) ---
full_grid <- spatial_master %>%
  group_split(location) %>%
  map_dfr(function(loc_df) {
    loc_name <- unique(loc_df$location)
    
    # Use specified limits for rows to be absolutely certain
    max_r <- case_when(
      loc_name == "fuchsberg" ~ 3,
      loc_name == "lachacker" ~ 6,
      loc_name == "bibon"     ~ 2,
      TRUE ~ max(loc_df$row)
    )
    max_c <- max(loc_df$column)
    
    expand_grid(
      location = loc_name,
      year     = c(2024, 2025),
      row      = 1:max_r,
      column   = 1:max_c
    )
  }) %>%
  left_join(spatial_master, by = c("location", "row", "column"))

# --- 4. Final Merges ---
df_pcr <- load_standardized("MYB_PCR.csv") %>% distinct(Genotype, .keep_all = TRUE)
df_co  <- load_standardized("co_status.csv") %>% distinct(Genotype, .keep_all = TRUE)

finalize <- function(df) {
  df %>%
    left_join(df_pcr, by = "Genotype") %>%
    left_join(df_co,  by = "Genotype") %>%
    mutate(myb = case_when(pop %in% c("gxr", "pxa") ~ 0, TRUE ~ myb),
           co_loc = case_when(pop == "bibon" ~ NA_real_, TRUE ~ co_loc))
}

replicated_full  <- full_grid %>% expand_grid(rep = 1:8) %>%
  left_join(df_firmness, by = c("Genotype", "year", "rep")) %>%
  left_join(df_rgb,      by = c("Genotype", "year", "rep")) %>%
  finalize()

juice_antho_full <- full_grid %>%
  left_join(df_antho, by = c("Genotype", "year")) %>%
  left_join(df_juice, by = c("Genotype", "year")) %>%
  finalize()

# --- 5. Export and Verify ---
write.csv(replicated_full,  paste0(base_path, "replicated_traits_FULL.csv"), row.names = FALSE, na = "NA")
write.csv(juice_antho_full, paste0(base_path, "juice_antho_traits_FULL.csv"), row.names = FALSE, na = "NA")

cat("\n--- FINAL GRID VERIFICATION ---\n")
full_grid %>% 
  group_by(location) %>% 
  summarise(Max_Row = max(row), Max_Col = max(column)) %>% 
  print()