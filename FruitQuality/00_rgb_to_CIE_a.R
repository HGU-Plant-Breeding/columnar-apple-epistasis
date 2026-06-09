# ==============================================================================
# FINAL APPLE REDNESS (a*) EXTRACTION
# Confirmed: 2025 = Semicolon (;) | 2024 = Comma (,)
# ==============================================================================

library(colorspace)

# 1. SET PATHS
path_2025 <- r"(C:\Users\nguevenc\Desktop\R_Working_Directory\FruitQuality\RGB pics\project-31-at-2026-03-20-12-11-27c140f9\inference_all\apple_measurements.csv)"
path_2024 <- r"(C:\Users\nguevenc\Desktop\R_Working_Directory\FruitQuality\RGB pics\project-29-at-2025-03-14-13-00-f7717b10\inference_2024\RGB_2024.csv)"
output_dir <- r"(C:\Users\nguevenc\Desktop\R_Working_Directory\FruitQuality)"

# 2. CONVERSION FUNCTION
# Calculates CIE a* from avg_R, avg_G, avg_B
calculate_a <- function(R, G, B) {
  sapply(seq_along(R), function(i) {
    if (any(is.na(c(R[i], G[i], B[i])))) return(NA)
    # Convert sRGB (0-255) to LAB and return the 'A' coordinate
    lab <- as(sRGB(R[i]/255, G[i]/255, B[i]/255), "LAB")
    round(coords(lab)[1, "A"], 2)
  })
}

# 3. PROCESS 2025 (SEMICOLON SEP)
if (file.exists(path_2025)) {
  df25 <- read.csv(path_2025, sep = ";", header = TRUE)
  df25$CIE_a <- calculate_a(df25$avg_R, df25$avg_G, df25$avg_B)
  
  write.csv(df25, file.path(output_dir, "RGB_traits_2025.csv"), row.names = FALSE)
  cat("2025 Success: Processed", nrow(df25), "apples (using Semicolon sep).\n")
}

# 4. PROCESS 2024 (COMMA SEP)
if (file.exists(path_2024)) {
  df24 <- read.csv(path_2024, sep = ",", header = TRUE)
  df24$CIE_a <- calculate_a(df24$avg_R, df24$avg_G, df24$avg_B)
  
  write.csv(df24, file.path(output_dir, "RGB_traits_2024.csv"), row.names = FALSE)
  cat("2024 Success: Processed", nrow(df24), "apples (using Comma sep).\n")
}

cat("\nAnalysis complete. Files are in the 'FruitQuality' folder.\n")