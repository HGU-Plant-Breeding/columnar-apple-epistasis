# ─────────────────────────────────────────────────────────────
# 0) Reset session
rm(list = ls()); graphics.off(); closeAllConnections()

# ─────────────────────────────────────────────────────────────
# 1) Libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(svglite)

# ─────────────────────────────────────────────────────────────
# 2) Paths
GT_path <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip/GT_filtered_numeric_transposed.csv"

pheno_file <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough/outlier_removed_raw/clean_in_len_2017_pop.csv"

OUT_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"

# ─────────────────────────────────────────────────────────────
# 3) Population colours

population_colors <- c(pxw = "#ff7f00",
                       pxa = "#e31a1c",
                       gxr = "#33a02c")

# ─────────────────────────────────────────────────────────────
# 4) Load marker matrix

GT <- fread(GT_path, data.table = FALSE)
colnames(GT)[1] <- "Genotype"

M <- as.matrix(GT[, -1])
rownames(M) <- GT$Genotype

# ─────────────────────────────────────────────────────────────
# 5) Compute heterozygosity per genotype

het <- data.frame(
  Genotype = rownames(M),
  H = rowMeans(M == 1, na.rm = TRUE)
)

# ─────────────────────────────────────────────────────────────
# 6) Load population labels

pop_df <- read.csv(pheno_file) |>
  select(Genotype, pop) |>
  distinct()

pop_df$pop <- tolower(pop_df$pop)

# ─────────────────────────────────────────────────────────────
# 7) Merge and keep the three populations

het <- het |>
  left_join(pop_df, by = "Genotype") |>
  filter(pop %in% c("gxr", "pxw", "pxa"))

het$pop <- factor(het$pop, levels = c("pxw", "pxa", "gxr"))

# ─────────────────────────────────────────────────────────────
# 8) Create violin plot

p <- ggplot(het, aes(x = pop, y = H, fill = pop)) +
  
  geom_violin(trim = FALSE,
              color = "black",
              linewidth = 0.4) +
  
  scale_fill_manual(values = population_colors,
                    name = "Population") +
  
  labs(x = "Population",
       y = "Individual heterozygosity") +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# ─────────────────────────────────────────────────────────────
# 9) Save SVG

svglite(file.path(OUT_DIR, "heterozygosity_populations.svg"),
        width = 6,
        height = 4)

print(p)

dev.off()