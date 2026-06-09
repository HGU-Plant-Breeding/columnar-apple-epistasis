################################### Non-GxE GBLUP ASReml-R #######################################
# Close all devices and delete all variables
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files


#Install packages
library(tidyverse)
library(data.table)
library(ggplot2)
library(asreml)
# library(ASExtras4)
library(ASRgenomics)
library(AGHmatrix)
library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)
library(yardstick)
library(broom)
library(ggpmisc)
library(dplyr)
library(tidyr)
library(stringr)
library(snpReady)
library(impute)

# Define population colors
pop_colors <- c(
  bibon = "#1f78b4",  # blue
  gxr   = "#33a02c",  # green
  pxa   = "#e31a1c",  # red
  pxw   = "#ff7f00"   # orange
)


c130_matrix <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip/GT_filtered_numeric_transposed.csv", sep = ",")
n <- nrow(c130_matrix)
rownames(c130_matrix) <- c130_matrix$X
c130_matrix <- c130_matrix[, -(1)]
colnames(c130_matrix) <- gsub("^AX\\.", "AX-", colnames(c130_matrix))
c130_matrix <- as.matrix(c130_matrix)
c130_df <- as.data.frame(c130_matrix)

###########################################################     FOR MYB Carrier Subpopulations extra MAF filtering step   ############

#### Anthocyanin
setwd("C:/Users/nguevenc/Desktop/R_Working_Directory/Fruit traits/Anthocyan")
antho <- read.csv("AnthocyanContents.csv", sep = ";")
antho <- antho[, -c(1:2)]
antho$Genotype <- as.character(antho$Genotype)

antho <- antho[antho$Genotype %in% rownames(c130_df),]
antho <- antho[!duplicated(antho$Genotype),]

juice <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Fruit traits/Juice samples/Apfelsaft_Ergebnisse Nuri_Güvencli.csv" , sep = ";")
juice <- juice[juice$Genotype %in% rownames(c130_df),]
juice <- juice[!duplicated(juice$Genotype),]

rgb <- read.csv("C:/Users/nguevenc/Desktop/R_Working_Directory/Fruit traits/RGB/fruit_traits.csv", sep = ";")
rgb <- rgb[!grepl("PXL", rgb$imageID), ]
rgb$imageID <- sub("\\.jpg$", "", rgb$imageID)  ; colnames(rgb)[1] <- "Genotype"
rgb <- rgb[rgb$Genotype %in% rownames(c130_df),]

rgb <- rgb %>%
  group_by(Genotype) %>%
  mutate(
    avg_R_m = mean(avg_R, na.rm = TRUE),
    avg_G_m = mean(avg_G, na.rm = TRUE),
    avg_B_m = mean(avg_B, na.rm = TRUE),
    area_m = mean(area, na.rm = TRUE),
    h_w_ratio_m = mean(h_w_ratio, na.rm = TRUE)
  ) %>%
  ungroup()
rgb <- rgb[,-c(2:8)]

rgb <- rgb[!duplicated(rgb$Genotype), ]

genotypes_allfruits <- intersect(rgb$Genotype, juice$Genotype)
genotypes_antho <- intersect(antho$Genotype, genotypes_allfruits)


rgb_filtered <- rgb[rgb$Genotype %in% genotypes_allfruits,]
juice_filtered <- juice[juice$Genotype %in% genotypes_allfruits,]

allfruits <- left_join(juice_filtered, rgb_filtered, by = "Genotype")




##########  Heatmap Juice+RGB ###############

library(Hmisc)     # for rcorr()
library(corrplot)  # for corrplot()

# 1) pull out only the numeric columns
num_df <- allfruits[ , setdiff(names(allfruits), "Genotype") ]

# 2) get Pearson r and p‐values
rc       <- rcorr(as.matrix(num_df), type="pearson")
corr_mat <- rc$r
p_mat    <- rc$P

# 3) make nice labels
orig <- colnames(corr_mat)
new  <- ifelse(
  grepl("_m$", orig),
  paste0(sub("_m$", "", orig), " (fruit RGB)"),
  paste0(orig,         " (fruit juice)")
)
rownames(corr_mat) <- new
colnames(corr_mat) <- new
rownames(p_mat)    <- new
colnames(p_mat)    <- new

# 4) mask lower‐triangle + diagonal of p_mat
#    so that melt(p.mat, na.rm=TRUE) will produce exactly the same
#    number of points as melt(corrMat, na.rm=TRUE)
p_mat[lower.tri(p_mat, diag=TRUE)] <- NA

# 5) draw it
corrplot(
  corr_mat,
  method     = "circle",                        # circles scaled by |r|
  type       = "upper",                         # strict upper‐triangle only
  order      = "original",
  diag       = FALSE,                           # drop the diagonal
  col        = colorRampPalette(c("navy","white","firebrick3"))(200),
  tl.col     = "black",  tl.srt = 45,           # trait labels in black
  addCoef.col= "black",  number.cex = 0.7,      # r‐values printed in black
  p.mat      = p_mat,                           # your masked p‐matrix
  sig.level  = 0.05,                            # α = 0.05
  insig      = "blank",                         # blank out p ≥ 0.05
  cl.pos     = "n",                             # drop the color legend
  main       = paste(
    strwrap("Correlation of juice and fruit RGB traits",
            width = 40),
    collapse = "\n"
  ),
  mar        = c(0,0,2,0)                       # room for the title
)



################ 

anthofruits <- right_join(allfruits, antho, by = "Genotype")
anthofruits <- anthofruits[!is.na(anthofruits$density)& !is.na(anthofruits$colour520),]
colnames(anthofruits)[19] <- "anthocyanin_m"



library(Hmisc)     # for rcorr()
library(corrplot)  # for corrplot()

# 1) pull out only the numeric columns
num_df <- anthofruits[ , setdiff(names(anthofruits), "Genotype") ]

# 2) get Pearson r and p‐values
rc       <- rcorr(as.matrix(num_df), type="pearson")
corr_mat <- rc$r
p_mat    <- rc$P

# 3) make nice labels
orig <- colnames(corr_mat)
new  <- ifelse(
  grepl("_m$", orig),
  paste0(sub("_m$", "", orig), " (fruit RGB)"),
  paste0(orig,         " (fruit juice)")
)
new[18] <- "anthocyanin"

rownames(corr_mat) <- new
colnames(corr_mat) <- new
rownames(p_mat)    <- new
colnames(p_mat)    <- new


# 4) mask lower‐triangle + diagonal of p_mat
#    so that melt(p.mat, na.rm=TRUE) will produce exactly the same
#    number of points as melt(corrMat, na.rm=TRUE)
p_mat[lower.tri(p_mat, diag=TRUE)] <- NA




# 5) draw it
# First, plot without the 'main' title
corrplot(
  corr_mat,
  method       = "circle",                        # circles scaled by |r|
  type         = "upper",                          # strict upper‐triangle only
  order        = "original",
  diag         = FALSE,                            # drop the diagonal
  col          = colorRampPalette(c("navy","white","firebrick3"))(200),
  tl.col       = "black",  tl.srt = 45,            # trait labels in black
  # Removed addCoef.col to hide numbers
  p.mat        = p_mat,                            # your masked p‐matrix
  sig.level    = 0.05,                             # α = 0.05
  insig        = "blank",                          # blank out p ≥ 0.05
  cl.pos       = "r",                              # Changed to "r" to show color legend on the right
  mar          = c(0,0,0,0)                        # Reset margins, or adjust as needed for overall plot
)

# Then add the title using the 'title()' function for more control
title(
  main = "Correlation of juice and fruit RGB traits with anthocyanin content",
  line = -0.5, # Experiment with this value (e.g., -1, 0, 0.5)
  cex.main = 1.2 # Adjust title size if needed
)

