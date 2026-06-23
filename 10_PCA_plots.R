# ── Setup ─────────────────────────────────────────────────────────────
rm(list = ls()); graphics.off(); closeAllConnections()
library(car)
library(tidyverse)
library(svglite)

BASE_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
IN_DIR   <- file.path(BASE_DIR, "deregression_output")
setwd(BASE_DIR)

population_colors <- c(pxw = "#ff7f00",
                       pxa = "#e31a1c",
                       gxr = "#33a02c")

# ── Read and prepare data ─────────────────────────────────────────────
files <- list.files(IN_DIR, pattern="^dGEBV_.*_pop\\.csv$", full.names=TRUE)

df <- map_dfr(files, read.csv, check.names = FALSE) %>%
  transmute(
    Year     = as.integer(Year),
    Trait    = recode(as.character(Trait),
                      "in_len"      = "INL",
                      "lsh_ang_1"   = "LSA",
                      "tr_brz_prop" = "TBZP",
                      "tr_coni"     = "TC",
                      "tr_len"      = "TL"),
    Genotype = as.character(Genotype),
    pop      = factor(pop, levels = c("pxw","pxa","gxr")),
    co_lab   = factor(ifelse(co_loc %in% c(1,"1","Mutant"),
                             "Mutant","Wildtype"),
                      levels = c("Wildtype","Mutant")),
    value    = as.numeric(dGEBV + pop_BLUE)
  ) %>%
  group_by(Year, Trait) %>%
  mutate(z = scale(value)[,1]) %>%
  ungroup()

# ── run PCA and clean inputs ──────────────────────────────────
run_pca_year <- function(dat_year) {
  
  wide <- dat_year %>%
    select(Genotype, Trait, z) %>%
    pivot_wider(names_from = Trait, values_from = z) %>%
    drop_na()
  
  meta <- dat_year %>%
    select(Genotype, pop, co_lab) %>%
    distinct()
  
  X <- wide %>% column_to_rownames("Genotype")
  meta <- meta %>% filter(Genotype %in% rownames(X)) %>%
    arrange(match(Genotype, rownames(X)))
  
  pr <- prcomp(X, center = FALSE, scale. = FALSE)
  
  scores <- as.data.frame(pr$x[,1:2])
  scores$Genotype <- rownames(scores)
  
  list(pr = pr,
       scores = left_join(scores, meta, by="Genotype"),
       loadings = as.data.frame(pr$rotation[,1:2]),
       var_exp = 100 * pr$sdev^2 / sum(pr$sdev^2))
}

# ── variance partitioning ─────────────────────────────────────
variance_partition <- function(scores, pc_name) {
  
  mod <- lm(reformulate(c("co_lab","pop"), response = pc_name),
            data = scores)
  
  a <- Anova(mod, type = 2)
  
  ss_total <- sum((scores[[pc_name]] -
                     mean(scores[[pc_name]]))^2)
  
  tibble(
    PC = pc_name,
    Term = rownames(a),
    Semi_R2 = a$`Sum Sq` / ss_total
  )
}

# ── biplot ────────────────────────────────────────────────────
plot_biplot <- function(res, year) {
  
  scores   <- res$scores
  loadings <- res$loadings
  var_exp  <- res$var_exp
  
  xr <- range(scores$PC1)
  yr <- range(scores$PC2)
  scale_arrow <- 0.8 * min(diff(xr), diff(yr))
  loadings <- loadings * scale_arrow
  loadings$Trait <- rownames(loadings)
  
  p <- ggplot(scores,
              aes(PC1, PC2,
                  color = pop,
                  shape = co_lab)) +
    geom_point(size = 3) +
    scale_color_manual(values = population_colors,
                       name = "Population") +
    scale_shape_manual(values = c(1,16),
                       name = "Co locus") +
    geom_segment(data = loadings,
                 aes(x=0,y=0,
                     xend=PC1,yend=PC2),
                 inherit.aes = FALSE,
                 arrow = arrow(length=unit(0.18,"cm"))) +
    geom_text(data = loadings,
              aes(PC1, PC2, label = Trait),
              inherit.aes = FALSE,
              fontface = "bold",
              vjust = -0.6) +
    labs(title = paste0("PCA (",year,")"),
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face="bold", hjust=0.5))
  
  svglite(file.path(BASE_DIR,
                    paste0("pca_biplot_",year,".svg")),
          width=8, height=7)
  print(p)
  dev.off()
}

# ── Main loop ─────────────────────────────────────────────────────────
results_all <- list()

for (yy in sort(unique(df$Year))) {
  
  dat_year <- filter(df, Year == yy)
  
  res <- run_pca_year(dat_year)
  
  plot_biplot(res, yy)
  
  vp1 <- variance_partition(res$scores, "PC1")
  vp2 <- variance_partition(res$scores, "PC2")
  
  results_year <- bind_rows(vp1, vp2) %>%
    mutate(Year = yy)
  
  results_all[[as.character(yy)]] <- results_year
}

final_results <- bind_rows(results_all)

write.csv(final_results,
          file.path(BASE_DIR, "pc_variance_partitioning.csv"),
          row.names = FALSE)

print(final_results)