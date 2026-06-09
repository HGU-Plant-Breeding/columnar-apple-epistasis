rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()
options(stringsAsFactors = FALSE, scipen = 999)

library(data.table)
library(ggplot2)
library(scales)
library(svglite)

BASE <- "C:/Users/nguevenc/Desktop/R_Working_Directory"
STW  <- file.path(BASE, "Lidar", "Single_trait_walkthrough")
IN   <- file.path(STW,  "sign_Kernel_VCs")
TRAITS <- c("in_len","lsh_ang_1","tr_brz_prop","tr_coni","tr_len")
trait_labs <- c(
  "in_len"      = "INL",
  "lsh_ang_1"   = "LSA",
  "tr_brz_prop" = "TBZP",
  "tr_coni"     = "TC",
  "tr_len"      = "TL")

YEARS  <- c(2017, 2018)
COMPS  <- c("A_sig","A_bg","Field spatial","Residual")
files <- list.files(IN, pattern="^vc_signKernel_.*\\.csv$", full.names=TRUE)


vc <- rbindlist(lapply(files, fread), fill = TRUE)
vc <- vc[Component %in% COMPS]

vc[, Prop_plot := as.numeric(gsub(",", ".", as.character(Prop), fixed=TRUE)) / 100]

vc[, Trait     := factor(Trait, levels=TRAITS)]
vc[, Year      := factor(as.integer(Year), levels=YEARS)]
vc[, Component := factor(Component, levels=COMPS)]
vc[, Model     := factor(Model)]

vc <- vc[is.finite(Prop_plot) & Prop_plot >= 0 & Prop_plot <= 1 & !is.na(Trait) & !is.na(Year) & !is.na(Model)]

fill_vals <- c("A_sig"="#5a68a1","A_bg"="#a6cee3","Field spatial"="#cc79a7","Residual"="#e69f00")
fill_labs <- c("A_sig"="Additive (GWAS-significant SNP kernel)",
               "A_bg"="Additive (background SNP kernel)",
               "Field spatial"="Field spatial",
               "Residual"="Residual")

p <- ggplot(vc, aes(Trait, Prop_plot, fill=Component)) +
  geom_col(width=0.75) +
  facet_grid(Year ~ Model) +
  scale_x_discrete(labels = trait_labs) +
  scale_fill_manual(values=fill_vals, labels=fill_labs, name="Variance component") +
  scale_y_continuous(labels=percent_format(accuracy=1), expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  labs(x=NULL, y=NULL) +
  theme_bw(base_size=16) +
  theme(
    axis.text.x        = element_text(angle=45, hjust=1, vjust=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill="gray95", color="black", linewidth=0.4),
    strip.text         = element_text(face="bold"),
    legend.position    = "right",
    plot.margin        = margin(t=6, r=14, b=8, l=8))
print(p)

out_svg <- file.path(STW, "vc_signKernel.svg")
ggsave(out_svg, p, width=12, height=7, device=svglite)