# ─────────────────────────────────────────────────────────────────────
# 0) Libraries
# ─────────────────────────────────────────────────────────────────────
rm(list = ls(all = TRUE)); graphics.off(); closeAllConnections()

library(data.table)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(viridisLite)
library(grid)
library(purrr)
library(svglite)   # <─ for SVG output

# ─────────────────────────────────────────────────────────────────────
# 1) Paths, traits & global switches
# ─────────────────────────────────────────────────────────────────────

BASE_STW      <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
TRAITS        <- c("INL", "LSA", "TBZP", "TC", "TL")
YEARS         <- c(2017, 2018)

SNP_MAP_FILE  <- "C:/Users/nguevenc/Desktop/R_Working_Directory/SNP Chip/GM_v11.csv"
CL_QTL_CSV    <- file.path(BASE_STW, "Coupel_Ledru_2022_QTLs.csv")
OUT_DIR       <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar"

OUT_CIRCOS_SVG          <- file.path(OUT_DIR, "circos.svg")
OUT_DENSITY_LEGEND_SVG  <- file.path(OUT_DIR, "legend_density.svg")
OUT_DISTANCE_LEGEND_SVG <- file.path(OUT_DIR, "legend_distance.svg")

BIN_SIZE      <- 1e6     # 1 Mb windows
VALID_CHRS    <- as.character(1:17)

# PCR Co-locus and Co-window definition
pcr_chr       <- "10"
pcr_pos       <- round(mean(c(27955554, 27955795)))
CO_WIN_START  <- 20000000L
CO_WIN_END    <- 36000000L

# Significance threshold (GWAS + Wald FDR)
FDR_CUTOFF    <- 0.05

# Trait colours
TRAIT_COLS <- c(
  "TL"   = "#1b9e77",
  "INL"  = "#d95f02",
  "TC"   = "#7570b3",
  "LSA"  = "#e7298a",
  "TBZP" = "#66a61e"
)

# CL track style
CL_LWD    <- 1.05
CL_BG_COL <- "#FFFFFF"

# ─────────────────────────────────────────────────────────────────────
# 2) Load all GWAS & Wald files into long format
# ─────────────────────────────────────────────────────────────────────

files_grid <- expand.grid(Trait = TRAITS, Year = YEARS, stringsAsFactors = FALSE)

gwas_list <- lapply(seq_len(nrow(files_grid)), function(i) {
  tr  <- files_grid$Trait[i]
  yr  <- files_grid$Year[i]
  dir_y  <- file.path(BASE_STW, sprintf("BLINK_dGEBV_%d", yr))
  gwas_f <- file.path(dir_y, sprintf("GAPIT.Association.GWAS_Results.BLINK.%s(NYC).csv", tr))
  df <- fread(gwas_f, data.table = FALSE, check.names = FALSE)
  df$Trait <- tr
  df$Year  <- yr
  df
})

wald_list <- lapply(seq_len(nrow(files_grid)), function(i) {
  tr   <- files_grid$Trait[i]
  yr   <- files_grid$Year[i]
  dir_y  <- file.path(BASE_STW, sprintf("BLINK_dGEBV_%d", yr))
  wald_f <- file.path(dir_y, sprintf("wald_%s_%d_SNPfixed_GBLUP.csv", tr, yr))
  df <- fread(wald_f, data.table = FALSE, check.names = FALSE)
  df$Trait <- tr
  df$Year  <- yr
  df
})

gwas_df <- rbindlist(gwas_list, fill = TRUE)
wald_df <- rbindlist(wald_list,  fill = TRUE)

# ─────────────────────────────────────────────────────────────────────
# 3) SNP map + CL QTLs
# ─────────────────────────────────────────────────────────────────────

snp_map <- fread(SNP_MAP_FILE, sep = ";")
snp_map$chrom <- as.character(snp_map$chrom)
snp_map <- snp_map %>% filter(chrom %in% VALID_CHRS)

genotype_map <- snp_map %>%
  group_by(chrom) %>%
  summarise(start = 0, end = max(pos, na.rm = TRUE), .groups = "drop") %>%
  rename(chr = chrom) %>%
  mutate(chr_num = as.integer(chr)) %>%
  arrange(chr_num, chr) %>%
  select(-chr_num) %>%
  as.data.frame()

CL_QTLs <- fread(CL_QTL_CSV, sep = ";")

# ─────────────────────────────────────────────────────────────────────
# 4) SNP density bins, excluding Co-window from density
# ─────────────────────────────────────────────────────────────────────

snp_map_density <- snp_map %>%
  filter(!(chrom == pcr_chr & pos >= CO_WIN_START & pos <= CO_WIN_END))

bin_one_chr <- function(chr, end_pos, df, bin = BIN_SIZE) {
  pos_vec <- df$pos[df$chrom == chr]
  brks <- seq(0, end_pos, by = bin)
  if (tail(brks, 1) < end_pos) brks <- c(brks, end_pos)
  h <- hist(pos_vec, breaks = brks, plot = FALSE)
  data.frame(
    chr   = chr,
    start = head(h$breaks, -1),
    end   = tail(h$breaks, -1),
    count = as.numeric(h$counts)
  )
}

density_bed <- do.call(
  rbind,
  lapply(seq_len(nrow(genotype_map)), function(i) {
    bin_one_chr(genotype_map$chr[i], genotype_map$end[i], snp_map_density, BIN_SIZE)
  })
)

col_fun <- colorRamp2(c(0, 25, 60, 100), viridisLite::magma(4, direction = -1))

# ─────────────────────────────────────────────────────────────────────
# 5) GWAS ticks per trait & year
# ─────────────────────────────────────────────────────────────────────

gwas_df <- gwas_df %>%
  mutate(
    Chr  = as.character(Chr),
    Pos  = as.numeric(Pos),
    Year = as.integer(Year)
  ) %>%
  filter(
    Chr   %in% VALID_CHRS,
    Trait %in% TRAITS,
    Year  %in% YEARS
  )

gwas_sig <- gwas_df %>%
  group_by(Trait, Year) %>%
  mutate(FDR_BH = p.adjust(P.value, method = "BH")) %>%
  ungroup() %>%
  filter(FDR_BH <= FDR_CUTOFF)

ticks_2017 <- gwas_sig %>% filter(Year == 2017) %>% transmute(chr = Chr, pos = Pos, Trait)
ticks_2018 <- gwas_sig %>% filter(Year == 2018) %>% transmute(chr = Chr, pos = Pos, Trait)

trait_cols_use <- sapply(TRAIT_COLS, adjustcolor, alpha.f = 0.95, USE.NAMES = TRUE)

# ─────────────────────────────────────────────────────────────────────
# 6) Coupel-Ledru ticks + distance to nearest GWAS hit
# ─────────────────────────────────────────────────────────────────────

cl_ticks <- CL_QTLs %>%
  transmute(
    CL_SNP = as.character(SNP),
    chr    = as.character(Chr),
    pos    = as.numeric(`physical position (bp)`)
  ) %>% filter(chr %in% VALID_CHRS)

sig_pos <- gwas_sig %>% distinct(chr = Chr, pos = Pos)

get_min_dist <- function(chr, pos, sig_df) {
  pos_sig <- sig_df$pos[sig_df$chr == chr]
  min(abs(pos_sig - pos))
}

cl_ticks$min_dist_bp <- mapply(
  get_min_dist,
  chr = cl_ticks$chr,
  pos = cl_ticks$pos,
  MoreArgs = list(sig_df = sig_pos)
)

cl_ticks <- cl_ticks %>% mutate(
  dist_kb  = min_dist_bp / 1e3,
  dist_cap = pmin(pmax(dist_kb, 150), 500)
)

dist_col_fun <- circlize::colorRamp2(c(150, 325, 500), c("#160b39", "#b21a6b", "#f98400"))


dens_legend <- ComplexHeatmap::Legend(
  col_fun   = col_fun,
  title     = "SNP density (1 Mb)",
  at        = c(0, 25, 60, 100),
  direction = "vertical",
  title_gp  = grid::gpar(fontface = "bold")
)

dist_legend <- ComplexHeatmap::Legend(
  col_fun   = dist_col_fun,
  title     = "Distance to nearest\nGWAS-hit (kb)",
  at        = c(150, 250, 350, 500),
  labels    = c("\u2264150", "250", "350", "\u2265500"),
  direction = "vertical",
  title_gp  = grid::gpar(fontface = "bold")
)

# ─────────────────────────────────────────────────────────────────────
# 7) Wald: Co:SNP and SNP:SNP interactions
# ─────────────────────────────────────────────────────────────────────

wald_df <- wald_df %>%
  mutate(
    Effect_clean = gsub("AX\\.", "AX-", gsub(" ", "", Effect)),
    Year         = as.integer(Year)
  ) %>%
  filter(
    Trait %in% TRAITS,
    Year  %in% YEARS,
    !is.na(FDR_BH),
    FDR_BH <= FDR_CUTOFF,
    grepl("AX-", Effect_clean)
  )

wald_co <- wald_df %>%
  filter(grepl("co_loc", Effect_clean)) %>%
  mutate(SNP = sub(".*(AX-[0-9]+).*", "\\1", Effect_clean)) %>%
  select(Trait, Year, SNP, FDR_BH)

wald_snpsnp <- wald_df %>%
  filter(
    !grepl("co_loc", Effect_clean),
    grepl("AX-[0-9]+:AX-[0-9]+", Effect_clean)
  ) %>%
  mutate(
    pair = sub(".*(AX-[0-9]+:AX-[0-9]+).*", "\\1", Effect_clean),
    SNP1 = sub(":(AX-[0-9]+)$", "", pair),
    SNP2 = sub("^(AX-[0-9]+):", "", pair)
  ) %>%
  select(Trait, Year, SNP1, SNP2, FDR_BH)

# Co:SNP links
co_links <- wald_co %>%
  mutate(SNP_clean = SNP) %>%
  inner_join(
    snp_map %>% select(name, chrom, pos),
    by = c("SNP_clean" = "name")
  ) %>%
  transmute(
    Trait  = Trait,
    Year   = Year,
    chr1   = pcr_chr,
    pos1   = as.numeric(pcr_pos),
    SNP    = SNP,
    chr2   = as.character(chrom),
    pos2   = as.numeric(pos),
    FDR_BH = FDR_BH
  ) %>%
  filter(chr2 %in% VALID_CHRS, is.finite(pos2))

# SNP:SNP links
snps1 <- wald_snpsnp %>%
  mutate(SNP1_clean = SNP1) %>%
  inner_join(
    snp_map %>% select(name, chrom, pos),
    by = c("SNP1_clean" = "name")
  ) %>%
  rename(chr1 = chrom, pos1 = pos)

snps2 <- wald_snpsnp %>%
  mutate(SNP2_clean = SNP2) %>%
  inner_join(
    snp_map %>% select(name, chrom, pos),
    by = c("SNP2_clean" = "name")
  ) %>%
  rename(chr2 = chrom, pos2 = pos)

snpsnp_links <- snps1 %>%
  select(Trait, Year, SNP1, SNP2, FDR_BH, chr1, pos1) %>%
  inner_join(
    snps2 %>% select(Trait, Year, SNP1, SNP2, chr2, pos2),
    by = c("Trait", "Year", "SNP1", "SNP2")
  ) %>%
  filter(
    chr1 %in% VALID_CHRS,
    chr2 %in% VALID_CHRS,
    is.finite(pos1),
    is.finite(pos2)
  )

# ─────────────────────────────────────────────────────────────────────
# 8) Plot circos (no title inside, just legend)
# ─────────────────────────────────────────────────────────────────────

plot_circos <- function() {
  circos.clear()
  op <- par(mai = rep(0.25, 4))
  on.exit(par(op), add = TRUE)
  
  circos.par(
    start.degree = 90,
    gap.after    = rep(2, nrow(genotype_map)),
    track.height = 0.05,
    track.margin = c(0.0008, 0.0008),
    cell.padding = c(0.01, 0.01, 0, 0)
  )
  
  circos.initialize(
    factors = genotype_map$chr,
    xlim    = genotype_map[, c("start", "end")]
  )
  
  # Track 1: Co PCR tick
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.015,
    bg.border = NA,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      if (sec == pcr_chr) {
        circos.segments(
          x0 = pcr_pos, y0 = 0,
          x1 = pcr_pos, y1 = 1,
          col = "red4", lwd = 2.0
        )
      }
    }
  )
  
  # Track 2: chromosome bars + labels
  circos.trackPlotRegion(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.045,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                  col = "grey85", border = NA)
      circos.text(mean(CELL_META$xlim), 0.5, labels = sec, facing = "inside",
                  niceFacing = TRUE, col = "black", cex = 0.75)
    }
  )
  
  # Track 3: SNP density heatmap
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.05,
    bg.border = NA,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      bed_chr <- density_bed[density_bed$chr == sec, , drop = FALSE]
      if (nrow(bed_chr)) {
        circos.rect(
          xleft   = bed_chr$start,
          ybottom = 0,
          xright  = bed_chr$end,
          ytop    = 1,
          col     = col_fun(bed_chr$count),
          border  = NA
        )
      }
    }
  )
  
  # Track 4: Coupel-Ledru ticks
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.04,
    bg.border = NA,
    bg.col = CL_BG_COL,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      df  <- cl_ticks %>% dplyr::filter(chr == sec)
      if (nrow(df)) {
        cols <- ifelse(
          is.finite(df$dist_cap),
          dist_col_fun(df$dist_cap),
          "#636363"
        )
        
        circos.segments(
          x0 = df$pos, y0 = 0,
          x1 = df$pos, y1 = 1,
          col = cols, lwd = CL_LWD
        )
      }
      
      if (CELL_META$sector.numeric.index == 1) {
        circos.text(
          x = mean(CELL_META$xlim), y = 0.5,
          labels = "Coupel-Ledru et al. 2022",
          facing = "bending.inside", niceFacing = TRUE,
          cex = 0.60, col = "#35586B", font = 2
        )
      }
    }
  )
  
  # Track 5: GWAS 2017
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.085,
    track.margin = c(0.0016, 0.0016),
    bg.col = "grey95",
    bg.border = NA,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      df  <- ticks_2017 %>% dplyr::filter(chr == sec)
      if (nrow(df)) {
        for (tr in unique(df$Trait)) {
          dft <- df %>% dplyr::filter(Trait == tr)
          circos.segments(
            x0 = dft$pos, y0 = 0,
            x1 = dft$pos, y1 = 1,
            col = trait_cols_use[[tr]], lwd = 0.85
          )
        }
      }
      
      if (CELL_META$sector.numeric.index == 1) {
        circos.text(
          x = mean(CELL_META$xlim), y = 0.5,
          labels = "2017",
          facing = "bending.inside", niceFacing = TRUE,
          cex = 0.70, col = "black"
        )
      }
    }
  )
  
  # Track 6: GWAS 2018
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.085,
    track.margin = c(0.0016, 0.0016),
    bg.col = "grey95",
    bg.border = NA,
    panel.fun = function(x, y) {
      sec <- CELL_META$sector.index
      df  <- ticks_2018 %>% dplyr::filter(chr == sec)
      if (nrow(df)) {
        for (tr in unique(df$Trait)) {
          dft <- df %>% dplyr::filter(Trait == tr)
          circos.segments(
            x0 = dft$pos, y0 = 0,
            x1 = dft$pos, y1 = 1,
            col = trait_cols_use[[tr]], lwd = 0.85
          )
        }
      }
      
      if (CELL_META$sector.numeric.index == 1) {
        circos.text(
          x = mean(CELL_META$xlim), y = 0.5,
          labels = "2018",
          facing = "bending.inside", niceFacing = TRUE,
          cex = 0.70, col = "black"
        )
      }
    }
  )
  
  # Inner links: Co:SNP
  for (tr in names(trait_cols_use)) {
    links_tr <- co_links %>% filter(Trait == tr)
    if (!nrow(links_tr)) next
    col_tr <- adjustcolor(trait_cols_use[[tr]], alpha.f = 0.85)
    for (i in seq_len(nrow(links_tr))) {
      circos.link(
        sector.index1 = as.character(links_tr$chr1[i]),
        point1        = as.numeric(links_tr$pos1[i]),
        sector.index2 = as.character(links_tr$chr2[i]),
        point2        = as.numeric(links_tr$pos2[i]),
        col           = col_tr,
        lwd           = 2.3
      )
    }
  }
  
  # Inner links: SNP:SNP
  for (tr in names(trait_cols_use)) {
    links_tr <- snpsnp_links %>% filter(Trait == tr)
    if (!nrow(links_tr)) next
    col_tr <- adjustcolor(trait_cols_use[[tr]], alpha.f = 0.40)
    for (i in seq_len(nrow(links_tr))) {
      circos.link(
        sector.index1 = as.character(links_tr$chr1[i]),
        point1        = as.numeric(links_tr$pos1[i]),
        sector.index2 = as.character(links_tr$chr2[i]),
        point2        = as.numeric(links_tr$pos2[i]),
        col           = col_tr,
        lwd           = 0.9
      )
    }
  }
  
  circos.clear()
  
  par(xpd = NA)
  sorted_traits <- sort(names(trait_cols_use))
  legend(
    "topleft",
    legend = sorted_traits,
    col    = as.vector(trait_cols_use[sorted_traits]),
    lwd    = 3,
    bty    = "n",
    cex    = 0.8,
    title  = "Traits"
  )
}

# ─────────────────────────────────────────────────────────────────────
# 9) Save SVG for circos
# ─────────────────────────────────────────────────────────────────────

svglite::svglite(OUT_CIRCOS_SVG, width = 10, height = 10)
plot_circos()
dev.off()

# ─────────────────────────────────────────────────────────────────────
# 10–11) Legend export helpers (SVG only)
# ─────────────────────────────────────────────────────────────────────

save_legend_svg <- function(legend_obj, svg_path, w = 4, h = 6) {
  svglite::svglite(svg_path, width = w, height = h)
  grid.newpage()
  grid.draw(legend_obj)
  dev.off()
}

save_legend_svg(dens_legend,  OUT_DENSITY_LEGEND_SVG,  w = 3.2, h = 4.5)
save_legend_svg(dist_legend,  OUT_DISTANCE_LEGEND_SVG, w = 3.2, h = 4.5)

# ─────────────────────────────────────────────────────────────────────
# 12) Final log
# ─────────────────────────────────────────────────────────────────────

cat("Ticks 2017 (FDR):", nrow(ticks_2017),
    "  Ticks 2018 (FDR):", nrow(ticks_2018), "\n")
cat("CL QTL ticks:", nrow(cl_ticks), "\n")
cat("FDR-significant Co:SNP interactions (Wald):", nrow(wald_co), "\n")
cat("FDR-significant SNP:SNP interactions (Wald):", nrow(wald_snpsnp), "\n")


# ─────────────────────────────────────────────────────────────────────
# 13) Save CSV Results
# ─────────────────────────────────────────────────────────────────────
write.csv(gwas_sig[,-(5:8)],     file.path(BASE_STW, "gwas_results.csv"), row.names = FALSE)
write.csv(co_links,     file.path(BASE_STW, "co_snp.csv"),        row.names = FALSE)
write.csv(snpsnp_links, file.path(BASE_STW, "snp_snp.csv"),       row.names = FALSE)

