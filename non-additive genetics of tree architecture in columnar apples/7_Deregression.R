# ──────────────────────────────────────────────────────────────────────────────
# Deregression from saved inputs → reliability, weights, dGEBVs (SE-based)
# ──────────────────────────────────────────────────────────────────────────────
rm(list=ls(all=TRUE)); graphics.off(); closeAllConnections()

DATA_DIR <- "C:/Users/nguevenc/Desktop/R_Working_Directory/Lidar/Single_trait_walkthrough"
IN_DIR   <- file.path(DATA_DIR, "h2_and_deregression_input")
OUT_DIR  <- file.path(DATA_DIR, "deregression_output")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TRAIT <- "in_len"   # in_len, lsh_ang_1, tr_brz_prop, tr_coni, tr_len
YEAR  <- 2017       # 2017, 2018

file_blup <- file.path(IN_DIR, sprintf("BLUPs_SEs_h2_%s_%d_pop.csv", TRAIT, YEAR))
blup <- read.csv(file_blup, check.names = FALSE)[, c("Genotype","predicted.value","std.error","AVSED","h2_Cullis", "pop","pop_BLUE","co_loc", "Gvar")]
names(blup) <- c("Genotype","GEBV","SE","AVSED","h2_Cullis","pop","pop_BLUE","co_loc", "Gvar")
Gvar  <- as.numeric(blup$Gvar)

rel    <- pmax(0, pmin(1, 1 - (blup$SE^2) / Gvar))
r      <- sqrt(rel)
dGEBV  <- ifelse(r > 0, blup$GEBV / r, NA_real_)
weight <- ifelse(rel < 1, rel / (1 - rel), NA_real_)

out <- data.frame(
  Trait      = TRAIT,
  Year       = YEAR,
  Genotype   = blup$Genotype,
  dGEBV      = dGEBV,
  weight     = weight,
  rel        = rel,
  GEBV       = blup$GEBV,
  SE         = blup$SE,
  AVSED      = blup$AVSED,
  h2_Cullis  = blup$h2_Cullis,
  pop        = blup$pop,
  pop_BLUE   = blup$pop_BLUE,
  co_loc = blup$co_loc,
  check.names = FALSE
)

outfile <- file.path(OUT_DIR, sprintf("dGEBV_%s_%d_pop.csv", TRAIT, YEAR))
write.csv(out, outfile, row.names = FALSE)
cat("✔ Wrote deregression output to:\n", outfile, "\n", sep = "")

cat(sprintf("N genotypes: %d | rel median = %.3f | h2 = %.3f | NA dGEBV: %d\n",
            nrow(out), median(out$rel, na.rm=TRUE),
            mean(out$h2_Cullis, na.rm=TRUE),
            sum(is.na(out$dGEBV))))
