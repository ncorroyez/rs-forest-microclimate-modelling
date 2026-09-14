# ==============================================================================
# Verify the four representative plots of Figure 1b against the design sample.
#
# WHY THIS EXISTS. `out_files/Chapter1/tables/typical_plots_per_cluster.csv` selects one
# plot per archetype for Figure 1b's point-cloud cross-sections, and no script in the repo
# writes it: which four plots the figure shows was unreproducible. This script does not
# re-derive the SELECTION (see the note below), but it makes the file auditable: every
# trait, coordinate and archetype label in it is checked against the canonical cLHS sample,
# so the file can no longer drift from the design it claims to sample.
#
# ON THE SELECTION RULE. The `pid` is the row index in the 400-plot sample
# (sprintf("S%04d", .I)), which reproduces exactly. The `dist` column does not: no
# centroid distance we tried returns these four plots, over z-scored or min-max scaled
# traits, in (LAI, Hmax, fCover), (+VCI), (+FPC1-3) or the FPC space alone. The selection
# very likely carried an extra constraint that is not recorded, most plausibly the
# availability of an extracted point cloud, since Figure 1b needs one and the LAS catalog
# lives on an external drive. Treat `dist` as provenance metadata, not as a reproducible
# quantity, and if the four plots ever need to change, record the rule here.
#
# Reads : out_files/Chapter1/tables/typical_plots_per_cluster.csv
#         out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds
#   Rscript scripts/c1_typical_plots_check.R
# Writes: nothing. It is a check: it prints, and stops non-zero on any mismatch.
#   Rscript scripts/c1_typical_plots_check.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table)})
source("R/cluster_relabel.R")
TYP <- "out_files/Chapter1/tables/typical_plots_per_cluster.csv"
S <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
S[, pid := sprintf("S%04d", .I)][, P := as.character(relabel_cluster(Cluster))]
T <- fread(TYP)

bad <- character()
for (i in seq_len(nrow(T))) {
  r <- T[i]; s <- S[pid == r$pid]
  if (!nrow(s)) { bad <- c(bad, sprintf("%s: pid not in the design sample", r$pid)); next }
  chk <- c(P = identical(r$P, s$P),
           LAI    = isTRUE(all.equal(r$LAI_1s,  s$LAI,    tolerance = 5e-3)),
           Hmax   = isTRUE(all.equal(r$Hmax,    s$Hmax,   tolerance = 0.5)),
           fCover = isTRUE(all.equal(r$fCover,  s$fCover, tolerance = 5e-3)),
           VCI    = isTRUE(all.equal(r$VCI,     s$VCI,    tolerance = 5e-3)),
           x = isTRUE(all.equal(r$x, s$x, tolerance = 1)),
           y = isTRUE(all.equal(r$y, s$y, tolerance = 1)))
  if (!all(chk)) bad <- c(bad, sprintf("%s (%s): mismatch on %s", r$pid, r$P,
                                       paste(names(chk)[!chk], collapse = ", ")))
  cat(sprintf("  %-3s %-6s LAI %.2f  Hmax %5.1f  fCover %.3f  VCI %.2f   %s\n",
              r$P, r$pid, s$LAI, s$Hmax, s$fCover, s$VCI, if (all(chk)) "ok" else "MISMATCH"))
}
if (length(setdiff(c("P1","P2","P3","P4"), T$P)))
  bad <- c(bad, "not all four archetypes are represented")
if (length(bad)) stop("typical_plots_per_cluster.csv does not match the design sample:\n  ",
                      paste(bad, collapse = "\n  "))
cat(sprintf("\nAll %d representative plots match the canonical design sample.\nDONE\n", nrow(T)))
