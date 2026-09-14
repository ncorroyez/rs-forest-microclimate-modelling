# ==============================================================================
# Table 1 (Section 2.6) — per-archetype trait means and step/SD ratios.
#
# The profile swap is DELIBERATELY NOT a column here: it is a measured effect, and
# Table 1 sits in Methods. It is kept in tab_table1_values.csv as the provenance
# record backing Section 3.2's prose and Table G1.
#
# WHY THIS EXISTS. Table 1 was inline markdown in both manuscripts with no generating
# script, so its values were unprotected: nothing would catch them drifting from the
# data. (The similarly named out_files/Chapter1/tables/tab_cluster_meansd.csv is NOT the
# source — it is the stale old-lineage version, P1 LAI 2.35 against the canonical 1.80.)
#
# Emits BOTH variants, because the two manuscripts differ: the full one carries the
# step/SD column, the short one drops it. Neither carries the profile swap.
#
# Reads : out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds  (design traits)
#         out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv        (profile swap)
#         R/cluster_relabel.R                                               (P1..P4 labels)
# Writes: out_files/Chapter1/tables/tab_table1_full.md
#         out_files/Chapter1/tables/tab_table1_short.md
#         out_files/Chapter1/tables/tab_table1_values.csv
# ==============================================================================
suppressPackageStartupMessages({library(data.table); source("R/cluster_relabel.R")})
LV <- c("P1","P2","P3","P4")
D_LAI <- 0.5; D_HMAX <- 1.0; D_FCOV <- 0.10          # the native perturbation steps

# Traits come from the units table itself, i.e. the 400 plots that were ACTUALLY
# simulated, rather than from re-drawing the cLHS subsample. Re-drawing reproduced the
# design to within a rounding step, but reading it guarantees the table describes the
# canopies the perturbation was run on.
U <- fread("out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv")
U[, P := factor(relabel_cluster(Cluster), levels = LV)]
S <- U[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
SW <- U[, .(swap = median(LAD_dt_all, na.rm = TRUE)), by = P]

T1 <- S[, .(LAI_m = mean(LAI), LAI_s = sd(LAI),
            fc_m  = mean(fCover), fc_s = sd(fCover),
            hm_m  = mean(Hmax), hm_s = sd(Hmax)), by = P][order(P)]
T1 <- merge(T1, SW, by = "P")[order(P)]
# Each fixed native step as a PERCENTAGE of the trait's own archetype mean. Expressed
# per unit, not per standard deviation: the design moved off SD normalization, so an
# SD-denominated column invited exactly the reading it was meant to prevent.
T1[, `:=`(pc_LAI = 100 * D_LAI / LAI_m, pc_fc = 100 * D_FCOV / fc_m, pc_hm = 100 * D_HMAX / hm_m)]
fwrite(T1, "out_files/Chapter1/tables/tab_table1_values.csv")

lab <- c(P1 = "P1 (open)", P2 = "P2", P3 = "P3", P4 = "P4 (dense)")
#' One Markdown row of the full Table 1 (with the step/SD column)
#' @param r one row of the per-archetype summary table
#' @return character, a pipe-delimited Markdown table row
row_full <- function(r) sprintf("| %s | %.2f ± %.2f | %.2f ± %.2f | %.1f ± %.1f | %.0f%% · %.0f%% · %.0f%% |",
  lab[as.character(r$P)], r$LAI_m, r$LAI_s, r$fc_m, r$fc_s, r$hm_m, r$hm_s,
  r$pc_LAI, r$pc_fc, r$pc_hm)
#' One Markdown row of the short Table 1 (no step/SD column)
#' @param r one row of the per-archetype summary table
#' @return character, a pipe-delimited Markdown table row
row_short <- function(r) sprintf("| %s | %.2f ± %.2f | %.2f ± %.2f | %.1f ± %.1f |",
  lab[as.character(r$P)], r$LAI_m, r$LAI_s, r$fc_m, r$fc_s, r$hm_m, r$hm_s)

hdr_full <- c("| Archetype | LAI (mean ± SD) | fCover (mean ± SD) | *H*~max~, m (mean ± SD) | fixed step, % of the archetype mean (LAI · cover · height) |",
              "|:--|:--:|:--:|:--:|:--:|")
hdr_short <- c("| Archetype | LAI (mean ± SD) | fCover (mean ± SD) | *H*~max~, m (mean ± SD) |",
               "|:--|:--:|:--:|:--:|")
mk <- c(hdr_full,  vapply(seq_len(nrow(T1)), function(i) row_full(T1[i]),  ""))
ms <- c(hdr_short, vapply(seq_len(nrow(T1)), function(i) row_short(T1[i]), ""))
mk <- gsub("(?<=[| ])-(?=[0-9])", "\u2212", mk, perl = TRUE)   # house style: Unicode minus
ms <- gsub("(?<=[| ])-(?=[0-9])", "\u2212", ms, perl = TRUE)
writeLines(mk, "out_files/Chapter1/tables/tab_table1_full.md")
writeLines(ms, "out_files/Chapter1/tables/tab_table1_short.md")
cat(paste(mk, collapse = "\n"), "\n\nDONE\n")
