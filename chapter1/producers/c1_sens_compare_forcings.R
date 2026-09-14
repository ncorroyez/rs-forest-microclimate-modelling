# ==============================================================================
# c1_sens_compare_forcings.R
#
# Does the Chapter 1 attribution survive the change of forcing?
#
# The attribution table shipped so far was extracted under in_files/FR-Blo_2021_v2.nc,
# whose boundary-layer height has a median of 49 m. MuSICA runs with an iterative
# atmospheric-boundary-layer coupling, so that depth drives the canopy-atmosphere
# exchange, and 49 m is an order of magnitude below a summer daytime value. The
# CHS41 hybrid carries the MERRA-2 PBLH instead (median 550 m). The runs were also
# wind-corrected, which the CHS41 runs no longer are: neither the radius sweep of
# this chapter nor the sub-canopy validation, which belongs to Chapter 3.
#
# This script compares the two tables on the quantities Chapter 1 actually claims:
# the ranking of the four levers, their magnitude, and whether the leaf-quantity
# lever still strengthens toward the closed canopy.
#
# Inputs : out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv        (FR-Blo, wind-corrected)
#          out_files/Chapter1/tables/sensitivity_perplot_units_v2_CHS41.csv  (CHS41, no wind correction)
# Output : out_files/Chapter1/tables/sens_forcing_comparison.csv
#   Rscript c1_sens_compare_forcings.R
# ==============================================================================

suppressPackageStartupMessages({ library(data.table) })
source("R/cluster_relabel.R")
LV <- c("P1", "P2", "P3", "P4")

A <- fread("out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv")
B <- fread("out_files/Chapter1/tables/sensitivity_perplot_units_v2_CHS41.csv")
for (D in list(A, B)) D[, P := as.character(relabel_cluster(Cluster))]

# The four levers, in the NATIVE steps of Figure 4, which is the only form in which
# they are comparable: the table stores per-unit sensitivities, so each is put back
# on its own step (LAI 0.5, Hmax 1 m, fCover 0.10) and the profile contrast is
# already in degrees. Both sides of a scalar step are averaged in absolute value.
STEP <- list(LAI = 0.5, Hmax = 1.0, fCover = 0.10)
LEV  <- list(LAI    = c("LAI_up_dt_all",  "LAI_dn_dt_all"),
             Hmax   = c("Hmax_up_dt_all", "Hmax_dn_dt_all"),
             fCover = c("fCov_up_dt_all", "fCov_dn_dt_all"),
             profil = c("LAD_dt_all", NA))

amp <- function(D, k) {
  cols <- LEV[[k]]
  if (is.na(cols[2])) return(abs(D[[cols[1]]]))          # contraste reel vs uniforme, en degC
  (abs(D[[cols[1]]]) + abs(D[[cols[2]]])) / 2 * STEP[[k]]
}

out <- rbindlist(lapply(names(LEV), function(k) {
  a <- amp(A, k); b <- amp(B, k)
  rbindlist(lapply(c("tous", LV), function(g) {
    ia <- if (g == "tous") rep(TRUE, nrow(A)) else A$P == g
    ib <- if (g == "tous") rep(TRUE, nrow(B)) else B$P == g
    data.table(levier = k, groupe = g,
               frblo = round(median(a[ia], na.rm = TRUE), 3),
               chs41 = round(median(b[ib], na.rm = TRUE), 3))
  }))
}))
out[, ecart := round(chs41 - frblo, 3)]
fwrite(out, "out_files/Chapter1/tables/sens_forcing_comparison.csv")

cat("\n=== Amplitude mediane du levier (|degC| par pas natif), par archetype ===\n")
print(dcast(out, groupe ~ levier, value.var = c("frblo", "chs41")))

cat("\n=== Classement des leviers (le resultat central du chapitre) ===\n")
for (g in c("tous", LV)) {
  ra <- out[groupe == g][order(-frblo)]$levier
  rb <- out[groupe == g][order(-chs41)]$levier
  cat(sprintf("  %-5s  FR-Blo : %-32s | CHS41 : %-32s | %s\n", g,
              paste(ra, collapse = " > "), paste(rb, collapse = " > "),
              if (identical(ra, rb)) "IDENTIQUE" else "CHANGE"))
}

cat("\n=== Le levier de quantite se renforce-t-il vers la canopee fermee ? ===\n")
for (col in c("frblo", "chs41")) {
  v <- out[levier == "LAI" & groupe %in% LV][[col]]
  cat(sprintf("  %-6s : P1 %.3f -> P4 %.3f  (%s)\n", col, v[1], v[4],
              if (v[4] > v[1]) "oui, se renforce" else "non"))
}
message("\nwrote out_files/Chapter1/tables/sens_forcing_comparison.csv")
