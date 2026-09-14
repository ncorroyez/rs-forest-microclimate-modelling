# ==============================================================================
# Vertical-profile balance of the REAL canopies: relative height of the LAD
# centroid over the 53 loggers, pooled and per archetype.
#
# WHY THIS EXISTS. Sections 4.1 and Appendix B quoted a centroid range of
# "0.42 to 0.54 of Hmax" with no script behind it. That range is close to the
# interquartile range (0.43-0.52) but was presented as the whole spread, and it
# hid the one real exception: the open P1 sits at a median of 0.24, distinctly
# bottom-heavy. Both statements are now generated from here.
#
# Reads : in_files/lad_z05/Blois_lad_z05_r25.csv                      (LAD profiles)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv       (P labels)
# Writes: out_files/Chapter1/tables/tab_profile_centroid_stats.csv
# ==============================================================================
suppressPackageStartupMessages({library(data.table); source("R/cluster_relabel.R")})
lad <- fread("in_files/lad_z05/Blois_lad_z05_r25.csv")
lyr <- grep("^LAD_Layer_", names(lad), value = TRUE)
hh  <- as.numeric(sub("LAD_Layer_", "", lyr))
Lm  <- as.matrix(lad[, ..lyr]); Lm[!is.finite(Lm)] <- 0
lad[, cen := (Lm %*% hh)[, 1] / pmax(rowSums(Lm), 1e-9)]
lad[, rel := cen / pmax(Hmax, 1e-9)]                       # 0-1, higher = top-heavy
C <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, P)]
D <- merge(lad[, .(id_plot, rel)], C, by = "id_plot")[is.finite(rel)]
S <- rbind(
  D[, .(group = "all loggers", n = .N, median = median(rel),
        q25 = quantile(rel, .25), q75 = quantile(rel, .75),
        min = min(rel), max = max(rel))],
  D[, .(n = .N, median = median(rel), q25 = quantile(rel, .25),
        q75 = quantile(rel, .75), min = min(rel), max = max(rel)), by = .(group = P)][order(group)])
fwrite(S, "out_files/Chapter1/tables/tab_profile_centroid_stats.csv")
print(S[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat(sprintf("\nquoted in the text: median %.2f, middle half %.2f to %.2f; open P1 median %.2f\n",
  S[group == "all loggers"]$median, S[group == "all loggers"]$q25,
  S[group == "all loggers"]$q75, S[group == "P1"]$median))
cat("DONE\n")
