# ==============================================================================
# Chapter 3 — paired bootstrap on the DIFFERENCE between two scenarios' ranking
# correlations, within an archetype. Cited in Methods 2.5 and used in 3.3, 4.1
# and 4.2 to say which product comparisons the samples actually resolve.
# Pure post-processing of junsep_chs41_perplot_macro-chs41.csv.
#   Rscript c3_rank_diff_bootstrap_CHS41.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table)})
source("Chapter3_config_CHS41.R")
TAB <- file.path(CFG_C3$out_dir, "tables", "junsep_chs41_perplot_macro-chs41.csv")
D <- fread(TAB); B <- 4000L; set.seed(1)

# `arch` is NOT called P: inside a data.table i-expression, `P == P` would
# resolve both sides to the column and silently select every row.
rank_diff <- function(arch, metric, a = "S2seul", b = "LiDARfixe") {
  sub <- D[P == arch]
  ids <- sort(unique(sub$id_plot))
  o  <- sub[scenario == b, setNames(get(paste0("obs_", metric, "_junsep")), id_plot)]
  sa <- sub[scenario == a, setNames(get(paste0("sim_", metric, "_junsep")), id_plot)]
  sb <- sub[scenario == b, setNames(get(paste0("sim_", metric, "_junsep")), id_plot)]
  ra <- cor(sa[ids], o[ids]); rb <- cor(sb[ids], o[ids])
  d  <- replicate(B, { k <- sample(ids, length(ids), TRUE)
                       cor(sa[k], o[k]) - cor(sb[k], o[k]) })
  d  <- d[is.finite(d)]; q <- quantile(d, c(.025, .975))
  data.table(P = arch, metric = metric, n = length(ids),
             r_a = round(ra, 3), r_b = round(rb, 3), diff = round(ra - rb, 3),
             lo = round(q[1], 3), hi = round(q[2], 3),
             resolved = (q[1] > 0) | (q[2] < 0))
}
R <- rbindlist(lapply(c("P1","P2","P3","P4"), function(p)
       rbindlist(lapply(c("slope","dtmax"), function(m) rank_diff(p, m)))))
print(R)
fwrite(R, file.path(CFG_C3$out_dir, "tables", "k_rankdiff_S2_vs_LiDAR_macro-chs41.csv"))
cat("\nDONE\n")
