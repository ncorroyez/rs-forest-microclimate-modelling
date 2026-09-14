# ==============================================================================
# Trait collinearity in the native-20 m design
# Trait correlation structure of the cLHS design — CANONICAL native20 sample.
#
# WHY THIS EXISTS. `scripts/make_trait_collinearity.R` builds the same panel from a
# mixture of `~/Documents/NC_Full` (an absolute path outside the repo), `lad_z05`, and
# `transfer/musica_version_benchmark/...`, i.e. neither the canonical design sample nor
# a self-contained input. The shipped Appendix F figure came from a third source again
# (`pipeline/12_corrplot_traits.R` on the superseded `floor05_v2` sample) and printed
# LAI-fCover 0.81, while the manuscript quotes the native20 values. This script draws
# the panel from the one sample the rest of the chapter uses, so text and figure agree.
#
# Reads : out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds  (400 plots)
# Writes: out_files/Chapter1/figures/fig_trait_collinearity_native20.{png,pdf}
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);source("scripts/_article_style.R")})
# FIVE panels, as Fig. E1's caption promises: one per archetype plus the pooled design.
# Before 2026-07-30 this drew the pooled matrix ONLY, so the caption described a figure
# that did not exist and quoted a P1 value (~0.45) that is not the data (0.79).
source("R/cluster_relabel.R")
S0 <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
V <- intersect(c("LAI","fCover","Hmax","VCI"), names(S0))
S0 <- S0[complete.cases(S0[, ..V])]
S0[, P := factor(relabel_cluster(Cluster), levels = c("P1","P2","P3","P4"))]
S <- na.omit(S0[, ..V])
#' Long-form correlation matrix of the trait set for one panel
#' @param d subset of the design sample (one archetype, or all 400)
#' @param nm panel label
#' @return data.table(a, b, r, panel)
grp <- function(d, nm) { m <- cor(d[, ..V]); z <- as.data.table(as.table(m))
  setnames(z, c("a","b","r")); z[, panel := nm][] }
D <- rbindlist(c(lapply(levels(S0$P), function(p) grp(S0[P == p], p)),
                 list(grp(S0, sprintf("all 400 (n = %d)", nrow(S0))))))
D[, panel := factor(panel, levels = unique(panel))]
D[, `:=`(a=factor(a, levels=V), b=factor(b, levels=rev(V)))]
M <- cor(S)
p <- ggplot(D, aes(a, b, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.1,
            colour = ifelse(abs(D$r) > 0.7, "white", "grey15"), fontface = 2) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, limits = c(-1, 1), name = expression(italic(r))) +
  facet_wrap(~ panel, nrow = 1) +
  coord_equal() + labs(x = NULL, y = NULL) + theme_article(11) +
  theme(axis.text = element_text(face = "bold"), panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave_article("out_files/Chapter1/figures/fig_trait_collinearity_native20", p, 13.5, 3.4)
d <- M[upper.tri(M)]
cat(sprintf("pooled n = %d | pairwise |r| from %.2f to %.2f\n", nrow(S), min(abs(d)), max(abs(d))))
cat("LAI-fCover per archetype: ");
for (p in levels(S0$P)) cat(sprintf("%s %.2f  ", p, cor(S0[P==p]$LAI, S0[P==p]$fCover))); cat("\n")
print(round(M, 3))
cat("DONE\n")
