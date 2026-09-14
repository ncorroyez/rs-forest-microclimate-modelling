# ==============================================================================
# Variant of scripts/c1_logslope_vs_structure.R answering one question: does the
# before/after comparison survive when each model series is plotted against the
# traits that ACTUALLY parameterised it, instead of the chapter's canonical
# 20 m trait table?
#
# WHY. The two lineages were parameterised on different footprints, identified
# from veget_height_top of the 53 output NetCDFs (residual = MuSICA's 0.5 m
# rounding):
#   before summer -> value of the 10 m pixel under the logger   (max diff 1.11 m)
#   after  summer -> MAX over a 10 m-radius disc of the 10 m raster (0.99 m)
#                    which is what build_hobo_inputs(mode = "buffer25") does
# The canonical table (native 20 m cell) sits 9.46 m from the first and 2.10 m
# from the second, so the standard figure carries a small abscissa mismatch.
#
# WHAT IS AND IS NOT ALIGNED HERE.
#   Hmax and LAI are genuine model inputs, so each series gets its own.
#   The after-summer LAI carries the scan-angle correction (LAI / sec_theta,
#   global factor 1.0513), because that is what the runs were fed.
#   VCI is NOT a model input in either lineage, so "its own traits" is undefined.
#   Its footprint is aligned anyway (cell vs 10 m disc) on the SAME normalized
#   raster, so the panel isolates footprint and not the definition change.
#   The field sensors are parameterised by nothing; they are drawn on the
#   after-summer footprint and labelled as such.
#
# Reads : out_files/Chapter1/tables/tab_logslope_vs_structure.csv (stage from the
#         main script: obs / pre / sim log-slopes per logger)
#         in_files/{lai_z1,max,fCover}_res_10_m.tif, in_files/vci_res_10_m_norm_Blois.tif
#         in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv  (scan-angle factor)
# Writes: out_files/Chapter1/figures/Fig_logslope_vs_structure_ownx.{png,pdf}
#   Rscript scripts/c1_logslope_vs_structure_ownx.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(grid)
  library(terra); library(sf)})
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
source("scripts/_article_style.R")

D <- fread("out_files/Chapter1/tables/tab_logslope_vs_structure.csv")
hp <- st_read(CFG$hobo_geojson, quiet = TRUE); hp <- hp[!hp$id_plot %in% CFG$ids_to_remove, ]
co <- st_coordinates(st_geometry(hp))[, 1:2]
r <- list(LAI = rast("in_files/lai_z1_res_10_m.tif"),
          Hmax = rast("in_files/max_res_10_m.tif"),
          VCI  = rast("in_files/vci_res_10_m_norm_Blois.tif"))
buf <- vect(st_buffer(hp, 10))
cell <- function(k) as.numeric(terra::extract(r[[k]], co)[, 1])
disc <- function(k, f) as.numeric(terra::extract(r[[k]], buf, fun = f, na.rm = TRUE)[, 2])
sec <- fread("in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv")[, mean(sec_theta, na.rm = TRUE)]

# traits as each lineage saw them
T_pre <- data.table(id_plot = hp$id_plot, Hmax = cell("Hmax"), LAI = cell("LAI"),  VCI = cell("VCI"))
T_sim <- data.table(id_plot = hp$id_plot, Hmax = disc("Hmax", "max"),
                    LAI = disc("LAI", "mean") / sec, VCI = disc("VCI", "mean"))
cat(sprintf("scan-angle factor applied to the after-summer LAI: %.4f\n", sec))
for (k in c("Hmax", "LAI", "VCI"))
  cat(sprintf("  %-4s : pre %.2f-%.2f (med %.2f) | after %.2f-%.2f (med %.2f)\n", k,
      min(T_pre[[k]]), max(T_pre[[k]]), median(T_pre[[k]]),
      min(T_sim[[k]]), max(T_sim[[k]]), median(T_sim[[k]])))

MET <- c(Hmax = "Maximum height (m)", VCI = "Vertical complexity index",
         LAI = "Leaf area index (one-sided)")
SRC <- c(obs = "Field sensors (HOBO, n = 53)",
         pre = "MuSICA v3.2.0, before summer",
         sim = "MuSICA v3.2.3 iter, after summer")
TR  <- list(obs = T_sim, pre = T_pre, sim = T_sim)     # obs drawn on the current footprint
L <- rbindlist(lapply(names(MET), function(m) rbindlist(lapply(names(SRC), function(k) {
  x <- merge(D[, .(id_plot, y = get(k))], TR[[k]][, .(id_plot, x = get(m))], by = "id_plot")
  data.table(metric = MET[[m]], x = x$x, y = x$y, src = SRC[[k]]) }))))
L <- L[is.finite(x) & is.finite(y)]
L[, `:=`(metric = factor(metric, levels = unname(MET)), src = factor(src, levels = unname(SRC)))]

PAL <- setNames(c("#128d84", "#8c6bb1", "#E69F00"), SRC)
p <- ggplot(L, aes(x, y, colour = src)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.9) +
  geom_point(alpha = 0.75, size = 1.9) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x", strip.position = "bottom") +
  labs(x = NULL, y = "log(slope)") + theme_article(12) +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.background = element_blank(),
        legend.margin = margin(t = -4, b = 0), legend.box.spacing = unit(4, "pt"))
ggsave_article("out_files/Chapter1/figures/Fig_logslope_vs_structure_ownx", p, 10, 4.2)

st <- L[, {f <- lm(y ~ poly(x, 2)); .(R2 = round(summary(f)$r.squared, 3),
        spearman = round(cor(x, y, method = "spearman"), 3))}, by = .(metric, src)]
cat("\n=== quadratic fit, own-abscissa variant ===\n"); print(st)
cat("DONE\n")
