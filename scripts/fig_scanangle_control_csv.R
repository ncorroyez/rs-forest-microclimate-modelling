# ==============================================================================
# Control of the scan-angle correction against the lidR-computed profiles: mean LAD
# profile with and without the correction, plus the distribution of the per-plot
# <sec theta> factor.
#
# Works from the already-computed CSVs, so it needs neither the LAS catalogue nor the
# external drive.
#
# Reads : in_files/lad_z05/Blois_lad_z05_r25.csv         (lidR standard profiles)
#         in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv  (corrected profiles + sec_theta)
# Writes: outputs/figures_pipeline/annex/fig_scanangle_control.{png,pdf}
#   Rscript scripts/fig_scanangle_control_csv.R
# ==============================================================================
# Control of the scan-angle correction against the lidR-computed profiles, from the
# already-computed CSVs: Blois_lad_z05_r25 (lidR standard, k=0.5, vertical beams) vs
# Blois_lad_z05_sacorr_r25 (scan-angle corrected, per-plot <sec theta>).
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("scripts/_article_style.R")
DZ <- 0.5
b <- fread("in_files/lad_z05/Blois_lad_z05_r25.csv")
s <- fread("in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv")
lyr <- grep("^LAD_Layer_", names(b), value = TRUE); h <- as.numeric(sub("LAD_Layer_", "", lyr))
#' Mean LAD per height layer across all plots of one profile table
#' @param dt table carrying the LAD_Layer_* columns
#' @return numeric vector, one mean per layer, in the order of `lyr`
mprof <- function(dt) colMeans(as.matrix(dt[, ..lyr]), na.rm = TRUE)
agg <- data.table(z = h, lidr = mprof(b), corr = mprof(s))[order(z)]
LAI_b <- rowSums(as.matrix(b[, ..lyr]), na.rm = TRUE) * DZ
LAI_s <- rowSums(as.matrix(s[, ..lyr]), na.rm = TRUE) * DZ
cat(sprintf("LAI lidR mean %.2f -> corrected %.2f (%+.1f%%, sd of per-plot ratio %.1f%%)\n",
            mean(LAI_b), mean(LAI_s), 100*(mean(LAI_s)/mean(LAI_b)-1), 100*sd(LAI_s/LAI_b)))
cat(sprintf("per-plot <sec theta>: mean %.4f range [%.3f, %.3f]\n", mean(s$sec_theta), min(s$sec_theta), max(s$sec_theta)))

pa <- ggplot(agg[z <= 32], aes(y = z)) +
  geom_path(aes(x = lidr, colour = "lidR (vertical beams)"), linewidth = 1) +
  geom_path(aes(x = corr, colour = "scan-angle corrected"), linewidth = 1, linetype = "22") +
  geom_point(aes(x = lidr, colour = "lidR (vertical beams)"), size = 1) +
  geom_point(aes(x = corr, colour = "scan-angle corrected"), size = 1, shape = 1) +
  scale_colour_manual(values = c("lidR (vertical beams)" = "grey30", "scan-angle corrected" = "#D7191C"), name = NULL) +
  labs(x = expression("Mean LAD  (m"^2*" m"^{-3}*")"), y = "Height (m)",
       subtitle = sprintf("Blois field plots (n=60): LAI %.2f -> %.2f (%+.1f%%)", mean(LAI_b), mean(LAI_s), 100*(mean(LAI_s)/mean(LAI_b)-1))) +
  theme_article(12) + theme(legend.position = "bottom")
pb <- ggplot(s, aes(x = sec_theta)) +
  geom_histogram(bins = 20, fill = "#D7191C", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(x = expression("<sec θ> per plot"), y = "Plots",
       subtitle = "near 1 = beams near-vertical") + theme_article(12)
p <- pa + pb + plot_layout(widths = c(2, 1.1)) + plot_annotation(tag_levels = "a")
OUT <- "outputs/figures_pipeline/annex"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(OUT, "fig_scanangle_control.png"), p, width = 9, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig_scanangle_control.pdf"), p, width = 9, height = 5, device = cairo_pdf)
cat("Saved fig_scanangle_control (png+pdf)\n")
