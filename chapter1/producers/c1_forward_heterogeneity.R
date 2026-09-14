# ==============================================================================
# Step-by-step validation vs TERRAIN / HETEROGENEITY covariates. As LiDAR
# variables are added cumulatively (pooled importance order LAI > LAD > fCover >
# Hmax), we track the correlation of the per-plot validation residual
# (observed - simulated mean ΔTmax, 53 HOBO) with each out-of-model covariate:
# elevation, slope, TWI, northness, canopy-height SD, rumple, gap fraction.
# Idea: once the LiDAR-structural signal is captured (mainly by LAI), the residual
# should be increasingly governed by what a 1-D plot model cannot see (terrain +
# within-pixel heterogeneity, esp. gap fraction). Uses cached coal_metrics + the
# Not_Masked covariates (no recompute).
#   PIPE_BRANCH=z05 Rscript c1_forward_heterogeneity.R
# Out: out_files/Chapter1/figures/FigSh_forward_heterogeneity.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
if (!exists("PIPE")) { Sys.setenv(PIPE_BRANCH = "z05"); source("pipeline/00_config.R") }

coal <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "coal_metrics.rds")))   # id_plot, bit, Tmax_all
dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, CFG$date_seq)); mean_macro <- mean(dm$Tmax_macro, na.rm = TRUE)
V    <- readRDS(file.path(PIPE$OUT_DATA, "ref_validation.rds"))
obs  <- as.data.table(V$obs_daily)[, .(obs = mean(Delta_obs, na.rm = TRUE)), by = id_plot]
cov  <- fread("outputs/figures_pipeline/annex/tab_residual_vs_topo.csv")
COVS <- c("elevation","slope","twi","northness","sd_height","rumple","gap_fraction")

# cumulative coalitions in pooled importance order LAI > LAD > fCover > Hmax
steps <- data.table(step = 0:4,
                    bit  = c("0000","1000","1001","1011","1111"),
                    add  = c("Baseline","+LAI","+LAD","+fCover","+Hmax (REF)"))

res <- rbindlist(lapply(seq_len(nrow(steps)), function(i) {
  b <- steps$bit[i]
  d <- merge(coal[bit == b, .(id_plot, sim = Tmax_all - mean_macro)], obs, by = "id_plot")
  d <- merge(d, cov[, c("id_plot", COVS), with = FALSE], by = "id_plot")
  d[, residual := obs - sim]
  rbindlist(lapply(COVS, function(c) {
    ok <- is.finite(d$residual) & is.finite(d[[c]])
    data.table(step = steps$step[i], add = steps$add[i], covariate = c,
               r = if (sum(ok) > 2) cor(d$residual[ok], d[[c]][ok]) else NA_real_, n = sum(ok))
  }))
}))
res[, covariate := factor(covariate, levels = COVS,
      labels = c("elevation","slope","TWI","northness","height SD","rumple","gap fraction"))]
fwrite(res, file.path(PIPE$OUT_TAB, "tab_forward_heterogeneity.csv"))
cat("=== r(residual, covariate) per step ===\n"); print(dcast(res, covariate ~ step, value.var = "r")[, lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)])

pal <- c("gap fraction"="#D7191C","height SD"="#2C7BB6","rumple"="#1A9850","elevation"="#7B3294",
         "slope"="#FDAE61","TWI"="#66C2A5","northness"="grey55")
lab4 <- res[step == 4]
p <- ggplot(res, aes(step, r, colour = covariate, group = covariate)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_text_repel(data = lab4, aes(label = covariate), size = 3, fontface = "bold", hjust = 0,
                  direction = "y", nudge_x = 0.15, segment.size = 0.2, show.legend = FALSE, max.overlaps = Inf) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(breaks = 0:4, labels = c("Baseline","+LAI","+LAD","+fCover","+Hmax (REF)"),
                     expand = expansion(mult = c(0.04, 0.22))) +
  labs(x = "Variables added (cumulative, pooled importance order)",
       y = "r (validation residual vs covariate, 53 plots)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave("out_files/Chapter1/figures/FigSh_forward_heterogeneity.png", p, width = 9, height = 5.6, dpi = 200, bg = "white")
cat("DONE -> out_files/Chapter1/figures/FigSh_forward_heterogeneity.png\n")
