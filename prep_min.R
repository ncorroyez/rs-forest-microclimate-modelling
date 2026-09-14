# Minimal, SAFE LAI prep: sources only library files (excludes h1_*/lovb_* which
# have top-level executable code — lovb_run_extensions.R can launch ~80 min of
# MuSICA on source). Rebuilds df_plots_lai.rds + ts_by_plot.rds only.
suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(mgcv); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]      # drop auto-executing analysis scripts
invisible(lapply(src, source))
source("Chapter3_config.R")
stopifnot(CFG_C3$dopt_variant == "common", CFG_C3$d_opt_m == 7)

cat("-- Loading LiDAR rasters + HOBO plots\n")
r_lidar  <- load_lidar_rasters(CFG_C3$in_dir, agg_factor = CFG_C3$agg_factor)
df_plots <- build_hobo_inputs(CFG_C3$hobo_geojson, r_lidar$stack, CFG_C3$ids_to_remove)
df_plots$plot_id <- sprintf("X%d_Y%d", round(df_plots$x), round(df_plots$y))
cat(sprintf("   %d HOBO plots\n", nrow(df_plots)))

cat("-- Building LAI correction table\n")
df_plots <- build_lai_correction_table(
  df_plots = df_plots, paths = make_lai_paths(CFG_C3),
  in_dir = CFG_C3$in_dir, agg_factor = CFG_C3$agg_factor)
saveRDS(df_plots, file.path(CFG_C3$lai_prep_dir, "df_plots_lai.rds"))

cat("-- Building S2 time series (incl. rf_pheno)\n")
ts_list <- build_all_s2_ts(df_plots, CFG_C3)
saveRDS(ts_list, file.path(CFG_C3$lai_prep_dir, "ts_by_plot.rds"))
n_per <- sapply(ts_list, function(x) if (is.null(x)) 0L else length(x))
cat(sprintf("   ts variants: %s\n", paste(names(n_per), n_per, sep="=", collapse=" | ")))
rp <- ts_list$rf_pheno
if (!is.null(rp) && length(rp) > 0) {
  pk <- sapply(rp, function(d) max(d$lai))
  cat(sprintf("   rf_pheno: %d plots, summer-peak mean=%.2f (ATBD ref=%.2f)\n",
              length(rp), mean(pk), mean(df_plots$LAI_S2_ATBD, na.rm = TRUE)))
} else cat("   rf_pheno: EMPTY\n")
cat("PREP_MIN_DONE\n")
