# ==============================================================================
# Chapter 3 — TEST: impact of feeding LAI vs 2×LAI to MuSICA (a few plots).
# STATIC_ALS (parametric phenology) so doubling lai_fn doubles everything
# (LAD integrates to 2×LAI, lai_max=2×LAI, parametric Leaf_area peaks at 2×LAI).
# Compares ΔTmax (sim vs HOBO) summer for LAI vs 2×LAI.
# Run from z_Example root:  Rscript c3_test_2xlai.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots

set.seed(1); df_sub <- df[df$id_plot %in% sample(df$id_plot, 6), ]
cat("test plots:", paste(df_sub$id_plot, collapse=", "), "\n")
cat("their LAI_ALS:", paste(round(df_sub$LAI_ALS,2), collapse=", "), "\n")

.fn <- function(col) function(pr) as.numeric(pr[[col]])
sc <- list(
  STATIC_ALS_1x = list(name="TEST_ALS_1x", lai_fn=.fn("LAI_ALS"), hmax_fn=.fn("Hmax"),
                       fcover_fn=.fn("fCover"), lad_fn=make_lad_real, phenology_fn=NULL),
  STATIC_ALS_2x = list(name="TEST_ALS_2x", lai_fn=function(pr) 2*as.numeric(pr$LAI_ALS),
                       hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=make_lad_real, phenology_fn=NULL))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
# fresh nc dir for the test so it always simulates
test_nc <- file.path(CFG_C3$out_dir, "nc_test2x"); dir.create(test_nc, showWarnings=FALSE, recursive=TRUE)
val <- validate_scenarios_at_hobos(df_sub, hd, sc, test_nc, dm, ds,
                                   CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
cat("\n=== ΔTmax sim vs HOBO, summer (6 plots) ===\n")
print(as.data.frame(val$metrics[,c("scenario","n","r2","rmse","bias")]), row.names=FALSE, digits=3)
cat("\n(bias = mean(sim - obs); positive = MuSICA under-buffers / too warm)\n")
# per-plot bias
d <- val$daily %>% group_by(scenario, id_plot) %>%
  summarise(bias=mean(Delta_sim-Delta_obs, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=scenario, values_from=bias)
cat("\n=== per-plot bias (1x vs 2x) ===\n"); print(as.data.frame(d), row.names=FALSE, digits=3)
cat("\nDONE\n")
