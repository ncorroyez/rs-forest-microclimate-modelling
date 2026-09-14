# Prime the 2×LAI nc cache: simulate the scenarios needed by the windowed /
# significance / leaf-out analyses that are not yet in nc/. Fullyear sim.
# Run from z_Example root:  Rscript c3_prime_2x.R
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots
all_sc <- make_all_scenarios_c3(prep$ts_by_plot, CFG_C3$list_year, CFG_C3$d_opt_m, CFG_C3$scenarios_mode)
want <- c("STATIC_ALS","DYN_RF","NAIVE_S2_FORMSH","DYN_S2_ATBD","DYN_S2_RESCALED")
sc <- all_sc[intersect(want, names(all_sc))]
cat("priming @2xLAI:", paste(names(sc), collapse=", "), "\n")
ds <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
nc_parent <- file.path(CFG_C3$out_dir, "nc")
val <- validate_scenarios_at_hobos(df, hd, sc, nc_parent, dm, ds,
                                   CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
cat("\nprimed. fullyear metrics:\n")
print(as.data.frame(val$metrics[,c("scenario","n","r2","rmse")]), row.names=FALSE, digits=3)
cat("\nDONE\n")
