# Simulate DYN_RF_PHENO_NODOPT into the common nc/ dir. SAFE sourcing (no h1_/lovb_).
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
source("Chapter3_config.R")
stopifnot(CFG_C3$dopt_variant == "common", CFG_C3$d_opt_m == 7)

prep    <- load_lai_prep(CFG_C3)
ts_list <- prep$ts_by_plot
sc_all  <- make_all_scenarios_c3(ts_list, list_year = CFG_C3$list_year,
                                 d_opt = CFG_C3$d_opt_m, mode = CFG_C3$scenarios_mode)
new <- sc_all["DYN_RF_PHENO_NODOPT"]
stopifnot(!is.null(new$DYN_RF_PHENO_NODOPT))
cat("Simulating:", paste(names(new), collapse=", "), "\n")

ds <- CFG_C3$date_seq
df_macro      <- extract_macro_daily(CFG_C3$forcing_file, ds)
df_hobo_daily <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, df_macro, CFG_C3$ids_to_remove)
validate_scenarios_at_hobos(
  df_hobo_inputs = prep$df_plots, df_hobo_daily = df_hobo_daily,
  scenarios = new, parent_dir = file.path(CFG_C3$out_dir, "nc"),
  df_macro = df_macro, date_seq = ds,
  forcing_file = CFG_C3$forcing_file, musica_cmd = CFG_C3$musica_cmd, force = FALSE)
cat(sprintf("DYN_RF_PHENO_NODOPT nc: %d\nSIM_DONE\n",
            length(list.files(file.path(CFG_C3$out_dir, "nc", "DYN_RF_PHENO_NODOPT"), pattern="\\.nc$"))))
