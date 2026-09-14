suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
stopifnot(CFG_C3$dopt_variant=="common", CFG_C3$d_opt_m==7)
prep <- load_lai_prep(CFG_C3)
sc_all <- make_all_scenarios_c3(prep$ts_by_plot, list_year=CFG_C3$list_year, d_opt=CFG_C3$d_opt_m, mode=CFG_C3$scenarios_mode)
new <- sc_all[c("DYN_RF_CLHS","DYN_RF_CLHS_DOPT")]
stopifnot(all(!sapply(new, is.null)))
cat("Simulating:", paste(names(new), collapse=", "), "\n")
ds <- CFG_C3$date_seq
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
validate_scenarios_at_hobos(prep$df_plots, hd, new, file.path(CFG_C3$out_dir,"nc"), dm, ds,
  CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
for(n in names(new)) cat(sprintf("%-18s %d nc\n", n, length(list.files(file.path(CFG_C3$out_dir,"nc",n),pattern="\\.nc$"))))
cat("SIM_DONE\n")
