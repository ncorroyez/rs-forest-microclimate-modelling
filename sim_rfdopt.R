# Simulate the two new RF-DOPT scenarios into the common nc/ dir (no clobber:
# new scenario dirs only). Config must be common/7.
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))
source("Chapter3_config.R")
stopifnot(CFG_C3$dopt_variant == "common", CFG_C3$d_opt_m == 7)

prep    <- load_lai_prep(CFG_C3)
ts_list <- prep$ts_by_plot
sc_all  <- make_all_scenarios_c3(ts_list, list_year = CFG_C3$list_year,
                                 d_opt = CFG_C3$d_opt_m, mode = CFG_C3$scenarios_mode)
new <- sc_all[c("STATIC_RF_DOPT", "DYN_RF_DOPT")]
cat("New scenarios to simulate:", paste(names(new), collapse = ", "), "\n")

ds <- CFG_C3$date_seq
df_macro      <- extract_macro_daily(CFG_C3$forcing_file, ds)
df_hobo_daily <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, df_macro, CFG_C3$ids_to_remove)
val <- validate_scenarios_at_hobos(
  df_hobo_inputs = prep$df_plots,
  df_hobo_daily  = df_hobo_daily,
  scenarios      = new,
  parent_dir     = file.path(CFG_C3$out_dir, "nc"),
  df_macro       = df_macro, date_seq = ds,
  forcing_file   = CFG_C3$forcing_file, musica_cmd = CFG_C3$musica_cmd, force = FALSE
)
cat("Done. nc counts:\n")
for (n in names(new))
  cat(sprintf("  %-16s %d nc\n", n,
              length(list.files(file.path(CFG_C3$out_dir, "nc", n), pattern = "\\.nc$"))))
