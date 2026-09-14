# ==============================================================================
# Chapter 3 — DRY RUN: targeted re-simulation of the empty (failed) DYN sims
# under v3.2.3-iter. Tests ONE scenario (DYN_S2_ATBD): move its empty nc aside,
# re-run only the missing plots, verify they come back valid (8760 steps), time it.
# The runner skips per-plot when the nc file exists, so moving the empties away is
# what makes it re-simulate exactly those 8.
# Run from z_Example root:  Rscript dry_run_resim_v323.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr)
  library(tidyr); library(stringr); library(rmusica); library(musica.tools)
})
invisible(lapply(grep("/(h1_|lovb_)", list.files("R", pattern="\\.R$", full.names=TRUE),
                      value=TRUE, invert=TRUE), source))
source("Chapter3_config.R")

TARGET <- "DYN_S2_ATBD"
ncp    <- file.path(CFG_C3$out_dir, "nc")
bak    <- file.path(CFG_C3$out_dir, "nc_empty_bak", TARGET)
dir.create(bak, recursive=TRUE, showWarnings=FALSE)

# 1) move empty/failed nc (< 100 KB) aside so the runner re-creates them
scen_dir <- file.path(ncp, TARGET)
allnc <- list.files(scen_dir, pattern="\\.nc$", full.names=TRUE)
empty <- allnc[file.info(allnc)$size < 100*1024]
cat(sprintf("Scenario %s: %d nc total, %d empty (<100KB) to re-sim:\n", TARGET, length(allnc), length(empty)))
cat("  ", paste(str_extract(basename(empty), "(?<=HOBO_).*(?=\\.nc)"), collapse=" "), "\n")
file.rename(empty, file.path(bak, basename(empty)))

# 2) setup (mirror Chapter3_main.R)
prep    <- load_lai_prep(CFG_C3); df_hobo <- prep$df_plots; ts_list <- prep$ts_by_plot
scenarios_c3 <- make_all_scenarios_c3(ts_list, list_year=CFG_C3$list_year,
                                      d_opt=CFG_C3$d_opt_m, mode=CFG_C3$scenarios_mode)
df_macro      <- extract_macro_daily(CFG_C3$forcing_file, CFG_C3$date_seq)
df_hobo_daily <- read_hobo_daily(CFG_C3$hobo_temp_csv, CFG_C3$date_seq, df_macro, CFG_C3$ids_to_remove)

# 3) re-run ONLY the target scenario, all 53 plots (present 45 skip, moved 8 re-sim)
t0 <- Sys.time()
val <- validate_scenarios_at_hobos(
  df_hobo_inputs = df_hobo, df_hobo_daily = df_hobo_daily,
  scenarios = scenarios_c3[TARGET], parent_dir = ncp,
  df_macro = df_macro, date_seq = CFG_C3$date_seq,
  forcing_file = CFG_C3$forcing_file, musica_cmd = CFG_C3$musica_cmd,
  force = FALSE, abl_flag = CFG_C3$abl_flag)
dt <- as.numeric(difftime(Sys.time(), t0, units="secs"))

# 4) verify the re-created nc are now valid (non-empty time dim)
cat(sprintf("\n=== VERIFY (elapsed %.0f s for %d re-sims, ~%.1f s/plot) ===\n",
            dt, length(empty), dt/max(1,length(empty))))
ok <- 0
for(f in file.path(scen_dir, basename(empty))){
  id <- str_extract(basename(f), "(?<=HOBO_).*(?=\\.nc)")
  nt <- tryCatch({ nc<-nc_open(f); n<-nc$dim$time$len; nc_close(nc); n }, error=function(e) -1)
  cat(sprintf("  %-8s time steps = %s\n", id, nt)); if(is.numeric(nt) && nt>1000) ok <- ok+1
}
cat(sprintf("\nRESULT: %d/%d re-sims valid (>1000 steps). %s\n", ok, length(empty),
            if(ok==length(empty)) "DRY RUN OK -> safe to scale to all DYN scenarios" else "SOME FAILED -> investigate before scaling"))
