# ==============================================================================
# c3_rerun_dyn_annual_53.R
#
# Re-run the dynamic Sentinel-2 scenario at all 53 Blois plots, on the annual
# leaf-area series rebuilt from the unmasked product by
# c3_rebuild_annual_notmasked.R.
#
# The shipped run covers 47 plots and leaves the open archetype P1 with two,
# too few to score, because the masked rasters carry no pixel there. The new
# series covers 53. Every plot is therefore re-run, not only the six that were
# missing: the series changed everywhere, so a mixed set of old and new
# NetCDFs would not be comparable.
#
# Nothing is overwritten. Output goes to a parallel directory and the shipped
# DYN_S2_ANNUAL run is left untouched.
#
# Usage: Rscript c3_rerun_dyn_annual_53.R <chunk> <n_chunks> [<n_plots>]
#   n_plots limits the run to the first plots of the chunk (dry run).
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ncdf4); library(lubridate); library(stringr)
  library(rmusica); library(musica.tools)
})
args    <- commandArgs(trailingOnly = TRUE)
chunk   <- if (length(args) >= 1) as.integer(args[1]) else 1L
K       <- if (length(args) >= 2) as.integer(args[2]) else 1L
n_limit <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_

src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
source("Chapter3_config_CHS41.R")

OUTDIR <- file.path(CFG_C3$out_dir, "nc_dyn_annual_nm", "DYN_S2_ANNUAL")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ts_new <- readRDS("out_files/Chapter3/lai_prep/ts_annual_notmasked.rds")
df     <- as.data.table(readRDS("out_files/Chapter3/lai_prep/df_plots_real53.rds"))
df[, pid := sprintf("X%d_Y%d", as.integer(x), as.integer(y))]
# run_musica_one and the phenology factory index plot_row as a plain data.frame
df     <- as.data.frame(df)
stopifnot(all(df$pid %in% names(ts_new)))

# The scenario is taken from the shipped factory, with the rebuilt series
# injected in place of `annual`. Nothing about the recipe is retyped here, so
# DYN_S2_ANNUAL keeps exactly the LAI source, LAD shape and phenology handling
# it had before; only the series behind the phenology changes.
ts_list <- readRDS("out_files/Chapter3/lai_prep/ts_by_plot.rds")
ts_list$annual <- ts_new
scen <- make_all_scenarios_c3(ts_list, list_year = CFG_C3$list_year,
                              d_opt = CFG_C3$d_opt_m,
                              mode = CFG_C3$scenarios_mode)[["DYN_S2_ANNUAL"]]
stopifnot(!is.null(scen), !is.null(scen$phenology_fn))

# Under abl_flag = "iter" the phenology must carry a row for the day before the
# forcing starts (2020/366), exactly as validate_scenarios_at_hobos() does it.
# Without that seed MuSICA stops immediately with code 16.
seed_pheno <- function(base_ph) {
  if (is.null(base_ph)) return(NULL)
  b <- as.data.frame(base_ph); s1 <- b[1, , drop = FALSE]
  if ("year" %in% names(s1))       s1$year <- 2020L
  if ("Julian_day" %in% names(s1)) s1$Julian_day <- 366L
  rbind(s1, b)
}

idx <- which(((seq_len(nrow(df)) - 1L) %% K) == (chunk - 1L))
if (!is.na(n_limit)) idx <- head(idx, n_limit)
cat(sprintf("chunk %d/%d : %d placettes\n", chunk, K, length(idx)))

t0 <- Sys.time()
for (i in idx) {
  prow <- df[i, ]                       # data.frame row, as the shipped driver passes
  out  <- file.path(OUTDIR, sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(out) && file.size(out) > 1e5) { cat("  ", prow$id_plot, "deja la\n"); next }
  sc_run    <- scen
  ph_seeded <- seed_pheno(scen$phenology_fn(prow))
  if (is.null(ph_seeded)) { cat("  ", prow$id_plot, ": PAS DE PHENOLOGIE\n"); next }
  sc_run$phenology_fn <- function(p) ph_seeded
  tryCatch(
    run_musica_one(plot_row = prow, scenario = sc_run, out_nc_file = out,
                   forcing_file = CFG_C3$forcing_file,
                   musica_cmd = CFG_C3$musica_cmd,
                   extra_setup = list("abl_flag" = sprintf('"%s"', CFG_C3$abl_flag))),
    error = function(e) cat("  ERREUR", prow$id_plot, ":", conditionMessage(e), "\n"))
  ok <- file.exists(out) && file.size(out) > 1e5
  cat(sprintf("  %-8s %s\n", prow$id_plot, if (ok) "ok" else "VIDE"))
}
cat(sprintf("chunk %d termine en %.1f min\n", chunk,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
