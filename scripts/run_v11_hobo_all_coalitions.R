# ==============================================================================
# V11 (model-3.2.3) — Run all 16 HOBO coalitions × 53 plots = 848 sims.
# Parallel, 4 workers. Output : out_files/musica_hobo_v11_model323/{bit}/
# REF (1111) already exists, will be skipped.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(parallel)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

NEW_BINARY <- here::here("in_files/model-3.2.3/musica")
stopifnot(file.exists(NEW_BINARY))
cli_h1("V11 HOBO — 16 coalitions × 53 plots (mclapply / fork)")
cli_alert("Binary : {NEW_BINARY}")

N_WORKERS <- 4

OUT_ROOT <- here::here("out_files/musica_hobo_v11_model323")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                             CFG$ids_to_remove,
                                             hobo_buffer_mode = "buffer25"))
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
df_hobo$id_plot <- hobo_pts$id_plot
df_hobo_df <- as.data.frame(df_hobo)
cli_alert("HOBO plots : {nrow(df_hobo_df)}")

fac_scs <- build_factorial_scenarios_archetypes(df_hobo)
cli_alert("Coalitions : {length(fac_scs)}")

# Build job list across all 16 coalitions × 53 plots
jobs <- list()
for (sc_nm in names(fac_scs)) {
  bit <- .ARCH_COAL_MAP[sc_nm]
  out_dir <- file.path(OUT_ROOT, bit)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(df_hobo_df))) {
    id <- df_hobo_df$id_plot[i]
    pr <- df_hobo_df[i, , drop = FALSE]
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", id))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs) + 1L]] <- list(pr = pr, sc = fac_scs[[sc_nm]],
                                       out_nc = out_nc,
                                       lbl = paste(bit, id))
  }
}
cli_alert("Jobs to run : {length(jobs)}")
if (length(jobs) == 0) {
  cli_alert_success("Nothing to do — all 848 sims already present.")
  quit(save = "no")
}

t0 <- Sys.time()
mclapply(jobs, function(job) {
  tryCatch(
    run_musica_one(job$pr, job$sc, job$out_nc, CFG$forcing_file, NEW_BINARY),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message))
  )
}, mc.cores = N_WORKERS, mc.preschedule = FALSE)
dt <- difftime(Sys.time(), t0, units = "mins")
cli_alert_success("HOBO coalitions done in {round(as.numeric(dt), 1)} min")
