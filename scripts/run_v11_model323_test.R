# ==============================================================================
# Test pipeline run with NEW MuSICA binary (model-3.2.3, May 2025).
#
#   - Binary    : in_files/model-3.2.3/musica
#   - Outputs   : out_files/musica_hobo_v11_model323/1111/  (REF coalition only)
#   - Approach  : single-plot smoke test ; if it succeeds, extend to all 53.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra)
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
cli_h1("V11 test — MuSICA model-3.2.3 (new binary)")
cli_alert("Binary : {NEW_BINARY}")
cli_alert("md5    : {tools::md5sum(NEW_BINARY)}")

OUT_ROOT <- here::here("out_files/musica_hobo_v11_model323")
OUT_REF  <- file.path(OUT_ROOT, "1111")
dir.create(OUT_REF, recursive = TRUE, showWarnings = FALSE)

# ---- Build HOBO inputs (buffer25 mode, legacy convention) -------------------
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                             CFG$ids_to_remove,
                                             hobo_buffer_mode = "buffer25"))
# add id_plot column from geojson order
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>%
  filter(!id_plot %in% CFG$ids_to_remove)
df_hobo$id_plot <- hobo_pts$id_plot
cli_alert("HOBO inputs : {nrow(df_hobo)} plots")

# ---- REF scenario (1111 = ARCH_r_r_r_r) -------------------------------------
fac_scs <- build_factorial_scenarios_archetypes(df_hobo)
sc_ref  <- fac_scs[["ARCH_r_r_r_r"]]
stopifnot(!is.null(sc_ref))

# ---- Single-plot smoke test -------------------------------------------------
# scenario lad_fn / lai_fn / etc. expect data.frame, not data.table
df_hobo_df <- as.data.frame(df_hobo)

test_id <- "41_01"
plot_row <- df_hobo_df[df_hobo_df$id_plot == test_id, , drop = FALSE]
stopifnot(nrow(plot_row) == 1)
out_nc <- file.path(OUT_REF, sprintf("musica_out_HOBO_%s.nc", test_id))

cli_h2("Single-plot smoke test : {test_id}")
cli_alert("Output : {out_nc}")

t0 <- Sys.time()
res <- tryCatch({
  run_musica_one(plot_row, sc_ref, out_nc, CFG$forcing_file, NEW_BINARY)
  list(ok = TRUE, err = NULL)
}, error = function(e) list(ok = FALSE, err = e$message))
dt <- difftime(Sys.time(), t0, units = "secs")

if (res$ok && file.exists(out_nc) && file.size(out_nc) > 1e5) {
  cli_alert_success("Single-plot run SUCCEEDED in {round(as.numeric(dt), 1)} s")
  cli_alert("NC size : {round(file.size(out_nc)/1e6, 2)} MB")
  nc <- nc_open(out_nc)
  vars <- names(nc$var)
  cli_alert("# variables : {length(vars)}")
  if ("Tair_z" %in% vars) {
    raw <- musica.tools::get_variable(nc, "Tair_z")
    cli_alert("Tair_z rows : {nrow(raw)} | unique nair : {length(unique(raw$nair))}")
  } else {
    cli_alert_warning("Tair_z NOT in output — variable name may have changed")
    cli_alert("First 20 vars : {paste(head(vars, 20), collapse=', ')}")
  }
  nc_close(nc)
} else {
  cli_alert_danger("Single-plot run FAILED : {res$err %||% '(no error msg)'}")
  if (file.exists(out_nc)) {
    cli_alert("Partial NC : {file.size(out_nc)} bytes")
  }
  stop("Aborting before fan-out.")
}

# ---- If single plot works, run all 53 --------------------------------------
cli_h2("Extending to all {nrow(df_hobo_df)} HOBO plots")

t_start <- Sys.time()
ok_count <- 0L; fail_count <- 0L
for (i in seq_len(nrow(df_hobo_df))) {
  id <- df_hobo_df$id_plot[i]
  pr <- df_hobo_df[i, , drop = FALSE]
  out_nc <- file.path(OUT_REF, sprintf("musica_out_HOBO_%s.nc", id))
  if (file.exists(out_nc) && file.size(out_nc) > 1e5) { ok_count <- ok_count + 1L; next }
  r <- tryCatch({
    run_musica_one(pr, sc_ref, out_nc, CFG$forcing_file, NEW_BINARY)
    TRUE
  }, error = function(e) { cat(sprintf("[%s] ERROR: %s\n", id, e$message)); FALSE })
  if (isTRUE(r) && file.exists(out_nc) && file.size(out_nc) > 1e5) ok_count <- ok_count + 1L
  else fail_count <- fail_count + 1L
  if (i %% 5 == 0) cli_alert("Progress : {i}/{nrow(df_hobo)} | ok={ok_count} fail={fail_count}")
}
cli_alert_success("Done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min : {ok_count} ok / {fail_count} fail")
