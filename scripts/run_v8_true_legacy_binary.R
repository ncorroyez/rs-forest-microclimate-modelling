# ==============================================================================
# V8 HOBO : run REF for all 53 HOBO sensors with the TRUE legacy MuSICA binary
# (/home/corroyez/Documents/musica/musica, md5 53072278, Nov 6 2024)
# Inputs : legacy metrics_results_25.csv + allometry_profiles.csv (radius 25 m)
# Expected : reproduce Oct 2025 45/8 buf/amp slope-based regime.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(future); library(furrr)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

LEGACY_BIN <- "/home/corroyez/Documents/musica/musica"
stopifnot(file.exists(LEGACY_BIN))
cli_h1("V8 : true legacy binary ({LEGACY_BIN})")
cli_alert("md5 = {tools::md5sum(LEGACY_BIN)}")

N_WORKERS <- 4
plan(multisession, workers = N_WORKERS)

metr  <- fread(here::here("out_files/radius_test/metrics_results_25.csv"))
allom <- fread(here::here("allometry_profiles.csv"))[radius == 25]

ids_keep <- setdiff(metr$id_plot, CFG$ids_to_remove)
cli_alert("HOBO sensors : {length(ids_keep)}")

hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE)

df_hobo <- list()
for (id in ids_keep) {
  m <- metr[id_plot == id]
  a <- allom[id_plot == id]
  if (nrow(m) == 0 || nrow(a) == 0) next
  lad_cols_raw <- grep("^h_", names(a), value = TRUE)
  z_real <- as.numeric(gsub("h_", "", lad_cols_raw))
  d <- as.numeric(a[1, ..lad_cols_raw]); d[is.na(d)] <- 0
  lad_named <- setNames(d, paste0("LAD_Layer_", z_real))
  cc <- sf::st_coordinates(hobo_pts[hobo_pts$id_plot == id, ])
  row <- as.data.frame(c(
    list(id_plot = id, x = cc[1, "X"], y = cc[1, "Y"],
         LAI = m$pai, Hmax = m$max, fCover = m$clumping),
    as.list(lad_named)))
  df_hobo[[id]] <- row
}
df_hobo <- bind_rows(df_hobo)
df_hobo$fCover[df_hobo$fCover < 0.5] <- 0.5
cli_alert("Built {nrow(df_hobo)} HOBO inputs from legacy CSVs")
cli_alert("LAI mean={round(mean(df_hobo$LAI),2)}, Hmax mean={round(mean(df_hobo$Hmax),2)} m")

# REF scenario : raw LAD from allometry (no rescale)
make_lad_raw_legacy <- function(p) {
  cols <- grep("^LAD_Layer_", names(p), value = TRUE)
  if (length(cols) == 0) return(NULL)
  z   <- as.numeric(sub("^LAD_Layer_", "", cols))
  d   <- as.numeric(p[1, cols, drop = FALSE])
  d[is.na(d)] <- 0
  ok  <- d > 0
  if (!any(ok)) return(NULL)
  data.frame(height = z[ok], density = d[ok])
}

sc_ref <- list(
  name      = "REF_legacybin",
  bit_code  = "1111",
  lai_fn    = function(p) as.numeric(p$LAI),
  hmax_fn   = function(p) as.numeric(p$Hmax),
  fcover_fn = function(p) as.numeric(p$fCover),
  lad_fn    = function(p, hmax = NULL, lai = NULL) make_lad_raw_legacy(p)
)

out_root <- here::here("out_files/musica_hobo_v8_truelegacy")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  out_nc <- file.path(out_root, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
  if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
  jobs[[length(jobs) + 1L]] <- list(pr = pr, out_nc = out_nc, lbl = pr$id_plot)
}
cli_alert("Jobs to run : {length(jobs)}")

t_start <- Sys.time()
future_walk(jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  sc_ref <- list(
    name      = "REF_legacybin",
    bit_code  = "1111",
    lai_fn    = function(p) as.numeric(p$LAI),
    hmax_fn   = function(p) as.numeric(p$Hmax),
    fcover_fn = function(p) as.numeric(p$fCover),
    lad_fn    = function(p, hmax = NULL, lai = NULL) {
      cols <- grep("^LAD_Layer_", names(p), value = TRUE)
      z   <- as.numeric(sub("^LAD_Layer_", "", cols))
      d   <- as.numeric(p[1, cols, drop = FALSE]); d[is.na(d)] <- 0
      ok  <- d > 0
      if (!any(ok)) return(NULL)
      data.frame(height = z[ok], density = d[ok])
    })
  tryCatch(
    run_musica_one(job$pr, sc_ref, job$out_nc,
                   CFG$forcing_file, "/home/corroyez/Documents/musica/musica"),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message)))
}, .options = furrr_options(seed = TRUE))
cli_alert_success("V8 done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min")

# Slope-based regime
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
DT_macro <- as.data.table(df_macro)[, .(date, Tmax_macro)]
regime_rows <- list()
for (f in list.files(out_root, pattern = "\\.nc$", full.names = TRUE)) {
  id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
  res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                        z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) next
  d <- as.data.table(res)
  if (!"Tmax_macro" %in% names(d)) d <- merge(d, DT_macro, by = "date")
  fit <- lm(Tmax_micro ~ Tmax_macro, data = d)
  regime_rows[[id]] <- data.frame(id_plot = id,
                                  slope   = coef(fit)[2],
                                  dT_mean = mean(d$Tmax_micro - d$Tmax_macro, na.rm = TRUE))
}
DT_reg <- as.data.table(do.call(rbind, regime_rows))
cli_h2("V8 (TRUE legacy binary) slope-based regime")
cat(sprintf("  Buffer (slope < 1)        : %d / %d\n",
            sum(DT_reg$slope < 1), nrow(DT_reg)))
cat(sprintf("  Amplification (slope > 1) : %d / %d\n",
            sum(DT_reg$slope > 1), nrow(DT_reg)))
cat(sprintf("  slope mean=%.3f median=%.3f range=[%.3f, %.3f]\n",
            mean(DT_reg$slope), median(DT_reg$slope),
            min(DT_reg$slope), max(DT_reg$slope)))
cat(sprintf("  ΔTmax mean=%+.3f median=%+.3f range=[%+.3f, %+.3f]\n",
            mean(DT_reg$dT_mean), median(DT_reg$dT_mean),
            min(DT_reg$dT_mean), max(DT_reg$dT_mean)))

saveRDS(DT_reg, file.path(out_root, "regime_summary.rds"))
fwrite(DT_reg, file.path(out_root, "regime_summary.csv"))
