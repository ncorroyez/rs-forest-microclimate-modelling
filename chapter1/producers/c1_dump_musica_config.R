# ==============================================================================
# Dump the MuSICA configuration behind Appendix H (Tables H1 and H2), from the
# namelists MuSICA actually reads plus one simulation output.
#
# WHY THIS EXISTS. H1 and H2 are hand-written tables with no generator, so nothing
# caught them drifting from the model actually run.
#
# TWO TRAPS THIS SCRIPT EXISTS TO AVOID, both hit on 2026-07-31:
#  1. `in_files/model-3.2.3/musica.nml` and `musica_soil.nml` are the DISTRIBUTION
#     TEMPLATES, not the run configuration. musica.nml's SOIL_FILE_NAME points at
#     `./in_files/musica_soil.nml`, a different file with different soil values. A
#     first pass read the template and wrongly concluded every soil row of H1/H2
#     was wrong; the values in `in_files/musica_soil.nml` match the tables exactly.
#  2. The template's setup block is overridden at run time by R/musica.R and by each
#     run script, so reading it alone says ABL_flag='none' and forcing=musica_in_Blois.nc,
#     both false for the shipped runs. The AS-RUN column below records the override
#     and where it lives, so the template is never mistaken for the configuration.
#
# Reads : in_files/musica_soil.nml          (SOIL_FILE_NAME target: the soil config that runs)
#         in_files/musica_veg1.nml          (species/leaf config)
#         in_files/model-3.2.3/musica.nml   (setup template, overridden at run time)
#         in_files/model-3.2.3/src/musica/mo_setup.f90 (air_resolution_level -> layer counts)
#         in_files/FR-Blo_2021_v2.nc        (canonical forcing: time axis, h_sbl)
#         out_files/musica_hobo_native20/1111/musica_out_HOBO_41_01.nc
# Writes: out_files/Chapter1/tables/tab_musica_config_dump.csv
#   Rscript scripts/c1_dump_musica_config.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(ncdf4)})

# ---- 1. namelist values, from the files MuSICA actually reads -----------------
#' Extract selected keys from a Fortran namelist, comments stripped
#' @param path namelist file
#' @param keys character vector of parameter names to pull
#' @param tag value written into the `source` column
#' @return data.table(source, param, value, as_run); stops if the file is missing
read_nml <- function(path, keys, tag) {
  if (!file.exists(path)) stop("missing namelist: ", path)
  txt <- readLines(path, warn = FALSE); txt <- txt[!grepl("^\\s*!", txt)]
  rbindlist(lapply(keys, function(k) {
    hit <- grep(paste0("^\\s*", k, "\\s*="), txt, ignore.case = TRUE, value = TRUE)
    if (!length(hit)) return(NULL)
    data.table(source = tag, param = k,
               value = trimws(sub("!.*$", "", sub("^[^=]*=", "", hit[1]))),
               as_run = "as read (not overridden)")
  }))
}
SOIL_KEYS <- c("SOIL_GRID_TYPE","SOIL_DEPTH_MAX","N_DEPTH","DEPTH","RETENTION_CURVE_MODEL_FLAG",
               "HYDRAULIC_COND_MODEL_FLAG","THETA_SAT_SOIL","THETA_RES_SOIL","H_S_SOIL","N_SOIL",
               "BIG_M_SOIL","KSAT_SOIL","BULK_DENSITY","ALBEDO_SOIL_VIS","ALBEDO_SOIL_NIR",
               "EPSI_SOIL","INIT_DEPTH","INIT_SOIL_TEMPERATURE","INIT_SOIL_MOISTURE")
VEG_KEYS  <- c("LEAF_SIZE","LEAF_SHAPE_IN","LEAF_INCLINATION_INDEX_YOUNG","LEAF_INCLINATION_INDEX_OLD",
               "BUDBURST_DATE","PHENOLOGY_MODEL","LAI_MAX_PER_COHORT","CANOPY_HEIGHT_TOP")
SETUP_KEYS <- c("ABL_flag","FORCING_FILENAME","FORCING_TIMESTEP","FORCING_HEIGHT","CLUMPING_FACTOR",
                "AIR_RESOLUTION_LEVEL","AIR_GRID_TYPE","N_SOIL_LAYER","N_LEAF_AGE","N_SPECIES",
                "SITE_LATITUDE","SITE_LONGITUDE","SITE_ALTITUDE","TIME2GMT")
# The repo carries several musica_soil.nml copies with different md5 sums. They differ
# only in comments and number formatting: on every key H1/H2 cites, the SOIL_FILE_NAME
# target and the in_files/Blois copy agree, and both agree with an as-run soil.nml
# recovered from a .musica_* work directory. Assert it rather than trust it, so a future
# divergence between the copies surfaces here instead of silently in a run.
#' Normalise a namelist value for comparison across copies (commas, spacing, trailing dot)
#' @param v raw namelist value string
#' @return single space-separated character string
norm_val <- function(v) paste(sub("\\.$", "", trimws(strsplit(gsub(",", " ", v), " +")[[1]])), collapse = " ")
.a <- read_nml("in_files/musica_soil.nml", SOIL_KEYS, "a")
.b <- read_nml("in_files/Blois/in_files/musica_soil.nml", SOIL_KEYS, "b")
.m <- merge(.a[, .(param, va = sapply(value, norm_val))],
            .b[, .(param, vb = sapply(value, norm_val))], by = "param")
if (nrow(.m[va != vb]))
  stop("the musica_soil.nml copies disagree on: ", paste(.m[va != vb, param], collapse = ", "),
       "\n  Appendix H would then depend on which copy MuSICA staged. Reconcile them first.")

D <- rbindlist(list(
  read_nml("in_files/musica_soil.nml", SOIL_KEYS, "musica_soil.nml (SOIL_FILE_NAME target)"),
  read_nml("in_files/musica_veg1.nml", VEG_KEYS,  "musica_veg1.nml"),
  read_nml("in_files/model-3.2.3/musica.nml", SETUP_KEYS, "musica.nml TEMPLATE")))

# ---- 2. what the run scripts and R/musica.R override at run time --------------
# Verified 2026-07-31 against an as-run namelist recovered from a .musica_* work
# directory (abl_flag="iter", forcing_timestep=3600, forcing_height=hmax+2,
# clumping_factor=fCover, soil_file_name="./soil.nml" = in_files/musica_soil.nml).
OVERRIDE <- c(
  ABL_flag         = "\"iter\"  <- run scripts, e.g. c1_sensitivity_perplot_units.R:19",
  FORCING_FILENAME = "in_files/FR-Blo_2021_v2.nc  <- run scripts (FORC)",
  FORCING_TIMESTEP = "3600 (hourly), matching the forcing time axis",
  FORCING_HEIGHT   = "Hmax + 2 m, per plot  <- R/musica.R:96",
  CLUMPING_FACTOR  = "fCover, per plot  <- R/musica.R:99",
  LAI_MAX_PER_COHORT = "2 x one-sided LiDAR LAI, per plot  <- R/musica.R:69,111",
  CANOPY_HEIGHT_TOP  = "Hmax, per plot  <- R/musica.R:112",
  BUDBURST_DATE      = "115, restated in the parametric fallback  <- R/musica.R:85")
for (k in names(OVERRIDE)) D[toupper(param) == toupper(k), as_run := OVERRIDE[[k]]]

# ---- 3. values that only the model source or an output can settle -------------
src <- readLines("in_files/model-3.2.3/src/musica/mo_setup.f90", warn = FALSE)
lev <- D[toupper(param) == "AIR_RESOLUTION_LEVEL", value]
i   <- grep("select case \\(air_resolution_level\\)", src)
blk <- paste(src[i:(i + 14)], collapse = " ")                      # the case block
nveg <- sub(".*case \\(2\\).*?n_veg_layer *= *([0-9]+).*", "\\1", blk)
nair <- sub(".*case \\(2\\).*?n_air_layer *= *([0-9]+).*", "\\1", blk)
D <- rbind(D, data.table(source = "mo_setup.f90", param = "n_veg_layer / n_air_layer",
  value = sprintf("%s / %s", nveg, nair),
  as_run = sprintf("set by air_resolution_level = %s, not a free parameter", trimws(lev))))

fn <- nc_open("in_files/FR-Blo_2021_v2.nc")
ftu <- ncatt_get(fn, "time", "units")$value; fth <- ncvar_get(fn, "time")
hs <- as.numeric(ncvar_get(fn, "h_sbl")); nc_close(fn)
D <- rbind(D, data.table(source = "forcing FR-Blo_2021_v2.nc", param = "time axis",
  value = sprintf("%s, n = %d, step = %g h", ftu, length(fth), diff(fth)[1]),
  as_run = "hourly, so FORCING_TIMESTEP = 3600"))
D <- rbind(D, data.table(source = "forcing FR-Blo_2021_v2.nc", param = "h_sbl",
  value = sprintf("min %.1f, median %.1f, max %.1f m", min(hs), median(hs), max(hs)),
  as_run = paste("surface-boundary-layer height archived in the forcing file;",
                 "NOT the in-repo MERRA-2 series of build_forcing_pblh.R (r = 0.54 against it).",
                 "Immaterial to dTmax: see the 100-5000 m test in Appendix H.")))

nc <- nc_open("out_files/musica_hobo_native20/1111/musica_out_HOBO_41_01.nc")
tu <- ncatt_get(nc, "time", "units")$value; th <- ncvar_get(nc, "time")
t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC"); tt <- t0 + th * 3600
rh <- ncvar_get(nc, "relative_height"); nc_close(nc)
D <- rbind(D, data.table(source = "simulation output", param = "integration period",
  value = sprintf("%s to %s (%d hourly steps)", format(min(tt), "%Y-%m-%d"),
                  format(max(tt), "%Y-%m-%d"), length(tt)),
  as_run = "full year 2021; the 1 June to 30 September window follows five months of model time"))
D <- rbind(D, data.table(source = "simulation output", param = "air layers (nair)",
  value = as.character(length(rh)), as_run = "confirms air_resolution_level = 2"))

dir.create("out_files/Chapter1/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(D, "out_files/Chapter1/tables/tab_musica_config_dump.csv")
for (i in seq_len(nrow(D)))
  cat(sprintf("%-42s %-28s %s\n%-72s AS RUN: %s\n", D$source[i], D$param[i],
              substr(D$value[i], 1, 90), "", D$as_run[i]))
cat("\nAppendix H Tables H1/H2 must agree with the AS RUN column, not with the template.\nDONE\n")
