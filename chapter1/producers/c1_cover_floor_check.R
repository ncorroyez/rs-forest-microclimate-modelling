# ==============================================================================
# Cover-floor control (Appendix C): does the 0.5 fCover floor explain the model's
# failure in open canopy? Re-scores the no-floor re-run against the floored
# baseline and the loggers, on the ALIGNED clock (R/dtmax_convention.R).
#
# WHY THIS EXISTS. The claim was in the manuscript with NO script behind it, and
# its observed value (+3.8 degC) was the PRE-alignment observation: the aligned
# figure is +2.73. The text also said "the four floored open plots" when seven
# loggers sit below the floor and only four were re-run. Both are fixed, and the
# experiment is now reproducible.
#
# Reads : out_files/musica_hobo_nofloor/1111  (4 re-run loggers, true cover)
#         out_files/musica_hobo_native20/1111 (floored baseline)
#         out_files/Chapter1/tables/tab_hobo_native20_validation.csv (aligned obs)
# Writes: out_files/Chapter1/tables/tab_cover_floor_check.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(terra);library(sf)})
src <- list.files("R","\\.R$",full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src, source))
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day"); Z <- 1
MREF <- macro_ref("in_files/FR-Blo_2021_v2.nc", ds)
NF <- "out_files/musica_hobo_nofloor/1111"; FL <- "out_files/musica_hobo_native20/1111"
ids <- sub("musica_out_HOBO_(.*)\\.nc", "\\1", list.files(NF, "\\.nc$"))
V  <- fread("out_files/Chapter1/tables/tab_hobo_native20_validation.csv")
g  <- st_read("in_files/data_Blois_utm31n.geojson", quiet=TRUE)
r0 <- rast("in_files_native20/fCover_res_10_m.tif"); g <- st_transform(g, crs(r0))
FC <- data.table(id_plot = g$id_plot, fcover_true = terra::extract(r0, vect(g))[,2])
n_floored <- length(intersect(FC[fcover_true < 0.5]$id_plot, V$id_plot))
#' Summer-mean DeltaTmax at 1 m for each re-run logger of one configuration
#' @param dir directory of per-logger NetCDFs (floored baseline or true-cover re-run)
#' @return named numeric, one DeltaTmax per id in `ids`, NA where the run is unreadable
sc <- function(dir) vapply(ids, function(id) {
  m <- micro_hourly_at(file.path(dir, sprintf("musica_out_HOBO_%s.nc", id)), Z)
  if (is.null(m)) NA_real_ else delta_tmax_mean(m, MREF, ds) }, numeric(1))
D <- data.table(id_plot = ids, fcover_true = FC[match(ids, id_plot)]$fcover_true,
                sim_floored = sc(FL), sim_truecover = sc(NF),
                obs = V[match(ids, V$id_plot)]$obs_dt)
fwrite(D, "out_files/Chapter1/tables/tab_cover_floor_check.csv")
cat(sprintf("loggers below the 0.5 floor: %d (of 53) | re-run at true cover: %d\n", n_floored, nrow(D)))
print(D[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat(sprintf("\nmean  sim floored %+.2f  ->  sim true cover %+.2f   |  observed %+.2f degC\n",
            mean(D$sim_floored), mean(D$sim_truecover), mean(D$obs)))
cat("DONE\n")
