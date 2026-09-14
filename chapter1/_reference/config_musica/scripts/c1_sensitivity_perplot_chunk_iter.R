# ==============================================================================
# iter-mode (yoyo) counterpart of c1_sensitivity_perplot_chunk_v323.R: same 400
# cLHS plots, same per-unit perturbations (ΔLAI/ΔHmax/ΔfCover + LAD real-vs-uniform),
# but MuSICA v3.2.3 with ABL_flag='iter' and REAL MERRA-2 PBLH as h_sbl
# (forcing musica_in_Blois_pblh.nc). iter needs the phenology day-before-sim
# (2020/366), injected via a per-run phenology_fn (one-sided LAI; run_musica_one
# doubles Leaf_area). 8 sims/plot, resumable.
#   Rscript c1_sensitivity_perplot_chunk_iter.R <chunk> <K>
# Out: out_files/Chapter3/tables/sensitivity_perplot_iter/part_<chunk>.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
  library(rmusica); library(musica.tools)
})
args  <- commandArgs(trailingOnly = TRUE)
chunk <- as.integer(args[1]); K <- as.integer(args[2])
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
CFG_C3$musica_cmd <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)  # v3.2.3 official
FORC <- "in_files/musica_in_Blois_pblh.nc"                                          # real MERRA-2 PBLH h_sbl
ABL  <- list("abl_flag" = '"iter"')

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N, min(.N, 100))], by = Cluster]; sub[, pid := sprintf("S%04d", .I)]
my_rows <- which(((seq_len(nrow(sub)) - 1) %% K) == (chunk - 1))
cat(sprintf("iter chunk %d/%d : %d plots\n", chunk, K, length(my_rows)))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
dm <- as.data.table(extract_macro_daily(FORC, ds))   # macro Tair unchanged by added h_sbl
Z_FIX <- 1.0; SHIFT <- 2L
metrics_one <- function(path) {
  if (!file.exists(path) || file.size(path) < 1000) return(NA_real_)
  nc <- try(nc_open(path), silent = TRUE); if (inherits(nc, "try-error")) return(NA_real_)
  on.exit(nc_close(nc))
  if (!all(c("Tair_z","relative_height","veget_height_top") %in% names(nc$var))) return(NA_real_)
  tu <- ncatt_get(nc, "time", "units")$value; t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time"); Tk <- ncvar_get(nc, "Tair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh; Z <- Z_FIX
  if (Z <= zl[1]) { ilo <- 1L; ihi <- 1L; w <- 0 }
  else if (Z >= zl[length(zl)]) { ilo <- length(zl); ihi <- ilo; w <- 0 }
  else { ilo <- max(which(zl <= Z)); ihi <- ilo + 1L; w <- (Z - zl[ilo]) / (zl[ihi] - zl[ilo]) }
  tvec <- t0 + dhours(th) - lubridate::hours(SHIFT); Tc <- ((1 - w) * Tk[ilo, ] + w * Tk[ihi, ]) - 273.15
  dd <- data.table(date = as.Date(floor_date(tvec, "hour")), Tc = Tc)[date %in% ds, .(Tmax = max(Tc)), by = date]
  m1 <- merge(dd, dm, by = "date"); mean(m1$Tmax - m1$Tmax_macro, na.rm = TRUE)
}

# phenology with the day-before-sim (2020/366) row, ONE-SIDED LAI (run_musica_one doubles Leaf_area)
mk_phen <- function(lai_onesided) function(p) {
  ph <- as.data.frame(calc_phenology(list.year = 2020:2022, nleafage = 1, budburst_date = 115,
        leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
        LAI_max_per_cohort = lai_onesided))
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]; d$Julian_day <- 366
  rbind(ph, d)
}

NCDIR <- "out_files/Chapter3/nc_sensitivity_perplot_iter"; dir.create(NCDIR, recursive = TRUE, showWarnings = FALSE)
DLAI <- 2; DHMAX <- 5; DFC <- 0.1; LAI_FLOOR <- 0.1; HMAX_FLOOR <- 3
run_get <- function(prow, tag, lai, hmax, fcov, ladf = make_lad_real) {
  nc <- file.path(NCDIR, sprintf("%s_%s.nc", prow$pid, tag))
  sc <- list(lai_fn = function(p) lai, hmax_fn = function(p) hmax, fcover_fn = function(p) fcov,
             lad_fn = ladf, phenology_fn = mk_phen(lai))
  run_musica_one(prow, sc, nc, FORC, CFG_C3$musica_cmd, extra_setup = ABL)
  if (!file.exists(nc) || file.size(nc) < 1000) return(NA_real_)
  metrics_one(nc)
}

pdir <- "out_files/Chapter3/tables/sensitivity_perplot_iter"; dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
partfile <- file.path(pdir, sprintf("part_%02d.csv", chunk)); if (file.exists(partfile)) file.remove(partfile)
nok <- 0L
for (j in seq_along(my_rows)) {
  prow <- as.data.frame(sub[my_rows[j]])
  r <- tryCatch({
    base <- run_get(prow, "base", prow$LAI, prow$Hmax, prow$fCover)
    unif <- run_get(prow, "unifLAD", prow$LAI, prow$Hmax, prow$fCover, ladf = make_lad_uniform)
    laip_v <- prow$LAI + DLAI;                 laim_v <- max(prow$LAI - DLAI, LAI_FLOOR)
    hmp_v  <- prow$Hmax + DHMAX;               hmm_v  <- max(prow$Hmax - DHMAX, HMAX_FLOOR)
    fcp_v  <- min(prow$fCover + DFC, 1);       fcm_v  <- max(prow$fCover - DFC, 0.5)
    stp_lai_add <- DLAI / 2;                   stp_lai_rem <- (prow$LAI - laim_v) / 2
    stp_hm_add  <- DHMAX;                      stp_hm_rem  <- (prow$Hmax - hmm_v)
    stp_fc_add  <- (fcp_v - prow$fCover) / DFC; stp_fc_rem <- (prow$fCover - fcm_v) / DFC
    laip <- run_get(prow, "LAIp", laip_v, prow$Hmax, prow$fCover)
    laim <- run_get(prow, "LAIm", laim_v, prow$Hmax, prow$fCover)
    hmp  <- run_get(prow, "Hmaxp", prow$LAI, hmp_v, prow$fCover)
    hmm  <- run_get(prow, "Hmaxm", prow$LAI, hmm_v, prow$fCover)
    fcpr <- run_get(prow, "fCovp", prow$LAI, prow$Hmax, fcp_v)
    fcmr <- run_get(prow, "fCovm", prow$LAI, prow$Hmax, fcm_v)
    nrm <- function(x, s) if (is.finite(x) && is.finite(s) && s > 1e-9) x / s else NA_real_
    data.table(pid = prow$pid, Cluster = prow$Cluster, LAI = prow$LAI, Hmax = prow$Hmax, fCover = prow$fCover,
               base = base,
               LAI_add  = nrm(laip - base, stp_lai_add), LAI_rem  = nrm(laim - base, stp_lai_rem),
               Hmax_add = nrm(hmp  - base, stp_hm_add/5), Hmax_rem = nrm(hmm  - base, stp_hm_rem/5),
               fCov_add = nrm(fcpr - base, stp_fc_add),   fCov_rem = nrm(fcmr - base, stp_fc_rem),
               dT_LAD   = base - unif)
  }, error = function(e) NULL)
  if (!is.null(r)) { fwrite(r, partfile, append = file.exists(partfile)); nok <- nok + 1L }
  if (j %% 5 == 0) cat(sprintf("  iter chunk %d: %d/%d\n", chunk, j, length(my_rows)))
}
cat(sprintf("iter chunk %d DONE (%d/%d rows)\n", chunk, nok, length(my_rows)))
