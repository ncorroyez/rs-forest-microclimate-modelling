# ==============================================================================
# Chapter 3 — the extinction-coefficient counterfactual, k = 0.5 -> k = 0.65.
#
# WHY. Chapter 1 states twice that a uniform 23 % reduction in retrieved leaf
# area "changes the simulated coupling", and refers the reader to this chapter
# for the evidence. No file carried that result. This generator produces it.
#
# WHAT. Beer-Lambert gives LAI = -ln(P_gap)/k, so LAI scales as 1/k. The whole
# Chapter 3 tree runs on lai_z1_res_10_m.tif, produced by myPAI(Z, zmin = 1)
# with the k = 0.5 default of 02_CODES/libraries/functions_lidar2.R (verified:
# the in_files copy is byte-identical to 03_RESULTS/Blois/Metrics/Raw). Moving
# to k = 0.65 therefore multiplies every leaf area by 0.5/0.65 = 0.769, a
# uniform 23.1 % reduction. Shape, Hmax and fCover are untouched: make_lad_real
# rescales the PAD profile to whatever LAI it is given, so overriding lai_fn is
# exactly the uniform rescale Chapter 1 describes.
#
# Only STATIC_ALS is run. The Sentinel-2 scenarios do not carry a LiDAR k.
# Same harness as gen_opt_CHS41_nowind.R: CHS41-Rmerge forcing as-is, no wind
# correction, ABL iter with the 2020/366 phenology seed. Writes to its own tree;
# nothing is deleted or overwritten.
#   DRY=1 Rscript gen_k065_CHS41.R   -> 3 plots (open/median/dense) -> nc_k065_DRY
#         Rscript gen_k065_CHS41.R   -> 53 plots -> nc_k065/STATIC_ALS_K065
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config_CHS41.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC <- CFG_C3$forcing_file; stopifnot(file.exists(FORC))
ABL  <- list("abl_flag" = '"iter"'); LY <- 2020:2022; NCORES <- 8L
DRY  <- nzchar(Sys.getenv("DRY"))

K_BASE <- 0.5; K_NEW <- 0.65; SCALE <- K_BASE / K_NEW
SCN    <- "STATIC_ALS_K065"
NCROOT <- if (DRY) file.path(CFG_C3$out_dir, "nc_k065_DRY") else file.path(CFG_C3$out_dir, "nc_k065")
dir.create(file.path(NCROOT, SCN), recursive = TRUE, showWarnings = FALSE)

inject_seed <- function(ph){ d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]
  if (nrow(d) != 1L) stop("missing 2020/365"); d$Julian_day <- 366L; rbind(ph, d) }
wrap_iter <- function(sc){ orig <- sc$phenology_fn; lai_fn <- sc$lai_fn
  sc$phenology_fn <- function(pr){ ph <- if (!is.null(orig)) orig(pr) else NULL
    if (is.null(ph)) ph <- calc_phenology(list.year = LY, nleafage = 1, budburst_date = 115,
        leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
        LAI_max_per_cohort = as.numeric(lai_fn(pr)))
    inject_seed(ph) }; sc }

df <- as.data.table(readRDS(file.path(CFG_C3$out_dir, "lai_prep", "df_plots_real53.rds")))
stopifnot(nrow(df) == 53L, "LAI_ALS" %in% names(df))

.fn <- function(col) function(pr) as.numeric(pr[[col]])
SC <- wrap_iter(list(name = SCN,
                     lai_fn    = function(pr) as.numeric(pr[["LAI_ALS"]]) * SCALE,
                     hmax_fn   = .fn("Hmax"), fcover_fn = .fn("fCover"),
                     lad_fn    = make_lad_real, phenology_fn = NULL))

cat(sprintf("k %.2f -> %.2f : leaf area x %.5f (a %.1f %% reduction)\n",
            K_BASE, K_NEW, SCALE, 100 * (1 - SCALE)))
cat(sprintf("LAI_ALS  median %.3f -> %.3f   range [%.3f, %.3f] -> [%.3f, %.3f]\n",
            median(df$LAI_ALS), median(df$LAI_ALS) * SCALE,
            min(df$LAI_ALS), max(df$LAI_ALS), min(df$LAI_ALS) * SCALE, max(df$LAI_ALS) * SCALE))

if (DRY) { o <- order(df$LAI_ALS); idx <- unique(c(o[1], o[ceiling(nrow(df)/2)], o[nrow(df)]))
} else idx <- seq_len(nrow(df))

run_one <- function(i){ prow <- as.data.frame(df[i, ])
  nc <- file.path(NCROOT, SCN, sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(nc) && file.size(nc) > 1e6) return(sprintf("skip %s", prow$id_plot))
  err <- NULL
  ok <- tryCatch({ run_musica_one(prow, SC, nc, FORC, BIN, extra_setup = ABL)
                   file.exists(nc) && file.size(nc) > 1e6 },
                 error = function(e){ err <<- conditionMessage(e); FALSE })
  sprintf("%s %s (%.2f MB)%s", if (ok) "OK" else "FAIL", prow$id_plot,
          if (file.exists(nc)) file.size(nc)/1e6 else 0,
          if (is.null(err)) "" else paste0("  ERR: ", err)) }

cat(sprintf("\n%s%s: %d plots on %d cores\n", SCN, if (DRY) " [DRY]" else "", length(idx), min(NCORES, length(idx))))
t0 <- Sys.time(); res <- unlist(mclapply(idx, run_one, mc.cores = min(NCORES, length(idx)), mc.preschedule = FALSE))
cat(sprintf("DONE in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(res, sep = "\n"); cat("\n")
