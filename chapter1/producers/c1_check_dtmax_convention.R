# ==============================================================================
# Guard: ΔTmax must use ONE convention everywhere in Chapter 1.
#
# WHY THIS EXISTS. Two extractors live in this repository. The canonical one is the
# time-matched convention B of R/dtmax_convention.R (the sub-canopy value is read AT
# the hour of the macroclimatic daily maximum; observations are additionally sampled
# one hour earlier, OBS_CLOCK_OFFSET_H, because the loggers are in UTC and the forcing
# is on solar time = UTC+1). The other is the LEGACY extract_deltatmax_one() in
# R/musica.R, which takes the sub-canopy daily maximum independently and applies a
# -2 h shift. On 2026-07-31 Fig. B1 was found to be the last shipped figure still on
# the legacy extractor, which is exactly the kind of drift this script now catches.
#
# Three checks, all cheap:
#   1. no script run by run_chapter1.R calls the legacy extractor;
#   2. the observation reference really is the simulation reference shifted by -1 h;
#   3. the headline validation still lands on the aligned numbers, not the unaligned
#      ones (an unaligned clock inflates the warm bias from about +0.98 to +1.51 degC).
#   Rscript scripts/c1_check_dtmax_convention.R
# Reads :
#         run_chapter1.R and every stage script it names   (grep for the legacy extractor)
#         in_files/FR-Blo_2021_v2.nc                       (both macro references)
#         out_files/Chapter1/tables/tab_hobo_native20_validation.csv (stage A1)
# Writes: nothing. It is a guard: it prints, and stops non-zero on any failure.
#   Rscript scripts/c1_check_dtmax_convention.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(lubridate)})
source("R/dtmax_convention.R")
fail <- character()

# ---- 1. no shipped stage may use the legacy extractor -------------------------
runner <- readLines(here::here("chapter1/orchestration/run_chapter1.R"), warn = FALSE)
staged <- unique(na.omit(regmatches(runner, regexpr("[A-Za-z0-9_/]+\\.R", runner))))
staged <- staged[file.exists(staged) & basename(staged) != here::here("chapter1/orchestration/run_chapter1.R")]
staged <- staged[basename(staged) != "c1_check_dtmax_convention.R"]   # this file names it on purpose
offend <- Filter(function(f) {
  txt <- readLines(f, warn = FALSE)
  any(grepl("extract_deltatmax_one", txt) & !grepl("^\\s*#", txt))   # calls, not comments
}, staged)
if (length(offend))
  fail <- c(fail, paste0("legacy extractor called by a shipped stage: ",
                         paste(offend, collapse = ", ")))
cat(sprintf("  [%s] %d shipped stages, none calls extract_deltatmax_one()\n",
            if (!length(offend)) "x" else " ", length(staged)))

# ---- 2. the observation reference is the simulation reference, shifted --------
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
ms <- macro_ref("in_files/FR-Blo_2021_v2.nc", ds)
mo <- macro_ref_obs("in_files/FR-Blo_2021_v2.nc", ds)
d  <- merge(ms[, .(date, t_sim = t_max)], mo[, .(date, t_obs = t_max)], by = "date")
off <- unique(as.numeric(difftime(d$t_obs, d$t_sim, units = "hours")))
if (!identical(off, as.numeric(OBS_CLOCK_OFFSET_H)))
  fail <- c(fail, sprintf("obs/sim reference offset is %s h, expected %d",
                          paste(off, collapse = "/"), OBS_CLOCK_OFFSET_H))
cat(sprintf("  [%s] observation reference = simulation reference %+d h, on all %d days\n",
            if (length(off) == 1) "x" else " ", OBS_CLOCK_OFFSET_H, nrow(d)))

# ---- 3. the headline validation is on the aligned clock ----------------------
V <- fread("out_files/Chapter1/tables/tab_hobo_native20_validation.csv")
bias <- mean(V$sim_dt - V$obs_dt, na.rm = TRUE)
r    <- cor(V$sim_dt, V$obs_dt, use = "complete.obs")
rsl  <- cor(V$sim_sl, V$obs_sl, use = "complete.obs")
ok3  <- abs(bias - 0.98) < 0.05 && abs(r - 0.50) < 0.03 && abs(rsl - 0.78) < 0.03
if (!ok3) fail <- c(fail, sprintf("headline drifted: bias %+.3f, r %.3f, slope r %.3f", bias, r, rsl))
cat(sprintf("  [%s] headline n = %d, bias %+.3f degC, r %.3f, slope r %.3f\n",
            if (ok3) "x" else " ", nrow(V), bias, r, rsl))

if (length(fail)) stop("ΔTmax convention check FAILED:\n  ", paste(fail, collapse = "\n  "))
cat("\nOne ΔTmax convention throughout.\nDONE\n")
