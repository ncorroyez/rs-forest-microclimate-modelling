# ==============================================================================
# Chapter 3 — generate the FULL scenario matrix at 53/53 under v3.2.3-iter, cleanly.
# Reuses the PROVEN wrapper machinery (run_v323iter_ch3_full/extra.R): 2020/366
# seed injection (wrap_iter), DOPT/near-bare fallbacks on df, iter binary + PBLH
# forcing, parallel + resumable (skips already-valid nc). Writes into cache B
# (out_files/Chapter3/nc_v323iter), which already holds 14 scenarios at 53/53;
# this fills in the remaining make_all_scenarios_c3 variants.
#
# Modes (env):
#   DRY=1  -> 3 representative scenarios x the 8 previously-failing plots into
#             nc_v323iter_drytest (verify valid nc + timing before the full run)
#   (none) -> all scenarios x 53 plots into nc_v323iter (resumable)
#   Rscript gen_all_scenarios_v323.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC))   # RELATIVE (proven iter forcing)
ABL  <- list("abl_flag" = '"iter"')
LY   <- 2020:2022
NCORES <- 10L
DRY  <- nzchar(Sys.getenv("DRY"))

inject_seed <- function(ph) {
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]
  if (nrow(d) != 1L) stop("missing 2020/365 row to seed 2020/366")
  d$Julian_day <- 366L; rbind(ph, d)
}
wrap_iter <- function(sc) {
  orig <- sc$phenology_fn; lai_fn <- sc$lai_fn
  sc$phenology_fn <- function(pr) {
    ph <- if (!is.null(orig)) orig(pr) else NULL
    if (is.null(ph)) ph <- calc_phenology(list.year = LY, nleafage = 1, budburst_date = 115,
          leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
          LAI_max_per_cohort = as.numeric(lai_fn(pr)))
    inject_seed(ph)
  }
  sc
}

# ---- data + fallbacks (identical to the proven wrapper) ----------------------
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot; dopt <- CFG_C3$d_opt_m
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
df$LAI_ALS_DOPT[is.na(df$LAI_ALS_DOPT)] <- df$LAI_ALS[is.na(df$LAI_ALS_DOPT)]     # near-bare -> LiDAR full
ropt <- mean(df$LAI_S2_DOPT, na.rm=TRUE)/mean(df$LAI_S2_ATBD, na.rm=TRUE)
df$LAI_S2_DOPT[is.na(df$LAI_S2_DOPT)] <- df$LAI_S2_ATBD[is.na(df$LAI_S2_DOPT)]*ropt
if ("LAI_S2_RESCALED" %in% names(df))                                             # near-bare -> LiDAR full
  df$LAI_S2_RESCALED[is.na(df$LAI_S2_RESCALED)] <- df$LAI_ALS[is.na(df$LAI_S2_RESCALED)]

# ---- full scenario matrix, seed-wrapped --------------------------------------
all_sc <- make_all_scenarios_c3(tb, LY, dopt, "full")
SC <- lapply(all_sc, wrap_iter)

NCROOT <- file.path(CFG_C3$out_dir, "nc_v323iter")
prob8  <- c("41_17","41_18","41_19","41_27","41_30","41_39","41_47","41_49")
if (DRY) {
  NCROOT <- file.path(CFG_C3$out_dir, "nc_v323iter_drytest")
  SC <- SC[intersect(c("DYN_ALS","DYN_S2_DOPT","STATIC_S2_DOPT"), names(SC))]
  idx <- which(df$id_plot %in% prob8)
  cat(sprintf("** DRY: %d scenarios x %d failing plots -> %s **\n", length(SC), length(idx), NCROOT))
} else {
  idx <- seq_len(nrow(df))
  cat(sprintf("FULL: %d scenarios x %d plots -> %s (resumable)\n", length(SC), length(idx), NCROOT))
}
for (s in names(SC)) dir.create(file.path(NCROOT, s), recursive=TRUE, showWarnings=FALSE)
tasks <- rbindlist(lapply(names(SC), function(s) data.table(scn=s, i=idx)))

run_one <- function(k) {
  s <- tasks$scn[k]; i <- tasks$i[k]; prow <- as.data.frame(df[i,])
  nc <- file.path(NCROOT, s, sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(nc) && file.size(nc) > 1000) return(sprintf("skip %s/%s", s, prow$id_plot))
  if (file.exists(nc)) file.remove(nc)      # drop a previous failed stub
  ok <- tryCatch({ run_musica_one(prow, SC[[s]], nc, FORC, BIN, extra_setup = ABL)
                   file.exists(nc) && file.size(nc) > 1000 }, error = function(e) FALSE)
  sprintf("%s %s/%s", if (ok) "OK" else "FAIL", s, prow$id_plot)
}
cat(sprintf("launching %d runs on %d cores\n", nrow(tasks), NCORES))
t0 <- Sys.time()
res <- unlist(mclapply(seq_len(nrow(tasks)), run_one, mc.cores = NCORES, mc.preschedule = FALSE))
dt <- as.numeric(difftime(Sys.time(), t0, units="mins"))
tab <- table(sub(" .*", "", res))
cat(sprintf("\nDONE in %.1f min | %s | ~%.1f s/run\n", dt, paste(names(tab), tab, sep="=", collapse=" "),
            60*dt/max(1,sum(!grepl("^skip", res)))))
fails <- grep("^FAIL", res, value=TRUE); if (length(fails)) { cat("FAILURES:\n"); cat(fails, sep="\n"); cat("\n") }
if (!DRY) saveRDS(res, file.path(NCROOT, "_run_status_full.rds"))
