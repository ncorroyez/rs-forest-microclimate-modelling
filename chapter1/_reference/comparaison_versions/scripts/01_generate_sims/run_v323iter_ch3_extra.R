# ==============================================================================
# Extra scenarios for the v3.2.3-iter robustness control: the narrative
# scenarios not in the 8-panel monthly figure but with an existing v3.2.0 nc:
#   DYN_RF, DYN_S2_ANNUAL, CONST_ALS, NAIVE_S2_FORMSH (from make_all_scenarios_c3),
#   FUSION_H (layered §3.5: short->dyn S2 opt, tall->dyn LiDAR; built like
#   c3_summer_dynamic_scenarios.R). Same iter setup + 2020/366 seed as the main run.
# Each scenario is reconstructed on the SAME series basis as its v3.2.0 nc (the
# masked tb=load_lai_prep basis), so v3.2.0 vs v3.2.3 differ only by the binary.
#   Rscript run_v323iter_ch3_extra.R   (writes out_files/Chapter3/nc_v323iter/<sc>/)
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC))
ABL  <- list("abl_flag" = '"iter"')
LY   <- 2020:2022
NCROOT <- file.path(CFG_C3$out_dir, "nc_v323iter"); dir.create(NCROOT, recursive=TRUE, showWarnings=FALSE)
NCORES <- 10L

inject_seed <- function(ph) {
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]
  if (nrow(d) != 1L) stop("missing 2020/365 row to seed 2020/366")
  d$Julian_day <- 366L; rbind(ph, d)
}
# generic: wrap any scenario's phenology_fn so the output carries the 2020/366
# seed; NULL (static fallback) -> parametric one-sided curve (run_musica_one doubles)
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

# ---- production series basis (masked tb), as the v3.2.0 nc were built ---------
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot; dopt <- CFG_C3$d_opt_m
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
df$LAI_ALS_DOPT[is.na(df$LAI_ALS_DOPT)] <- df$LAI_ALS[is.na(df$LAI_ALS_DOPT)]
ropt <- mean(df$LAI_S2_DOPT, na.rm=TRUE)/mean(df$LAI_S2_ATBD, na.rm=TRUE)
df$LAI_S2_DOPT[is.na(df$LAI_S2_DOPT)] <- df$LAI_S2_ATBD[is.na(df$LAI_S2_DOPT)]*ropt

# build the full orchestrator set with list_year=2020:2022 (so DYN phenology spans 2020)
all_sc <- make_all_scenarios_c3(tb, LY, dopt, "full")
WANT <- c("DYN_RF","DYN_S2_ANNUAL","CONST_ALS","NAIVE_S2_FORMSH")
SC <- all_sc[WANT]

# ---- FUSION_H (layered §3.5), built on the same masked basis ------------------
.fn <- function(col) function(pr) as.numeric(pr[[col]]); lad <- make_lad_real
mk  <- function(ts) make_phenology_fn_factory(ts, LY)
get_atbd <- function(p) tb$atbd[[p]]
flat_ts  <- function(p, val){ a<-get_atbd(p); doy<-if(is.null(a)||!nrow(a))1:365 else a$doy; data.frame(doy=doy, lai=rep(val,length(doy))) }
ropt_p   <- setNames(df$LAI_S2_DOPT/df$LAI_S2_ATBD, df$pid)
safe_ser <- function(p, scale, fb){ a<-get_atbd(p)
  if (is.null(a) || !all(c("doy","lai") %in% names(a)) || nrow(a)==0) return(data.frame(doy=1:365, lai=rep(fb,365)))
  data.frame(doy=a$doy, lai=a$lai*scale) }
ts_opt  <- setNames(lapply(seq_len(nrow(df)), function(i) safe_ser(df$pid[i], ropt_p[[df$pid[i]]], df$LAI_S2_DOPT[i])), df$pid)
dyn_als <- function(p, val){ a<-tb$annual[[p]]; if(is.null(a)||all(is.na(a$lai))) flat_ts(p,val) else a }
ts_fus  <- setNames(lapply(seq_len(nrow(df)), function(i){ p<-df$pid[i]
  if (df$Hmax[i] < dopt) ts_opt[[p]] else dyn_als(p, df$LAI_ALS[i]) }), df$pid)
SC$FUSION_H <- list(name="FUSION_H", lai_fn=.fn("LAI_ALS"), hmax_fn=.fn("Hmax"),
                    fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_fus))

SC <- lapply(SC, wrap_iter)
cat(sprintf("extra scenarios: %s\n", paste(names(SC), collapse=", ")))
for (s in names(SC)) dir.create(file.path(NCROOT, s), recursive=TRUE, showWarnings=FALSE)

SMOKE <- nzchar(Sys.getenv("SMOKE"))
idx <- if (SMOKE) 1L else seq_len(nrow(df))   # SMOKE: 1 plot per scenario
tasks <- rbindlist(lapply(names(SC), function(s) data.table(scn=s, i=idx)))
if (SMOKE) { NCROOT <- file.path(CFG_C3$out_dir, "nc_v323iter_smoke"); for (s in names(SC)) dir.create(file.path(NCROOT,s),recursive=TRUE,showWarnings=FALSE); cat("** SMOKE MODE **\n") }
cat(sprintf("launching %d runs (%d scenarios x %d plots) on %d cores\n", nrow(tasks), length(SC), nrow(df), NCORES))
run_one <- function(k) {
  s <- tasks$scn[k]; i <- tasks$i[k]; prow <- as.data.frame(df[i,])
  nc <- file.path(NCROOT, s, sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(nc) && file.size(nc) > 1000) return(sprintf("skip %s/%s", s, prow$id_plot))
  if (file.exists(nc)) file.remove(nc)
  ok <- tryCatch({ run_musica_one(prow, SC[[s]], nc, FORC, BIN, extra_setup = ABL)
                   file.exists(nc) && file.size(nc) > 1000 }, error = function(e) FALSE)
  sprintf("%s %s/%s", if (ok) "OK" else "FAIL", s, prow$id_plot)
}
t0 <- Sys.time()
res <- unlist(mclapply(seq_len(nrow(tasks)), run_one, mc.cores = NCORES, mc.preschedule = FALSE))
dt <- as.numeric(difftime(Sys.time(), t0, units="mins"))
tab <- table(sub(" .*", "", res))
cat(sprintf("\nDONE in %.1f min | %s\n", dt, paste(names(tab), tab, sep="=", collapse=" ")))
fails <- grep("^FAIL", res, value=TRUE); if (length(fails)) { cat("FAILURES:\n"); cat(fails, sep="\n"); cat("\n") }
