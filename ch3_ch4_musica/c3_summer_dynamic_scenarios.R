# ==============================================================================
# Chapter 3 — SIMPLE summer scenarios with the CORRECT temporal logic:
#   S2 scenarios are DYNAMIC (daily interpolated/smoothed S2 LAI series, tb$*),
#   LiDAR scenarios are STATIC (single ALS date). Validated on 53 HOBO, summer.
# Builds/runs the 5 not yet simulated (DYN_S2_ATBD & STATIC_ALS already have nc).
#   1 DYN_S2_ATBD (exists)  2 DYN_S2_DOPT  3 STATIC_ALS (exists)  4 STATIC_ALS_DOPT
#   5 DYN_ATBD_rfull (tb$atbd*rA)  6 DYN_ATBD_rdopt (tb$atbd*bT)
#   7 FUSION_H (short -> dynamic S2 opt ; tall -> static LiDAR)
# Run from z_Example root:  Rscript c3_summer_dynamic_scenarios.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(data.table); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot
dopt <- CFG_C3$d_opt_m
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))

# ---- fill the 6 open/sparse plots (NA d_opt products), STATIC side ----
df$LAI_ALS_DOPT[is.na(df$LAI_ALS_DOPT)] <- df$LAI_ALS[is.na(df$LAI_ALS_DOPT)]
ropt <- mean(df$LAI_S2_DOPT, na.rm=TRUE) / mean(df$LAI_S2_ATBD, na.rm=TRUE)
df$LAI_S2_DOPT[is.na(df$LAI_S2_DOPT)] <- df$LAI_S2_ATBD[is.na(df$LAI_S2_DOPT)] * ropt
rA <- mean(df$LAI_ALS)/mean(df$LAI_S2_ATBD)
bT <- mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_ATBD)

# ---- fill the dynamic S2-opt series (tb$dopt) for plots where it is missing ----
get_atbd <- function(p) tb$atbd[[p]]
for (p in df$pid) {
  x <- tb$dopt[[p]]
  if (is.null(x) || all(is.na(x$lai))) {
    a <- get_atbd(p); if (!is.null(a)) tb$dopt[[p]] <- data.frame(doy=a$doy, lai=a$lai*ropt)
  }
}
flat_ts <- function(p, val){ a <- get_atbd(p); doy <- if(is.null(a)||!nrow(a)) 1:365 else a$doy; data.frame(doy=doy, lai=rep(val, length(doy))) }

all_sc <- make_all_scenarios_c3(tb, CFG_C3$list_year, dopt, CFG_C3$scenarios_mode)
.fn  <- function(col) function(pr) as.numeric(pr[[col]])
mk   <- function(ts) make_phenology_fn_factory(ts, CFG_C3$list_year)
lad  <- make_lad_real

# DYNAMIC "S2 opt" = the forest-tuned opt LUT magnitude (per plot) carried by the
# ATBD temporal shape: opt_dyn(t) = atbd(t) * (LAI_S2_DOPT[plot] / LAI_S2_ATBD[plot]).
ropt_p <- setNames(df$LAI_S2_DOPT / df$LAI_S2_ATBD, df$pid)
# NULL-safe daily series: if a plot has no S2 series (6 open plots), use a flat
# series at its static value so it still runs (53 everywhere).
safe_ser <- function(p, scale, fallback_val){
  a <- get_atbd(p)
  if (is.null(a) || !all(c("doy","lai") %in% names(a)) || nrow(a)==0)
    return(data.frame(doy=1:365, lai=rep(fallback_val, 365)))
  data.frame(doy=a$doy, lai=a$lai*scale)
}
ts_opt   <- setNames(lapply(seq_len(nrow(df)), function(i) safe_ser(df$pid[i], ropt_p[[df$pid[i]]], df$LAI_S2_DOPT[i])), df$pid)
ts_rfull <- setNames(lapply(seq_len(nrow(df)), function(i) safe_ser(df$pid[i], rA, df$LAI_S2_ATBD[i]*rA)), df$pid)
ts_rdopt <- setNames(lapply(seq_len(nrow(df)), function(i) safe_ser(df$pid[i], bT, df$LAI_S2_ATBD[i]*bT)), df$pid)
# S2 OPT rescaled by global ratio to the LiDAR d_opt magnitude
bod <- mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_DOPT)
ts_opt_rdopt <- setNames(lapply(df$pid, function(p){ a<-ts_opt[[p]]; data.frame(doy=a$doy, lai=a$lai*bod) }), df$pid)
# fully DYNAMIC fusion: short -> dynamic S2 opt (LUT) ; tall -> dynamic LiDAR
# (tb$annual = LAI_ALS * S2 phenology shape). Fallback to flat LAI_ALS if missing.
dyn_als <- function(p, val){ a<-tb$annual[[p]]; if(is.null(a)||all(is.na(a$lai))) flat_ts(p,val) else a }
ts_fusion<- setNames(lapply(seq_len(nrow(df)), function(i){ p<-df$pid[i]
  if (df$Hmax[i] < dopt) ts_opt[[p]] else dyn_als(p, df$LAI_ALS[i]) }), df$pid)

sc <- list(
  DYN_S2_ATBD     = all_sc[["DYN_S2_ATBD"]],
  DYN_S2_OPT      = list(name="DYN_S2_OPT",     lai_fn=.fn("LAI_S2_DOPT"),  hmax_fn=.fn("Hmax"),
                         fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_opt)),
  STATIC_ALS      = all_sc[["STATIC_ALS"]],
  STATIC_ALS_DOPT = all_sc[["STATIC_ALS_DOPT"]],
  DYN_ATBD_rfull  = list(name="DYN_ATBD_rfull", lai_fn=.fn("LAI_ALS"),      hmax_fn=.fn("Hmax"),
                         fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_rfull)),
  DYN_ATBD_rdopt  = list(name="DYN_ATBD_rdopt", lai_fn=.fn("LAI_ALS_DOPT"), hmax_fn=.fn("Hmax"),
                         fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_rdopt)),
  FUSION_H        = list(name="FUSION_H",       lai_fn=.fn("LAI_ALS"),      hmax_fn=.fn("Hmax"),
                         fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_fusion)),
  DYN_OPT_rdopt   = list(name="DYN_OPT_rdopt", lai_fn=.fn("LAI_S2_DOPT"), hmax_fn=.fn("Hmax"),
                         fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_opt_rdopt)))
cat(sprintf("rA=%.3f bT=%.3f ropt=%.3f | plots Hmax<d_opt (fusion uses S2 opt): %d\n",
            rA, bT, ropt, sum(df$Hmax<dopt)))

ds  <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm  <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd  <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
ncp <- file.path(CFG_C3$out_dir, "nc")
val <- validate_scenarios_at_hobos(df, hd, sc, ncp, dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
cat("\nDONE dynamic-scenario MuSICA runs\n")
