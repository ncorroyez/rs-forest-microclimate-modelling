# ==============================================================================
# Forward-inclusion coalition sims under v3.2.3 YOYO/ITER + STATION forcing (+h_sbl).
# 53 HOBO × 16 coalitions (2^4 lattice of LAI/Hmax/fCover/LAD), isolated dir.
# ITER mode (abl_flag='iter', h_sbl) + phenology seed 2020/366 (required by iter).
# Then extract ΔTmax(bit,logger) vs station macro → coal_metrics_v323station.rds.
#   NMAX=1 Rscript c1_forward_v323_sims.R   (gate) ; then full.
# ==============================================================================
suppressMessages({ library(here); library(data.table); library(terra); library(parallel)
                   library(ncdf4); library(lubridate); library(rmusica); library(musica.tools) })
source(here::here("pipeline/00_config.R"))
MB   <- normalizePath("in_files/model-3.2.3/musica", mustWork=TRUE)
FORC <- "in_files/FR-Blo_2021_v2.nc"                                       # FR-Blo (has own h_sbl → iter)
OUTD <- "out_files/musica_hobo_z05_v323station_iter"                   # isolated coalition dir (iter)
ABL  <- list("abl_flag" = '"iter"')
stopifnot(file.exists(FORC))
# phenology seed with the 2020/366 day-before row (iter needs it); equivalent to
# the default fallback (2*LAI) + the extra seed row. LAI = the coalition's one-sided LAI.
mk_phen <- function(lai1) { force(lai1)          # force NOW: else lazy-eval captures the last loop pr
  function(p) {
  ph <- as.data.frame(calc_phenology(list.year=2020:2022, nleafage=1, budburst_date=115,
        leaf_age_max_in=0.56, relative_age_firstmax=0.10, relative_age_lastmax=0.75, LAI_max_per_cohort=lai1))
  d <- ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]; d$Julian_day <- 366; rbind(ph, d) } }

bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
df_floor <- as.data.table(readRDS(PIPE$CLUSTER_SAMPLE))
fac_scs  <- build_factorial_scenarios_archetypes(df_floor, fcov_b = PIPE$FCOV_BASELINE)
rasters  <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo  <- as.data.frame(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove, hobo_buffer_mode="buffer25"))  # 10 m PAD
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5
NMAX <- as.integer(Sys.getenv("NMAX", nrow(df_hobo)))
df_hobo <- df_hobo[seq_len(min(NMAX, nrow(df_hobo))), , drop=FALSE]

jobs <- list()
for (i in seq_len(nrow(df_hobo))) { pr <- df_hobo[i,,drop=FALSE]
  for (bit in PIPE$BITS) {
    od <- file.path(OUTD, bit); dir.create(od, recursive=TRUE, showWarnings=FALSE)
    onc <- file.path(od, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(onc) && file.size(onc) > 1e6) next
    sc <- fac_scs[[bit2scn[bit]]]
    sc$phenology_fn <- mk_phen(sc$lai_fn(pr))          # iter needs the 2020/366 seed
    jobs[[length(jobs)+1L]] <- list(pr=pr, sc=sc, out_nc=onc, lbl=sprintf("%s/%s",bit,pr$id_plot))
  } }
cat(sprintf("coalition sims to run: %d (of %d) [ITER]\n", length(jobs), length(PIPE$BITS)*nrow(df_hobo)))
if (length(jobs) > 0) {
  t0 <- Sys.time()
  mclapply(jobs, function(j) tryCatch(run_musica_one(j$pr, j$sc, j$out_nc, FORC, MB, extra_setup=ABL),
           error=function(e) cat(sprintf("[%s] ERR: %s\n", j$lbl, e$message))), mc.cores=4, mc.preschedule=FALSE)
  cat(sprintf("ran %d sims in %.1f min\n", length(jobs), as.numeric(difftime(Sys.time(),t0,units="mins"))))
}
cnt <- sapply(PIPE$BITS, function(b) length(list.files(file.path(OUTD,b),"\\.nc$")))
print(data.frame(bit=PIPE$BITS, n=cnt))

# ---- extract ΔTmax(bit, logger) vs station macro -----------------------------
if (all(cnt >= nrow(df_hobo))) {
  ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day"); Z<-1; SH<-2L
  dm <- as.data.table(extract_macro_daily(FORC, ds))
  dTmax <- function(p){ if(!file.exists(p)||file.size(p)<1e5) return(NA_real_)
    nc<-try(nc_open(p),silent=TRUE); if(inherits(nc,"try-error")) return(NA_real_); on.exit(nc_close(nc))
    if(!all(c("Tair_z","relative_height","veget_height_top")%in%names(nc$var))) return(NA_real_)
    tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC")
    th<-ncvar_get(nc,"time"); Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height")
    vh<-stats::median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE); zl<-rh*vh
    if(Z<=zl[1]){il<-1L;ih<-1L;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1L;w<-(Z-zl[il])/(zl[ih]-zl[il])}
    tv<-t0+dhours(th)-lubridate::hours(SH); Tc<-((1-w)*Tk[il,]+w*Tk[ih,])-273.15
    dd<-data.table(date=as.Date(floor_date(tv,"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date]
    merge(dd,dm,by="date")[,mean(Tmax-Tmax_macro,na.rm=TRUE)] }
  ids <- df_hobo$id_plot
  rows <- list()
  for (bit in PIPE$BITS) for (id in ids) {
    v <- dTmax(file.path(OUTD, bit, sprintf("musica_out_HOBO_%s.nc", id)))
    rows[[length(rows)+1L]] <- data.table(id_plot=id, bit=bit, Delta_sim=v) }
  coal <- rbindlist(rows)
  saveRDS(coal, "out_files/Chapter1/coal_metrics_v323station_iter.rds")
  cat(sprintf("\nWROTE coal_metrics_v323station_iter.rds (%d rows, %d loggers x %d bits)\n",
              nrow(coal), uniqueN(coal$id_plot), uniqueN(coal$bit)))
}
