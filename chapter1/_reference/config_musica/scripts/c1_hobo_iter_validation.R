# ==============================================================================
# HOBO validation in YOYO (iter) mode: run the 53 HOBO at REF (real z05 canopy)
# with v3.2.3 + ABL_flag='iter' + real MERRA-2 PBLH, then compare per-plot r/bias
# vs HOBO obs against the existing legacy v3.2.0 (r~0.93) and v3.2.3-none (r~0.63)
# REF runs. Decides whether the yoyo rescues the field validation (it shouldn't:
# iter≈none on cLHS). REF only (53 sims), reuses z05 inputs.
#   Rscript c1_hobo_iter_validation.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(sf)
  library(terra); library(parallel); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("pipeline/00_config.R")
MB   <- normalizePath("in_files/model-3.2.3/musica", mustWork=TRUE)
FORC <- "in_files/musica_in_Blois_pblh.nc"

# ---- z05 HOBO REF inputs (mirror pipeline/03z_musica_z05.R lines 26-44) -------
z05 <- fread("in_files/lad_z05/Blois_lad_z05_r25.csv")
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
hob_r <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
hp <- sf::st_read(CFG$hobo_geojson, quiet=TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
hob_r$id_plot <- hp$id_plot
fcov <- hob_r[, .(id_plot, fCover)]
df <- merge(z05, fcov, by="id_plot"); df[fCover < 0.5, fCover := 0.5]; df <- as.data.frame(df)
cat(sprintf("HOBO plots: %d\n", nrow(df)))

mk_phen <- function(lai1) function(p) {
  ph <- as.data.frame(calc_phenology(list.year=2020:2022, nleafage=1, budburst_date=115,
        leaf_age_max_in=0.56, relative_age_firstmax=0.10, relative_age_lastmax=0.75, LAI_max_per_cohort=lai1))
  d <- ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]; d$Julian_day <- 366; rbind(ph, d)
}
OUT <- "out_files/musica_hobo_z05_iter/1111"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# ---- run REF iter for all 53 (resumable) -------------------------------------
jobs <- lapply(seq_len(nrow(df)), function(i) df[i,,drop=FALSE])
invisible(mclapply(jobs, function(pr) {
  nc <- file.path(OUT, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
  if (file.exists(nc) && file.size(nc) > 1e6) return(NULL)
  sc <- list(lai_fn=function(p) pr$LAI, hmax_fn=function(p) pr$Hmax, fcover_fn=function(p) pr$fCover,
             lad_fn=make_lad_real, phenology_fn=mk_phen(pr$LAI))
  tryCatch(run_musica_one(pr, sc, nc, FORC, MB, extra_setup=list("abl_flag"='"iter"')),
           error=function(e) cat(sprintf("ERR %s: %s\n", pr$id_plot, e$message)))
}, mc.cores=4, mc.preschedule=FALSE))
cat(sprintf("iter REF nc: %d/53\n", length(list.files(OUT,"\\.nc$"))))

# ---- consistent ΔTmax extraction (1 m, -2 h, JJAS, vs macro) -----------------
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- as.data.table(extract_macro_daily(CFG$forcing_file, ds)); Z<-1; SH<-2L
dT <- function(p){ if(!file.exists(p)||file.size(p)<1e5) return(NA_real_)
  nc<-try(nc_open(p),silent=TRUE); if(inherits(nc,"try-error")) return(NA_real_); on.exit(nc_close(nc))
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%names(nc$var))) return(NA_real_)
  tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC")
  th<-ncvar_get(nc,"time"); Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height")
  vh<-stats::median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE); zl<-rh*vh
  if(Z<=zl[1]){il<-1L;ih<-1L;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1L;w<-(Z-zl[il])/(zl[ih]-zl[il])}
  tv<-t0+dhours(th)-lubridate::hours(SH); Tc<-((1-w)*Tk[il,]+w*Tk[ih,])-273.15
  dd<-data.table(date=as.Date(floor_date(tv,"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date]
  merge(dd,dm,by="date")[,mean(Tmax-Tmax_macro,na.rm=TRUE)] }

VDIR <- c("v3.2.0"="out_files/musica_hobo_z05/1111",
          "v3.2.3 none"="out_files/musica_hobo_z05_v323/1111",
          "v3.2.3 iter"="out_files/musica_hobo_z05_iter/1111")
ids <- df$id_plot
sim <- sapply(VDIR, function(d) sapply(ids, function(id) dT(file.path(d, sprintf("musica_out_HOBO_%s.nc", id)))))
# observed per-plot ΔTmax
V <- readRDS("outputs/figures_pipeline_z05/data/ref_validation.rds")
obs <- as.data.table(V$obs_daily)[, .(obs=mean(Delta_obs,na.rm=TRUE)), by=id_plot]
obs <- obs[match(ids, id_plot)]$obs

cat("\n=== HOBO validation per version (53 plots, ΔTmax) ===\n")
for (v in colnames(sim)) {
  s <- sim[,v]; ok <- is.finite(s)&is.finite(obs)
  cat(sprintf("  %-12s r=%.3f  bias=%+.2f  RMSE=%.2f  n=%d\n", v,
      cor(s[ok],obs[ok]), mean(s[ok]-obs[ok]), sqrt(mean((s[ok]-obs[ok])^2)), sum(ok)))
}
cat("\n(legacy v3.2.0 r~0.93 is the validated reference; does iter beat v3.2.3 none r~0.63?)\n")
fwrite(data.table(id_plot=ids, obs=obs, sim), "out_files/Chapter3/tables/tab_hobo_iter_validation.csv")
cat("DONE -> tab_hobo_iter_validation.csv\n")
