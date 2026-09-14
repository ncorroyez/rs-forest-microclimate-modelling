# ==============================================================================
# Chapter 3 — premise MuSICA numbers under v3.2.3-iter (cache B). STATIC_ALS
# (LiDAR full) summer per-plot buffering slope lm(Tsim ~ Tmacro, 1 m, no shift):
#   (1) r(sim_slope, total LAI)         [article v3.2.0 = -0.94]
#   (2) r(sim_slope, obs_slope)         [article v3.2.0 =  0.96]
#   (3) r(obs_slope, total LAI)         [binary-independent, ~ -0.92]
# Run from z_Example root:  Rscript c3_premise_musica_v323.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(dplyr); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
ncp <- "out_files/Chapter1/nc_v323iter"; scn <- "STATIC_ALS"
ds_all <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")

cov <- as.data.table(readRDS("out_files/Chapter3/lai_prep/df_plots_lai.rds"))[
  , .(id_plot=as.character(id_plot), LAI_ALS)]
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, hr:=floor_date(time,"hour")]
slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_

hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), time, Tobs=t_hobo)]
hb[, hr:=floor_date(time,"hour")]; hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
OBS <- hb[, .(so=slp(Tobs,Tmacro)), by=id_plot]

read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  as.data.table(res)[as.Date(time) %in% ds_all, .(time, Tsim=Tair_sim)] }
fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
SIM <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL)
  sm <- merge(s[,.(time,Tsim,hr=floor_date(time,"hour"))], macro_h[,.(hr,Tmacro)], by="hr")
  data.table(id_plot=id, ss=slp(sm$Tsim,sm$Tmacro)) }), fill=TRUE)

b <- merge(merge(SIM, OBS, by="id_plot"), cov, by="id_plot")
b <- b[is.finite(ss) & is.finite(so) & is.finite(LAI_ALS)]
cat(sprintf("\nPREMISE v3.2.3-iter (STATIC_ALS, summer, n=%d)\n", nrow(b)))
cat(sprintf("  (1) r(sim_slope, LAI_ALS)  = %+.2f   [v3.2.0 = -0.94]\n", cor(b$ss, b$LAI_ALS)))
cat(sprintf("  (2) r(sim_slope, obs_slope)= %+.2f   [v3.2.0 =  0.96]\n", cor(b$ss, b$so)))
cat(sprintf("  (3) r(obs_slope, LAI_ALS)  = %+.2f   [binary-independent]\n", cor(b$so, b$LAI_ALS)))
fwrite(b, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table_premise_slope_v323.csv")
cat("DONE -> Table_premise_slope_v323.csv\n")
