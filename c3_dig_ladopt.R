# ==============================================================================
# Chapter 3 — dig into DYN_ALS_LADOPT (top-d_opt LAD concentration doubles the
# summer ΔTmax ranking). Reads nc_genuine53. (1) LADOPT effect per magnitude
# (full vs top-d_opt LAD, same total LAI). (2) SD-recovery (does concentration
# amplify the between-plot ΔTmax spread toward obs?). (3) jackknife R² of
# DYN_ALS_LADOPT (influential plots?). (4) ranking with the DYN_S2 family.
#   Rscript c3_dig_ladopt.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
d3 <- file.path(CFG_C3$out_dir, "nc_genuine53"); d1 <- "out_files/Chapter1/nc_v323iter"
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
ncf<-nc_open(CFG_C3$forcing_file); macro<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
macro<-macro[as.Date(time)%in%ds]; mday<-macro[,.(Tmax_macro=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time))]
OBS<-merge(hb[,.(Tmax_obs=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)),by=id_plot]
readsc<-function(dir,scn){ fs<-list.files(file.path(dir,scn),pattern="\\.nc$",full.names=TRUE)
  rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
    r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
    s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time))]
    sd<-merge(s[,.(Tmax_sim=max(Tsim,na.rm=TRUE)),by=date],mday,by="date"); data.table(id_plot=id,ds=mean(sd$Tmax_sim-sd$Tmax_macro,na.rm=TRUE))}),fill=TRUE) }
getpp<-function(scn,dir=d3) merge(readsc(dir,scn),OBS,by="id_plot")

R2<-function(x,y) cor(x,y)^2
cat("=== (1) LADOPT effect: full-LAD vs top-d_opt LAD (SAME total LAI) ===\n")
pairs<-list(c("DYN_ALS","DYN_ALS_LADOPT"), c("STATIC_S2_ATBD","STATIC_S2_ATBD_LADOPT"),
            c("STATIC_S2_DOPT","STATIC_S2_DOPT_LADOPT"), c("DYN_S2_DOPT","DYN_S2_DOPT_LADOPT"))
sdo<-sd(OBS$do)
for(p in pairs){ a<-getpp(p[1]); b<-getpp(p[2])
  cat(sprintf("  %-22s R2=%.3f SDrec=%.2f  ->  %-24s R2=%.3f SDrec=%.2f   (dR2=%+.3f)\n",
    p[1],R2(a$ds,a$do),sd(a$ds)/sdo, p[2],R2(b$ds,b$do),sd(b$ds)/sdo, R2(b$ds,b$do)-R2(a$ds,a$do))) }
cat(sprintf("  [obs SD(ΔTmax)=%.2f °C]\n", sdo))

cat("\n=== (2) jackknife R² of DYN_ALS_LADOPT (drop-1 plot) ===\n")
b<-getpp("DYN_ALS_LADOPT"); full<-R2(b$ds,b$do)
jk<-sapply(seq_len(nrow(b)),function(i) R2(b$ds[-i],b$do[-i]))
infl<-order(abs(jk-full),decreasing=TRUE)[1:5]
cat(sprintf("  full R2=%.3f | drop-1 range [%.3f, %.3f]\n", full, min(jk), max(jk)))
cat("  5 most influential plots (id: R2_without):\n")
for(i in infl) cat(sprintf("    %-8s %.3f\n", b$id_plot[i], jk[i]))

cat("\n=== (3) ranking: DYN_ALS_LADOPT vs the DYN_S2 family (summer ΔTmax) ===\n")
fam<-c("DYN_ALS_LADOPT","DYN_ALS","STATIC_ALS","DYN_S2_ATBD","DYN_S2_DOPT","DYN_S2_DOPT_LADOPT","DYN_S2_RESCALED","STATIC_S2_ATBD")
out<-rbindlist(lapply(fam,function(s){d<-getpp(s);data.table(scenario=s,n=nrow(d),dT_R2=round(R2(d$ds,d$do),3),
  SDrec=round(sd(d$ds)/sdo,2),bias=round(mean(d$ds-d$do),3),RMSE=round(sqrt(mean((d$ds-d$do)^2)),3))}))
print(out[order(-dT_R2)])
cat("\nDONE\n")
