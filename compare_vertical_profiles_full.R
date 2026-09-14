# ==============================================================================
# Vertical in-canopy profiles (T, wind, RH, VPD), archetypes P1-P4 x binaries
# v3.2.0/v3.2.3, across time-of-day {day 12-14h, night 22-03h, global 24h} and
# period {JJAS, hot10 = macro daily Tmax>=p90}. One comprehensive figure
# facet_grid(period+tod ~ variable) + a full gradient table (understory-top).
#   RH/VPD from wair_z (mixing ratio mol/mol) and Tair_z via Tetens, P=101325 Pa.
#   Rscript compare_vertical_profiles_full.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(ggplot2)
  library(dplyr); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
YEAR<-2021L; MM<-6:9; P0<-101325
TOD<-list(day=12:14, night=c(22,23,0,1,2,3), global=0:23)
OUT_TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
OUT_FIG<-"/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
NCDIRS<-c("v3.2.0"=file.path(CFG_C3$out_dir,"nc","STATIC_ALS"),
          "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter","STATIC_ALS"))

# hot days (macro daily Tmax >= p90, JJAS)
ncf<-nc_open(CFG_C3$forcing_file)
mh<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
mh<-mh[year(time)==YEAR&month(time)%in%MM]; mday<-mh[,.(mx=max(Tm)),by=.(d=as.Date(time))]
hotd<-mday[mx>=quantile(mx,.90),d]; PER<-list(JJAS="all", hot10=hotd)

prep<-load_lai_prep(CFG_C3); df<-as.data.table(prep$df_plots)
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
ho<-unique(as.character(hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)]$id_plot))
cand<-df[id_plot%in%ho&!(id_plot%in%CFG_C3$ids_to_remove)][order(LAI_ALS)]
sel<-cand[round(quantile(seq_len(nrow(cand)),c(.10,.40,.65,.90)))]; sel[,P:=sprintf("P%d",.I)]
sel[,lab:=sprintf("%s (LAI %.1f, Hmax %.0fm)",P,LAI_ALS,Hmax)]
cat("Archetypes:\n"); print(sel[,.(P,id_plot,LAI=round(LAI_ALS,1),Hmax=round(Hmax))])

# read one nc fully, return profiles for all tod x period
prof_all<-function(f){ nc<-nc_open(f); rh<-as.numeric(ncvar_get(nc,"relative_height"))
  hmax<-as.numeric(ncvar_get(nc,"veget_height_top"))[1]; tv<-force_utc_nc(f,"time")
  Ta<-ncvar_get(nc,"Tair_z"); Wa<-ncvar_get(nc,"wair_z"); Wi<-ncvar_get(nc,"wind_z"); nc_close(nc)
  base<-year(tv)==YEAR & month(tv)%in%MM
  out<-list()
  for(pn in names(PER)){ pdays<-PER[[pn]]; inper<-if(identical(pdays,"all")) base else base & (as.Date(tv)%in%pdays)
    for(tn in names(TOD)){ k<-which(inper & hour(tv)%in%TOD[[tn]]); if(!length(k))next
      Tc<-rowMeans(Ta[,k,drop=FALSE])-273.15; r<-rowMeans(Wa[,k,drop=FALSE]); wind<-rowMeans(Wi[,k,drop=FALSE])
      e<-(r/(1+r))*P0; es<-611.2*exp(17.62*Tc/(243.12+Tc))
      out[[paste(pn,tn)]]<-data.table(period=pn,tod=tn,relh=rh,height=rh*hmax,
                                      Tair=Tc,wind=wind,RH=100*e/es,VPD=(es-e)/1000) } }
  rbindlist(out) }

PROF<-rbindlist(lapply(seq_len(nrow(sel)),function(k) rbindlist(lapply(names(NCDIRS),function(v){
  f<-file.path(NCDIRS[v],sprintf("musica_out_HOBO_%s.nc",sel$id_plot[k])); if(!file.exists(f))return(NULL)
  p<-prof_all(f); p[,`:=`(P=sel$P[k],version=v)] }))))
PROF[,tod:=factor(tod,levels=c("day","night","global"))]
PROF[,period:=factor(period,levels=c("JJAS","hot10"))]
fwrite(PROF, file.path(OUT_TAB,"vertical_profiles_full.csv"))

# ---- gradient table: understory (relh~0.03) - canopy top (relh~1.0) ----------
grad<-PROF[,{ ub<-.SD[which.min(abs(relh-0.03))]; tp<-.SD[which.min(abs(relh-1.0))]
  .(dTair=ub$Tair-tp$Tair, dWind=ub$wind-tp$wind, dRH=ub$RH-tp$RH, dVPD=ub$VPD-tp$VPD) }, by=.(P,version,period,tod)]
cat("\n=== Vertical gradient (understory - canopy top) — all conditions ===\n")
print(grad[order(period,tod,P,version),.(period,tod,P,version,
  dTair=round(dTair,2),dWind=round(dWind,2),dRH=round(dRH,1),dVPD=round(dVPD,2))])
fwrite(grad, file.path(OUT_TAB,"vertical_profiles_full_gradient.csv"))

# ---- comprehensive figure: facet_grid(period+tod ~ variable) -----------------
PL<-melt(PROF[relh<=1.5], id.vars=c("relh","P","version","period","tod"),
         measure.vars=c("Tair","wind","RH","VPD"), variable.name="var")
PL[,var:=factor(var,levels=c("Tair","wind","RH","VPD"),
                labels=c("T (degC)","Wind (m/s)","RH (%)","VPD (kPa)"))]
PL[,rowlab:=factor(paste(period,tod), levels=c("JJAS day","JJAS night","JJAS global",
                                               "hot10 day","hot10 night","hot10 global"))]
g<-ggplot(PL, aes(value, relh, colour=P, linetype=version))+
  geom_hline(yintercept=1, colour="grey75", linewidth=.3)+ geom_path(linewidth=.55)+
  facet_grid(rowlab ~ var, scales="free_x")+
  scale_colour_manual(values=PAL_CLUSTER, name="Archetype")+
  scale_linetype_manual(values=c("v3.2.0"="solid","v3.2.3 iter"="22"), name="Binary")+
  labs(x="Value", y="Relative height z / Hmax  (1.0 = canopy top)",
       subtitle=paste0("Vertical profiles by time-of-day x period. Colour = archetype P1 open -> P4 dense; line = binary. ",
                       paste(sel$lab,collapse="  |  "))) +
  theme_article(9)+theme(legend.position="right")
ggsave_article(file.path(OUT_FIG,"Fig_vertical_profiles_full"), g, 12, 12)
cat("\nDONE -> Fig_vertical_profiles_full + vertical_profiles_full{,_gradient}.csv\n")
