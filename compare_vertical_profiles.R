# ==============================================================================
# Vertical in-canopy profiles of T, wind, RH, VPD, compared BETWEEN archetypes
# P1-P4 (open->dense) AND between MuSICA v3.2.0 vs v3.2.3-iter. Midday (12-14h
# UTC) mean over JJAS 2021 — the condition with the strongest canopy gradients.
#
# nc gives Tair_z (K), wind_z (m/s), wair_z (water-vapour mixing ratio mol/mol)
# on 15 levels at relative heights z/Hmax (up to 1.77, i.e. above canopy top=1).
#   e   = (r/(1+r)) * P         [Pa]   (P = 101325 Pa, lowland, fixed)
#   es  = 611.2*exp(17.62*Tc/(243.12+Tc))   (Tetens)
#   RH  = 100*e/es ;  VPD = (es-e)/1000  [kPa]
# Plotted vs relative height (1.0 = canopy top) so gradient shape is comparable
# across archetypes of different stature.
#   Rscript compare_vertical_profiles.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(ggplot2)
  library(dplyr); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
YEAR<-2021L; MM<-6:9; MIDDAY<-12:14; P0<-101325
OUT_TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
OUT_FIG<-"/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
NCDIRS<-c("v3.2.0"=file.path(CFG_C3$out_dir,"nc","STATIC_ALS"),
          "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter","STATIC_ALS"))

# ---- 4 archetype plots (same selection as the time-series figures) -----------
prep<-load_lai_prep(CFG_C3); df<-as.data.table(prep$df_plots)
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
ho<-unique(as.character(hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)]$id_plot))
cand<-df[id_plot%in%ho&!(id_plot%in%CFG_C3$ids_to_remove)][order(LAI_ALS)]
sel<-cand[round(quantile(seq_len(nrow(cand)),c(.10,.40,.65,.90)))]
sel[,P:=sprintf("P%d",.I)]; sel[,lab:=sprintf("%s (LAI %.1f, Hmax %.0fm)",P,LAI_ALS,Hmax)]
cat("Archetypes:\n"); print(sel[,.(P,id_plot,LAI_ALS=round(LAI_ALS,1),Hmax=round(Hmax,1))])

# ---- time index for the nc (build once) --------------------------------------
prof_one<-function(f){ nc<-nc_open(f)
  tt<-as.POSIXct(force_utc_nc(f,"time")[1],tz="UTC")  # not used; build from forcing-like
  rh<-as.numeric(ncvar_get(nc,"relative_height")); hmax<-as.numeric(ncvar_get(nc,"veget_height_top"))[1]
  # time vector from nc time dim via musica.tools helper
  tvec<-force_utc_nc(f,"time"); keep<-which(year(tvec)==YEAR & month(tvec)%in%MM & hour(tvec)%in%MIDDAY)
  Ta<-ncvar_get(nc,"Tair_z")[,keep]; Wa<-ncvar_get(nc,"wair_z")[,keep]; Wi<-ncvar_get(nc,"wind_z")[,keep]; nc_close(nc)
  Tc<-rowMeans(Ta,na.rm=TRUE)-273.15; r<-rowMeans(Wa,na.rm=TRUE); wind<-rowMeans(Wi,na.rm=TRUE)
  e<-(r/(1+r))*P0; es<-611.2*exp(17.62*Tc/(243.12+Tc))
  data.table(relh=rh, height=rh*hmax, Tair=Tc, wind=wind, RH=100*e/es, VPD=(es-e)/1000) }

PROF<-rbindlist(lapply(seq_len(nrow(sel)),function(k){ rbindlist(lapply(names(NCDIRS),function(v){
  f<-file.path(NCDIRS[v],sprintf("musica_out_HOBO_%s.nc",sel$id_plot[k])); if(!file.exists(f))return(NULL)
  p<-prof_one(f); p[,`:=`(P=sel$P[k],lab=sel$lab[k],version=v)] })) }))
fwrite(PROF, file.path(OUT_TAB,"vertical_profiles_midday.csv"))

# ---- gradient table: understory (lowest level ~1-2 m) minus canopy-top (relh~1)
grad<-PROF[, { ub<-.SD[which.min(abs(relh-0.03))]; tp<-.SD[which.min(abs(relh-1.0))]
  .(Tair_grad=ub$Tair-tp$Tair, wind_grad=ub$wind-tp$wind, RH_grad=ub$RH-tp$RH, VPD_grad=ub$VPD-tp$VPD,
    understory_VPD=ub$VPD, top_VPD=tp$VPD) }, by=.(P,version)]
cat("\n=== Vertical gradient (understory - canopy top), midday JJAS ===\n")
print(grad[order(P,version),.(P,version,dTair=round(Tair_grad,2),dWind=round(wind_grad,2),
  dRH=round(RH_grad,1),dVPD=round(VPD_grad,2))])
fwrite(grad, file.path(OUT_TAB,"vertical_profiles_gradient.csv"))

# ---- figure: 4 variables x relative height, colour=archetype, linetype=version
PL<-melt(PROF, id.vars=c("relh","height","P","lab","version"),
         measure.vars=c("Tair","wind","RH","VPD"), variable.name="var")
PL[,var:=factor(var,levels=c("Tair","wind","RH","VPD"),
                labels=c("Air temperature (degC)","Wind speed (m/s)","Relative humidity (%)","VPD (kPa)"))]
g<-ggplot(PL[relh<=1.5], aes(value, relh, colour=P, linetype=version))+
  geom_hline(yintercept=1, colour="grey75", linewidth=.3)+
  geom_path(linewidth=.6)+ geom_point(size=.8, alpha=.6)+
  facet_wrap(~var, nrow=1, scales="free_x")+
  scale_colour_manual(values=PAL_CLUSTER, name="Archetype")+
  scale_linetype_manual(values=c("v3.2.0"="solid","v3.2.3 iter"="22"), name="Binary")+
  labs(x="Value", y="Relative height z / Hmax  (1.0 = canopy top)",
       subtitle=paste0("Midday (12-14h) mean vertical profiles, JJAS 2021. Colour = archetype P1 open -> P4 dense; line = MuSICA binary.",
                       "\nGrey line = canopy top. ",paste(sel$lab,collapse="  |  "))) +
  theme_article(10)+theme(legend.position="right")
ggsave_article(file.path(OUT_FIG,"Fig_vertical_profiles_midday"), g, 13, 5.2)
cat("\nDONE -> Fig_vertical_profiles_midday + vertical_profiles_{midday,gradient}.csv\n")
