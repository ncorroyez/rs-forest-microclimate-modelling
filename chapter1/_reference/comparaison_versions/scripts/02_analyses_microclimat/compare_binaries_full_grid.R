# ==============================================================================
# FULL comparison grid: MuSICA v3.2.0 vs v3.2.3-iter (STATIC_ALS) vs observed.
# Metrics x time-of-day x period, both RANKING and AMPLITUDE, plus bias/RMSE and
# Williams test (v3.2.0 vs v3.2.3, shared=obs). JJAS 2021, n=53 plots.
#
#   Families : slope = lm(T_micro ~ T_macro) coefficient (Gril coupling slope)
#              offset = mean(T_micro - T_macro)          (buffering magnitude)
#   Time-of-day : all (24h) | day (10-15h) | night (22-03h)
#   Period      : JJAS (Jun-Sep) | hot10 (macro daily Tmax >= p90, ~10% days)
#   Extras (global) : DTR = mean(Tmax-Tmin) ; equilibrium = intercept/(1-slope)
#   Per cell, sim vs obs across 53 plots:
#     ranking  : Pearson r, Spearman rho
#     amplitude: SD-recovery = SD(pred)/SD(obs)
#     error    : bias = mean(pred-obs), RMSE
#     Williams : t,p for r(obs,v3.2.0) vs r(obs,v3.2.3)
#   Rscript compare_binaries_full_grid.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(ggplot2)
  library(dplyr); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
YEAR<-2021L; MM<-6:9
TOD <- list(all=0:23, day=10:15, night=c(22,23,0,1,2,3))
OUT_TAB<-"/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/tables"
OUT_FIG<-"/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/figures"
NCDIRS<-c("v3.2.0"=file.path(CFG_C3$out_dir,"nc","STATIC_ALS"),
          "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter","STATIC_ALS"))
slp<-function(y,x){k<-is.finite(x)&is.finite(y); if(sum(k)<24)return(NA_real_); as.numeric(coef(lm(y[k]~x[k]))[2])}
intc<-function(y,x){k<-is.finite(x)&is.finite(y); if(sum(k)<24)return(NA_real_); as.numeric(coef(lm(y[k]~x[k]))[1])}

# macro hourly + hot days
ncf<-nc_open(CFG_C3$forcing_file)
mh<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),T_macro=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
mh<-mh[year(time)==YEAR&month(time)%in%MM][,`:=`(hr=floor_date(time,"hour"),H=hour(time),d=as.Date(time))]
mday<-mh[,.(Tmac_max=max(T_macro),Tmac_min=min(T_macro)),by=d]
hotd<-mday[Tmac_max>=quantile(Tmac_max,.90),d]
PER<-list(JJAS=unique(mh$d), hot10=hotd)
cat(sprintf("JJAS days=%d | hot10 days=%d\n", length(PER$JJAS), length(PER$hot10)))

prep<-load_lai_prep(CFG_C3); df<-as.data.table(prep$df_plots)
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&year(time)==YEAR&month(time)%in%MM,
       .(id_plot=as.character(id_plot),T=t_hobo,hr=floor_date(time,"hour"))]
hb<-merge(hb,mh[,.(hr,T_macro,H,d)],by="hr")

# per-plot, per-source long metric table
cells<-function(r){ out<-list()
  for(pn in names(PER)){ rp<-r[d%in%PER[[pn]]]
    for(tn in names(TOD)){ rt<-rp[H%in%TOD[[tn]]]
      out[[paste0("slope.",tn,".",pn)]]<-slp(rt$T,rt$T_macro)
      out[[paste0("offset.",tn,".",pn)]]<-mean(rt$T-rt$T_macro,na.rm=TRUE) }
    rd<-rp[,.(Tmax=max(T),Tmin=min(T)),by=d]; rd<-merge(rd,mday,by="d")
    out[[paste0("DTR.all.",pn)]]<-mean(rd$Tmax-rd$Tmin,na.rm=TRUE)
    ra<-rp  # equilibrium from all-hours lm
    out[[paste0("equil.all.",pn)]]<-{s<-slp(ra$T,ra$T_macro); i<-intc(ra$T,ra$T_macro); i/(1-s)} }
  as.data.table(out) }

P<-rbindlist(lapply(df$id_plot,function(id){
  ro<-hb[id_plot==id,.(T,T_macro,H,d)]; rows<-cbind(data.table(id_plot=id,src="obs"),cells(ro))
  for(v in names(NCDIRS)){ f<-file.path(NCDIRS[v],sprintf("musica_out_HOBO_%s.nc",id)); if(!file.exists(f))next
    s<-as.data.table(get_tair_at_z(nc_open(f),1.0));s<-s[year(time)==YEAR&month(time)%in%MM][,hr:=floor_date(time,"hour")]
    s<-merge(s,mh[,.(hr,T_macro,H,d)],by="hr"); rr<-s[,.(T=Tair_sim,T_macro,H,d)]
    rows<-rbind(rows, cbind(data.table(id_plot=id,src=ifelse(v=="v3.2.0","v320","v323")),cells(rr))) }
  rows }), fill=TRUE)
fwrite(P, file.path(OUT_TAB,"compare_full_grid_plotdata.csv"))

# ---- aggregate: per metric cell, sim vs obs ----------------------------------
mcols<-setdiff(names(P),c("id_plot","src"))
W<-dcast(melt(P,id.vars=c("id_plot","src"),measure.vars=mcols),
         id_plot+variable~src, value.var="value")
williams<-function(rxy,rxz,ryz,n){Rd<-1-rxy^2-rxz^2-ryz^2+2*rxy*rxz*ryz
  t<-(rxy-rxz)*sqrt((n-1)*(1+ryz)/(2*((n-1)/(n-3))*Rd+((rxy+rxz)^2/4)*(1-ryz)^3)); c(t=t,p=2*pt(-abs(t),n-3))}
agg<-W[,{ o<-obs; f<-is.finite(o)
  res<-list()
  for(m in c("v320","v323")){ p<-get(m); k<-f&is.finite(p)
    res[[paste0("r_",m)]]<-cor(o[k],p[k]); res[[paste0("rho_",m)]]<-cor(o[k],p[k],method="spearman")
    res[[paste0("SDrec_",m)]]<-sd(p[k])/sd(o[k]); res[[paste0("bias_",m)]]<-mean(p[k]-o[k]); res[[paste0("rmse_",m)]]<-sqrt(mean((p[k]-o[k])^2)) }
  kk<-f&is.finite(v320)&is.finite(v323); wl<-williams(cor(o[kk],v320[kk]),cor(o[kk],v323[kk]),cor(v320[kk],v323[kk]),sum(kk))
  c(res, list(will_t=wl["t"], will_p=wl["p"])) }, by=variable]
agg[,c("metric","tod","period"):=tstrsplit(variable,".",fixed=TRUE)]
setcolorder(agg,c("metric","tod","period"))
fwrite(agg, file.path(OUT_TAB,"compare_full_grid_stats.csv"))

cat("\n=== RANKING (Spearman rho) ===\n")
print(agg[order(metric,tod,period),.(metric,tod,period,rho_v320=round(rho_v320,2),rho_v323=round(rho_v323,2),
  better=fifelse(rho_v323>rho_v320,"v3.2.3","v3.2.0"), will_p=round(will_p,3))])
cat("\n=== AMPLITUDE (SD recovery) + bias ===\n")
print(agg[order(metric,tod,period),.(metric,tod,period,SDrec_v320=round(SDrec_v320,2),SDrec_v323=round(SDrec_v323,2),
  bias_v320=round(bias_v320,2),bias_v323=round(bias_v323,2))])

# ---- figure: two heatmaps (ranking rho ; amplitude SDrec) --------------------
L<-melt(agg, id.vars=c("metric","tod","period"),
        measure.vars=c("rho_v320","rho_v323","SDrec_v320","SDrec_v323"))
L[,c("stat","model"):=tstrsplit(as.character(variable),"_")]
L[,model:=factor(model,labels=c("v3.2.0","v3.2.3"))]
L[,row:=factor(paste(metric,tod),levels=unique(paste(metric,tod)))]
L[,col:=factor(paste(model,period),levels=c("v3.2.0 JJAS","v3.2.0 hot10","v3.2.3 JJAS","v3.2.3 hot10"))]
mk<-function(st,ttl,lim,mid){ ggplot(L[stat==st],aes(col,row,fill=value))+geom_tile(colour="white")+
  geom_text(aes(label=sprintf("%.2f",value)),size=2.6)+
  scale_fill_gradient2(low="#B2182B",mid="white",high="#2166AC",midpoint=mid,limits=lim,oob=scales::squish,name=NULL)+
  labs(x=NULL,y=NULL,subtitle=ttl)+theme_article(9)+
  theme(axis.text.x=element_text(angle=20,hjust=1),legend.position="right",panel.grid=element_blank()) }
g<-mk("rho","(a) Ranking — Spearman rho (sim vs obs)",c(-1,1),0) +
   patchwork::plot_spacer() +
   mk("SDrec","(b) Amplitude — SD recovery = SD(pred)/SD(obs)  (1 = perfect)",c(0,1.3),1) +
   patchwork::plot_layout(widths=c(1,0.04,1))
ggsave_article(file.path(OUT_FIG,"Fig_compare_full_grid"), g, 14, 6.5)
cat("\nDONE -> Fig_compare_full_grid + compare_full_grid_{plotdata,stats}.csv\n")
