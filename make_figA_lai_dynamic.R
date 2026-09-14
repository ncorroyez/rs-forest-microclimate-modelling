# Annex: LAI per scenario per plot over summer (6 scenarios), S2 series rebuilt from
# NOT_MASKED (all 53 plots dynamic). LiDAR static.
suppressPackageStartupMessages({ library(terra); library(data.table); library(ggplot2); library(stringr); library(rmusica); library(musica.tools) })
src <- list.files("R","\\.R$",full.names=TRUE); src<-src[!grepl("/(h1_|lovb_)",src)]; invisible(lapply(src,source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
prep<-load_lai_prep(CFG_C3); df<-as.data.table(prep$df_plots)
df[, pid:=sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT:=LAI_ALS]
ropt<-mean(df$LAI_S2_DOPT,na.rm=TRUE)/mean(df$LAI_S2_ATBD,na.rm=TRUE); df[is.na(LAI_S2_DOPT), LAI_S2_DOPT:=LAI_S2_ATBD*ropt]
rA<-mean(df$LAI_ALS)/mean(df$LAI_S2_ATBD); bod<-mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_DOPT); ropt_p<-setNames(df$LAI_S2_DOPT/df$LAI_S2_ATBD, df$pid)
nmdir<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files<-list.files(nmdir, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
stk<-rast(files); dates<-as.Date(str_extract(basename(files),"\\d{4}-\\d{2}-\\d{2}"))
pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs=crs(stk)); vals<-as.data.frame(terra::extract(stk,pts))[,-1,drop=FALSE]
long<-rbindlist(lapply(seq_along(dates), function(i) data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
atbd<-smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
DOY<-152:273; rows<-list(); add<-function(id,scn,doy,lai) rows[[length(rows)+1]]<<-data.table(id_plot=id,scenario=scn,doy=doy,LAI=lai)
for(i in seq_len(nrow(df))){ p<-df$pid[i]; id<-as.character(df$id_plot[i]); a<-atbd[[p]]; if(is.null(a)) next
  A<-as.data.table(a)[doy%in%DOY]; optser<-A$lai*ropt_p[[p]]
  add(id,"S2 ATBD (dyn)",A$doy,A$lai); add(id,"S2 opt (dyn)",A$doy,optser)
  add(id,"S2 ATBD x ratio->full (dyn)",A$doy,A$lai*rA); add(id,"S2 opt x ratio->d_opt (dyn)",A$doy,optser*bod)
  add(id,"LiDAR full (static)",DOY,df$LAI_ALS[i]); add(id,"LiDAR d_opt (static)",DOY,df$LAI_ALS_DOPT[i])
}
m<-rbindlist(rows,fill=TRUE)
lev<-c("S2 ATBD (dyn)","S2 opt (dyn)","S2 ATBD x ratio->full (dyn)","S2 opt x ratio->d_opt (dyn)","LiDAR full (static)","LiDAR d_opt (static)")
m[, scenario:=factor(scenario,levels=lev)]
g<-ggplot(m,aes(doy,LAI,group=id_plot))+geom_line(alpha=0.16,colour="#1A9850",linewidth=0.3)+facet_wrap(~scenario,ncol=3)+
  labs(x="Day of year (summer 2021)",y="LAI (one-sided, m² m⁻²)")+theme_article(11)
outdir<-"/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"FigA_scenario_lai_timeseries"),g,9.5,4.6)
ggsave_article("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures/FigA_scenario_lai_timeseries",g,9.5,4.6)
cat("DONE annex (Not_Masked, plots:", length(unique(m$id_plot)),")\n")
