# ==============================================================================
# Fig 5 (v3.2.3 iter + station): forward-inclusion cumulative-r curve vs 53 HOBO.
# GREEDY forward per cluster: at each step add the trait that most improves r
# against the observed ΔTmax (station-referenced). bits = LAI,Hmax,fCover,LAD.
# Reads coal_metrics_v323station_iter.rds (id_plot,bit,Delta_sim).
#   Rscript c1_forward_v323_curve.R
# Out: Fig_forward_curve_v323iter_station.{png,pdf} + tab
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(lubridate); library(data.table); library(ggplot2)
  src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
  source("pipeline/00_config.R"); source("scripts/_article_style.R") })
coal<-as.data.table(readRDS("out_files/Chapter1/coal_metrics_v323station_iter.rds"))   # id_plot,bit,Delta_sim
clu<-fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[,.(id_plot,P)]
ds<-CFG$date_seq; ds<-ds[format(ds,"%m")%in%c("06","07","08","09")]

# observed ΔTmax vs STATION macro (per logger)
nc<-nc_open("out_files/MuSICA_in_CHS41-Blois_2021-station_Rmerge.nc")
tu<-ncatt_get(nc,"time","units")$value; th<-ncvar_get(nc,"time"); t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
macD<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15); nc_close(nc)
macD<-macD[as.Date(time)%in%ds][,.(Tmax_mac=max(Tmac)),by=.(date=as.Date(time))]
hobo<-as.data.table(read.csv(CFG$hobo_temp_csv)); hobo[,datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hobo<-hobo[position_sensor=="a" & as.Date(datetime)%in%ds & !id_plot%in%CFG$ids_to_remove]
obs<-hobo[,.(Tmax_mic=max(t_hobo,na.rm=TRUE)),by=.(id_plot,date=as.Date(datetime))]
obs<-merge(obs,macD,by="date")[,.(Delta_obs=mean(Tmax_mic-Tmax_mac,na.rm=TRUE)),by=id_plot]

pos<-c(LAI=1L,Hmax=2L,fCover=3L,LAD=4L); mk<-function(a){b<-rep("0",4L);b[a]<-"1";paste(b,collapse="")}
rr<-function(bt,ids){ d<-merge(coal[bit==bt & id_plot%in%ids],obs,by="id_plot")
  ok<-is.finite(d$Delta_sim)&is.finite(d$Delta_obs); if(sum(ok)<3) return(NA_real_); cor(d$Delta_sim[ok],d$Delta_obs[ok]) }

greedy<-function(rowname,ids){
  active<-integer(0); rem<-1:4; out<-list(list(step=0L,add="Baseline",bit="0000",r=rr("0000",ids)))
  for(k in 1:4){ cand<-sapply(rem,function(a) rr(mk(c(active,a)),ids))
    best<-rem[which.max(cand)]; active<-c(active,best); rem<-setdiff(rem,best)
    tr<-names(pos)[pos==best]
    out[[k+1L]]<-list(step=k, add=if(k==4)sprintf("+%s (REF)",tr)else sprintf("+%s",tr), bit=mk(active), r=max(cand,na.rm=TRUE)) }
  d<-rbindlist(out); d[,rowlab:=rowname]; d }

rows<-rbindlist(c(lapply(paste0("P",1:4),function(cl) greedy(cl,clu[P==cl,id_plot])),
                  list(greedy("All",clu$id_plot))))
rows[,rowlab:=factor(rowlab,levels=c(paste0("P",1:4),"All"))]
cat("=== greedy forward order + cumulative r (v3.2.3 iter + station) ===\n")
print(rows[,.(rowlab,step,add,r=round(r,3))])
fwrite(rows,"out_files/Chapter1/tables/tab_forward_v323iter_station.csv")

pal<-c(P1=PAL_CLUSTER[["P1"]],P2=PAL_CLUSTER[["P2"]],P3=PAL_CLUSTER[["P3"]],P4=PAL_CLUSTER[["P4"]],All="grey30")
p<-ggplot(rows,aes(step,r,colour=rowlab,group=rowlab))+
  geom_line(linewidth=0.9)+geom_point(size=2.4)+
  ggrepel::geom_text_repel(aes(label=sub("^\\+","",sub(" \\(REF\\)","",add))),size=3,show.legend=FALSE,segment.color="grey80",max.overlaps=20)+
  scale_colour_manual(values=pal,name=NULL)+
  scale_x_continuous(breaks=0:4,labels=c("Baseline","+1","+2","+3","REF (+4)"))+
  labs(x="Traits added (greedy, most-improving first)", y=expression("Cumulative validation "*italic(r)*" vs 53 HOBO"),
       subtitle="Forward inclusion, v3.2.3 iter + station forcing. Which traits recover the observed ΔTmax pattern.")+
  theme_article(12)+theme(legend.position="right")
ggsave_article("out_files/Chapter1/figures/Fig_forward_curve_v323iter_station",p,8,5.2)
cat("DONE -> Fig_forward_curve_v323iter_station\n")
