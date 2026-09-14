# Fig 6 (H2 controlled test) + Fig 7 (within-canopy vertical T gradient, real vs uniform),
# CHS41-Rmerge/no-wind. Fig6 from h2_controlled_chs41.csv; Fig7 from nc_archetype_chs41/.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(ggplot2);library(patchwork)})
source("R/cluster_relabel.R")
FIG<-"out_files/Chapter1/figures"; dir.create(FIG,showWarnings=FALSE,recursive=TRUE)
th<-theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position="top")

# ---- Fig 6: H2 controlled (ΔTmax vs centre of mass, by LAI) ----
H<-fread("out_files/Chapter1/tables/h2_controlled_chs41.csv")
H[,LAI:=factor(LAI)]
g6<-ggplot(H,aes(com,dtmax,colour=LAI,group=LAI))+
  geom_hline(yintercept=0,linetype="dotted",colour="grey50")+
  geom_line(linewidth=.6)+geom_point(size=2)+
  scale_colour_viridis_d(option="C",end=.9,name="one-sided LAI")+
  labs(x="LAD centre of mass (fraction of Hmax; 0 bottom-heavy, 1 top-heavy)",
       y=expression(Delta*italic(T)[max]*" (°C)"))+th
ggsave(file.path(FIG,"Fig6_h2_controlled_chs41.png"),g6,width=7,height=5,dpi=300,bg="white")

# ---- Fig 7: vertical T profile per archetype, real vs uniform, daytime JJAS mean ----
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
prof<-function(f){nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  on.exit(nc_close(nc)); v<-names(nc$var)
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%v))return(NULL)
  tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC"); tt<-t0+dhours(ncvar_get(nc,"time"))
  Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height"); vh<-median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE)
  keep<-as.Date(tt)%in%ds & hour(tt)>=10 & hour(tt)<=16
  Tmean<-rowMeans(Tk[,keep,drop=FALSE],na.rm=TRUE)-273.15
  data.table(z_m=rh*vh, Tair=Tmean)}
arch<-rbindlist(lapply(paste0("P",1:4),function(p) rbindlist(lapply(c("real","unif"),function(a){
  d<-prof(sprintf("out_files/Chapter1/nc_archetype_chs41/%s_%s.nc",p,a))
  if(is.null(d))return(NULL); d[,`:=`(P=p,arm=ifelse(a=="real","real","uniform"))]}))),fill=TRUE)
arch[,P:=factor(P,levels=paste0("P",1:4))]
g7<-ggplot(arch,aes(Tair,z_m,colour=P,linetype=arm))+
  geom_hline(yintercept=1,linetype="dotted",colour="grey55")+   # 1 m logger readout
  geom_path(linewidth=.7)+
  scale_colour_manual(values=PAL_CLUSTER,name=NULL)+
  scale_linetype_manual(values=c(real=1,uniform=2),name=NULL)+
  facet_wrap(~P,scales="free",nrow=1)+
  labs(x="Air temperature (°C), daytime JJAS mean", y="Height (m)")+th
ggsave(file.path(FIG,"Fig7_vertical_Tprofile_chs41.png"),g7,width=12,height=4,dpi=300,bg="white")

# do the real and uniform profiles converge at 1 m? report the gap at ~1 m
cat("=== real vs uniform Tair gap at ~1 m, by archetype ===\n")
g1m<-arch[,.(z=z_m[which.min(abs(z_m-1))]),by=.(P,arm)]
w<-dcast(merge(arch,g1m,by=c("P","arm"))[abs(z_m-z)<1e-6],P~arm,value.var="Tair")
w[,gap_1m:=round(real-uniform,3)]; print(w)
cat("DONE -> Fig6, Fig7\n")
