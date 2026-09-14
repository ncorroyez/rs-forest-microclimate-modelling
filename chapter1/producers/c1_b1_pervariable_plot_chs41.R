# Fig B1 (per-variable within-canopy vertical T gradient), CHS41-Rmerge/no-wind.
# RE-AUTHORED 2026-09-14: the shipped FigB1_pervariable_gradient_chs41.png was plotted
# interactively on 2026-09-08 and its script was lost (see PRODUCER_FIGURE_MAP_2026-09-14.md).
# The b3 sims (c1_b3_pervariable_chs41.R) write nc_b3_chs41/{P}_{base,LAIp,fCovp,Hmaxp}.nc;
# this reuses the vertical-profile extraction of c1_fig6_fig7_chs41.R (Fig 7) on those arms.
# Data are faithful; validate the aesthetic against the shipped png (md5 will not match).
# Reads : out_files/Chapter1/nc_b3_chs41/{P}_{base,LAIp,fCovp,Hmaxp}.nc
# Writes: out_files/Chapter1/figures/FigB1_pervariable_gradient_chs41.png  (hand-copy to chapter1/manuscript/figures/)
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(ggplot2)})
source("R/cluster_relabel.R")
FIG<-Sys.getenv("CH1_CHS41_FIGDIR","out_files/Chapter1/figures"); dir.create(FIG,showWarnings=FALSE,recursive=TRUE)
th<-theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position="top")

# daytime (10-16 h) JJAS-mean vertical air-temperature profile from one MuSICA output
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
prof<-function(f){nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  on.exit(nc_close(nc)); v<-names(nc$var)
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%v))return(NULL)
  tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC"); tt<-t0+dhours(ncvar_get(nc,"time"))
  Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height"); vh<-median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE)
  keep<-as.Date(tt)%in%ds & hour(tt)>=10 & hour(tt)<=16
  Tmean<-rowMeans(Tk[,keep,drop=FALSE],na.rm=TRUE)-273.15
  data.table(z_m=rh*vh, Tair=Tmean)}

# arms: base + the three single-variable native steps (real LAD throughout)
ARMS<-c(base="base", LAIp="+Leaf area (+0.5)", fCovp="+Cover (+10 pts)", Hmaxp="+Height (+1 m)")
grid<-rbindlist(lapply(paste0("P",1:4),function(p) rbindlist(lapply(names(ARMS),function(a){
  d<-prof(sprintf("out_files/Chapter1/nc_b3_chs41/%s_%s.nc",p,a))
  if(is.null(d))return(NULL); d[,`:=`(P=p,arm=ARMS[[a]])]}))),fill=TRUE)
grid[,P:=factor(P,levels=paste0("P",1:4))]
grid[,arm:=factor(arm,levels=unname(ARMS))]
PAL_ARM<-c("base"="grey35","+Leaf area (+0.5)"="#1B9E77","+Cover (+10 pts)"="#D95F02","+Height (+1 m)"="#7570B3")

g<-ggplot(grid,aes(Tair,z_m,colour=arm))+
  geom_hline(yintercept=1,linetype="dotted",colour="grey55")+   # 1 m logger readout
  geom_path(linewidth=.7)+
  scale_colour_manual(values=PAL_ARM,name=NULL)+
  facet_wrap(~P,scales="free",nrow=1)+
  labs(x="Air temperature (°C), daytime JJAS mean", y="Height (m)")+th
ggsave(file.path(FIG,"FigB1_pervariable_gradient_chs41.png"),g,width=13,height=4,dpi=300,bg="white")

# report the near-1 m readout per arm, by archetype
cat("=== Tair at ~1 m per arm, by archetype ===\n")
g1<-grid[,.SD[which.min(abs(z_m-1))],by=.(P,arm)][,.(P,arm,z1=round(z_m,2),Tair_1m=round(Tair,3))]
print(g1)
cat("DONE -> FigB1_pervariable_gradient_chs41.png\n")
