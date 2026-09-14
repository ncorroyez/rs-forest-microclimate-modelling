# Fig 6 v2 (prototype): 2D response surface of ΔTmax over the Gaussian profile
# parameters — peak height mu and vertical diffuseness sigma — faceted by LAI.
# Overlays the realised sigma_rel band of real cLHS profiles (analysis a).
# Reads out_files/Chapter1/tables/h2_gaussian_grid_chs41.csv.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("R/cluster_relabel.R")
R<-fread("out_files/Chapter1/tables/h2_gaussian_grid_chs41.csv")
R[,LAIf:=factor(sprintf("one-sided LAI %g",LAI),levels=sprintf("one-sided LAI %g",sort(unique(LAI))))]
# realised sigma_rel band from real profiles (same metric as analysis a)
s<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
s<-s[is.finite(LAI)&is.finite(Hmax)&Hmax>0]; ladc<-grep("LAD_Layer_",names(s),value=TRUE); zc<-as.numeric(gsub("LAD_Layer_","",ladc))
M<-as.matrix(s[,..ladc]);M[is.na(M)]<-0; sig<-numeric(nrow(s))
for(i in seq_len(nrow(s))){d<-M[i,];zr<-zc/s$Hmax[i];sw<-sum(d);if(sw<=0){sig[i]<-NA;next};mu<-sum(d*zr)/sw;sig[i]<-sqrt(sum(d*(zr-mu)^2)/sw)}
band<-quantile(sig,c(.25,.75),na.rm=TRUE)
lim<-max(abs(R$dtmax),na.rm=TRUE)
g<-ggplot(R,aes(mu,sigma,fill=dtmax))+
  geom_raster(interpolate=TRUE)+
  geom_hline(yintercept=band,linetype="dashed",colour="grey20",linewidth=.3)+
  geom_text(aes(label=sprintf("%+.2f",dtmax)),size=2.4,colour="grey15")+
  scale_fill_gradient2(low="#1A9850",mid="white",high="#D7191C",midpoint=0,limits=c(-lim,lim),
    name=expression(Delta*italic(T)[max]*" (°C)"))+
  facet_wrap(~LAIf,nrow=1)+
  labs(x=expression("peak height "*mu*" (fraction of "*italic(H)[max]*")"),
       y=expression("vertical diffuseness "*sigma*" (fraction of "*italic(H)[max]*")"))+
  theme_bw(base_size=11)+theme(panel.grid=element_blank(),legend.position="right",
    strip.text=element_text(face="bold"))
ggsave("out_files/Chapter1/figures/Fig6_h2_gaussian_chs41.png",g,width=13.6,height=3.5,dpi=300,bg="white")
cat("realised sigma_rel IQR band:",round(band,3),"\n")
cat("=== sigma-slope vs mu-slope per LAI (mean |ΔTmax| change per full axis span) ===\n")
for(l in sort(unique(R$LAI))){d<-R[LAI==l]
  se<-d[,max(dtmax)-min(dtmax),by=mu]$V1; me<-d[,max(dtmax)-min(dtmax),by=sigma]$V1
  cat(sprintf("LAI %g: sigma-effect(mean over mu) %.3f | mu-effect(mean over sigma) %.3f °C\n",l,mean(se),mean(me)))}
cat("saved Fig6v2_gaussian_musigma_chs41.png\n")
