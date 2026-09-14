# ==============================================================================
# Chapter 3 — Fig 5 (premise), full version: buffering slope ~ total LiDAR LAI for
# BOTH the observed HOBO loggers and the MuSICA STATIC_ALS simulation. Resolves the
# sign question: cor(slope, LAI) should be negative for both (more leaf area → lower
# slope → more buffering); the "MuSICA r = 0.92" in the text is the sim-vs-obs slope
# agreement, reported separately. Styled to _article_style.R, written to the article
# figures folder. Run from z_Example root:  Rscript fig5_premise_full.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(dplyr); library(lubridate); library(stringr)
  library(data.table); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
prep <- load_lai_prep(CFG_C3)
ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
era5 <- as.data.table(build_era5_hourly(CFG_C3$forcing_file, ds))   # time, Tair_era5
era5[, hr := floor_date(time,"hour")]
slp <- function(y,x){ v<-var(x,na.rm=TRUE); if(is.na(v)||v==0) NA else cov(y,x,use="complete.obs")/v }

# observed
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=floor_date(as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC"),"hour")]
hb <- hb[position_sensor=="a" & as.Date(time) %in% ds & !(id_plot %in% CFG_C3$ids_to_remove),
         .(t_hobo=mean(t_hobo,na.rm=TRUE)), by=.(id_plot=as.character(id_plot), hr=time)]
hb <- merge(hb, era5[,.(hr,Tair_era5)], by="hr")
obs <- hb[, .(slope=slp(t_hobo,Tair_era5), src="Observed (HOBO)"), by=id_plot]

# MuSICA STATIC_ALS (1 m)
ncp <- file.path(CFG_C3$out_dir, "nc", "STATIC_ALS")
fs <- list.files(ncp, pattern="\\.nc$", full.names=TRUE)
sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)")
  nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(r)||!nrow(r))return(NULL)
  d<-as.data.table(r)[as.Date(time) %in% ds, .(hr=floor_date(time,"hour"), Tsim=Tair_sim)]; d[,id_plot:=id][] }), fill=TRUE)
sim <- merge(sim, era5[,.(hr,Tair_era5)], by="hr")
mus <- sim[, .(slope=slp(Tsim,Tair_era5), src="MuSICA (STATIC_ALS)"), by=id_plot]

lai <- as.data.table(prep$df_plots)[, .(id_plot=as.character(id_plot), LAI_ALS)]
B <- merge(rbindlist(list(obs,mus)), lai, by="id_plot")
B <- B[is.finite(slope) & is.finite(LAI_ALS)]

r_obs <- B[src=="Observed (HOBO)", cor(slope, LAI_ALS)]
r_mus <- B[src=="MuSICA (STATIC_ALS)", cor(slope, LAI_ALS)]
# sim-vs-obs slope agreement across plots
w <- dcast(B, id_plot+LAI_ALS ~ src, value.var="slope")
r_simobs <- cor(w[["Observed (HOBO)"]], w[["MuSICA (STATIC_ALS)"]], use="complete.obs")
cat(sprintf("r(slope_obs ~ LAI) = %.2f\n", r_obs))
cat(sprintf("r(slope_MuSICA ~ LAI) = %.2f\n", r_mus))
cat(sprintf("r(slope_MuSICA, slope_obs) [sim-vs-obs agreement] = %.2f\n", r_simobs))

B[, src:=factor(src, levels=c("Observed (HOBO)","MuSICA (STATIC_ALS)"))]
g <- ggplot(B, aes(LAI_ALS, slope, colour=src, fill=src)) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.8, alpha=0.15) +
  geom_point(size=2, alpha=0.85) +
  geom_hline(yintercept=1, linetype="dotted", colour="grey40") +
  scale_colour_manual(values=unname(PAL_SENSOR), name=NULL) +
  scale_fill_manual(values=unname(PAL_SENSOR), name=NULL) +
  labs(x=expression("LiDAR LAI (one-sided,"~m^2/m^2*")"),
       y="buffering slope (sub-canopy ~ macroclimate)") +
  annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=1.4, size=3.2,
           label=sprintf("r = %.2f (obs)\nr = %.2f (MuSICA)", r_obs, r_mus)) +
  theme_article(11) + theme(legend.position="top")
ggsave_article("/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3/Fig5_premise_slope_LAI", g, 6.9, 5)
cat("DONE -> Fig5_premise_slope_LAI\n")
