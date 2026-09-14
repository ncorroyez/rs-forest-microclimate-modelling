# ==============================================================================
# Scatters validation MuSICA vs HOBO, par plot (53), pour ΔTmax et slope micro/macro,
# chacun TOUS LES JOURS et 10% LES + CHAUDS, v3.2.0 (legacy) vs v3.2.3 iter (yoyo).
# Observé calculé depuis les HOBO ; simulé depuis les nc REF (z05).
#   Rscript c1_cmp_scatter_metrics.R
# Out: FigCmp_scatter_metrics_v320_v323iter.png (+ comparaison_versions/figures)
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(ggplot2)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("pipeline/00_config.R")
Z_FIX <- 1.0; SHIFT <- 2L; HOTQ <- 0.90; ds <- CFG$date_seq
dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, ds)); setnames(dm, "Tmax_macro","Tmax_mac")
hot  <- dm[Tmax_mac >= quantile(Tmax_mac, HOTQ, na.rm=TRUE), date]
era5 <- as.data.table(build_era5_hourly(CFG$forcing_file, ds))   # time, Tair_era5

# ---- OBSERVÉ par plot (HOBO) ----
hobo <- as.data.table(read.csv(CFG$hobo_temp_csv))[position_sensor=="a" & !(id_plot %in% CFG$ids_to_remove)]
hobo[, datetime := as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hobo <- hobo[as.Date(datetime) %in% ds]; hobo[, time := floor_date(datetime,"hour")]
hh <- hobo[, .(t_obs = mean(t_hobo, na.rm=TRUE)), by=.(id_plot, time)]
hh[, date := as.Date(time)]
obs_daily <- hh[, .(Tmax=max(t_obs, na.rm=TRUE)), by=.(id_plot,date)]
obs_daily <- merge(obs_daily, dm, by="date")
hh <- merge(hh, era5, by="time")
obs <- hh[, {
  dd <- obs_daily[id_plot==.BY$id_plot]
  .(Tmax_all = mean(dd$Tmax - dd$Tmax_mac, na.rm=TRUE),
    Tmax_hot = mean(dd[date %in% hot, Tmax - Tmax_mac], na.rm=TRUE),
    slope_all = as.numeric(coef(lm(t_obs ~ Tair_era5, .SD))[2]),
    slope_hot = if (nrow(.SD[date %in% hot])>10) as.numeric(coef(lm(t_obs ~ Tair_era5, .SD[date %in% hot]))[2]) else NA_real_)
}, by=id_plot]

# ---- SIMULÉ par plot, par version (nc REF) — métriques identiques ----
metrics_one <- function(path) {
  out <- c(Tmax_all=NA_real_,Tmax_hot=NA_real_,slope_all=NA_real_,slope_hot=NA_real_)
  if (!file.exists(path)) return(out)
  nc <- tryCatch(nc_open(path), error=function(e) NULL); if (is.null(nc)) return(out); on.exit(nc_close(nc))
  if (!all(c("Tair_z","relative_height","veget_height_top") %in% names(nc$var))) return(out)
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu), tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE)
  zl <- rh*vh; Z <- Z_FIX
  if (Z<=zl[1]){il<-1L;ih<-1L;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1L;w<-(Z-zl[il])/(zl[ih]-zl[il])}
  Tc <- ((1-w)*Tk[il,] + w*Tk[ih,]) - 273.15
  tv <- t0 + dhours(th) - lubridate::hours(SHIFT)
  dd <- data.table(date=as.Date(floor_date(tv,"hour")), Tc=Tc)[date %in% ds, .(Tmax=max(Tc)), by=date]
  m1 <- merge(dd, dm, by="date")
  out["Tmax_all"]<-mean(m1$Tmax-m1$Tmax_mac,na.rm=TRUE); out["Tmax_hot"]<-mean(m1[date%in%hot,Tmax-Tmax_mac],na.rm=TRUE)
  mic <- data.table(time=floor_date(t0+dhours(th),"hour"), Tmic=Tc); mm <- merge(mic[as.Date(time)%in%ds], era5, by="time")
  out["slope_all"]<-as.numeric(coef(lm(Tmic~Tair_era5,mm))[2])
  mh<-mm[as.Date(time)%in%hot]; if(nrow(mh)>10) out["slope_hot"]<-as.numeric(coef(lm(Tmic~Tair_era5,mh))[2])
  out
}
VERS <- c("v3.2.0"="out_files/musica_hobo_z05/1111", "v3.2.3 iter"="out_files/musica_hobo_z05_iter/1111")
sim <- rbindlist(lapply(names(VERS), function(v) {
  rbindlist(lapply(obs$id_plot, function(id) {
    m <- metrics_one(file.path(VERS[[v]], sprintf("musica_out_HOBO_%s.nc", id)))
    data.table(id_plot=id, version=v, Tmax_all=m[["Tmax_all"]], Tmax_hot=m[["Tmax_hot"]],
               slope_all=m[["slope_all"]], slope_hot=m[["slope_hot"]])
  }))
}))
D <- melt(sim, id.vars=c("id_plot","version"), variable.name="mp", value.name="sim")
O <- melt(obs, id.vars="id_plot", variable.name="mp", value.name="obs")
D <- merge(D, O, by=c("id_plot","mp")); D <- D[is.finite(sim)&is.finite(obs)]
MPLAB <- c(Tmax_all="ΔTmax · tous", Tmax_hot="ΔTmax · chauds10%", slope_all="slope · tous", slope_hot="slope · chauds10%")
D[, mp := factor(mp, levels=names(MPLAB), labels=MPLAB)]
D[, version := factor(version, levels=names(VERS))]
lab <- D[, .(R2=cor(sim,obs)^2, rmse=sqrt(mean((sim-obs)^2)), n=.N), by=.(mp,version)]
lab[, txt := sprintf("R² = %.2f\nRMSE = %.2f", R2, rmse)]
cat("=== R² / RMSE par métrique×période×version ===\n"); print(lab[,.(mp,version,R2=round(R2,2),RMSE=round(rmse,2),n)][order(mp,version)])
# range COMMUN x=y par ligne (mp), partagé entre les 2 versions
lim <- D[, .(lo=min(c(sim,obs)), hi=max(c(sim,obs))), by=mp]
lab <- merge(lab, lim, by="mp")
vers <- levels(D$version)
blank <- rbindlist(lapply(vers, function(v) {
  rbind(lim[, .(mp, version=factor(v,levels=vers), obs=lo, sim=lo)],
        lim[, .(mp, version=factor(v,levels=vers), obs=hi, sim=hi)]) }))
p <- ggplot(D, aes(obs, sim)) +
  geom_blank(data=blank, aes(obs, sim)) +                         # force x-range == y-range par panneau
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey55") +
  geom_smooth(method="lm", se=FALSE, colour="#D7791B", linewidth=0.6) +
  geom_point(alpha=0.6, size=1.3, colour="#2C3E50") +
  geom_text(data=lab, aes(x=lo, y=hi, label=txt), hjust=0, vjust=1, size=2.7, colour="grey20") +
  facet_wrap(~ mp + version, ncol=2, scales="free") +
  labs(x="Observé (HOBO)", y="Simulé (MuSICA, 1 m)",
       title="Validation MuSICA vs HOBO par plot — ΔTmax & slope micro/macro (tous jours / 10% chauds)",
       subtitle="53 loggers. v3.2.0 (legacy validé) vs v3.2.3 iter (yoyo, vraie BLH). 1:1 tireté, ajustement orange.") +
  theme_bw(base_size=11) +
  theme(aspect.ratio=1,                                          # panneaux carrés
        panel.grid.minor=element_blank(), strip.text=element_text(face="bold", size=9),
        plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))
out <- "out_files/Chapter3/figures/FigCmp_scatter_metrics_v320_v323iter.png"
ggsave(out, p, width=7, height=13, dpi=200, bg="white")
file.copy(out, "Chapitre1/comparaison_versions/figures/", overwrite=TRUE)
cat("DONE ->", out, "\n")
