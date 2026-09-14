# ==============================================================================
# CONSOLIDATED RECAP — mechanistic (MuSICA v3.2.0 & v3.2.3-iter) vs statistical
# (Gril slope-and-equilibrium) microclimate buffering, JUNE-SEPTEMBER (JJAS),
# with advanced statistics. n=53 plots, all 53/53.
#
# Targets (observed HOBO, per plot):
#   gril_slope   = fixed slope of lmer(T_micro ~ T_macro + (1|Month))   [Gril]
#   gril_equil   = median month intercept/(1-slope)                     [Gril]
#   dTmax_summer = mean(Tmax_micro - Tmax_macro)  over JJAS
#   dTmax_hot    = same on hot days (macro daily Tmax >= p90)
# Predictors:
#   MECHANISTIC  : MuSICA STATIC_ALS v3.2.0 and v3.2.3-iter (pure forward)
#   STATISTICAL  : Gril LOO regression (canopy ; canopy+topo)
# Advanced stats: bootstrap 95% CI on Pearson r and RMSE (B=2000, plot resample);
#   Spearman (rank); Williams test for dependent correlations (is v3.2.0 better
#   than v3.2.3? is mechanistic vs statistical?); variance partitioning of the
#   observed slope between canopy and topography.
#   Rscript gril_recap_advanced_JJAS.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(lme4); library(ggplot2); library(patchwork)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

YEAR <- 2021L; MM <- 6:9                  # June-September
MET  <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Deciduous_Only"
OUT_TAB <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
OUT_FIG <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
NCDIRS <- c("v3.2.0"=file.path(CFG_C3$out_dir,"nc","STATIC_ALS"),
            "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter","STATIC_ALS"))
set.seed(42)

# ---- macro hourly + daily Tmax + hot days ------------------------------------
ncf <- nc_open(CFG_C3$forcing_file)
mh  <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                  T_macro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
mh  <- mh[year(time)==YEAR & month(time)%in%MM][, `:=`(hr=floor_date(time,"hour"), d=as.Date(time))]
md  <- mh[, .(Tmac_max=max(T_macro)), by=d]; hotd <- md[Tmac_max>=quantile(Tmac_max,.90), d]

# ---- per-plot metric extractors ----------------------------------------------
gril <- function(dt){ dt<-dt[is.finite(T_micro)&is.finite(T_macro)]
  if(uniqueN(dt$Month)<2||nrow(dt)<200) return(c(slope=NA,equil=NA))
  m<-suppressMessages(suppressWarnings(lmer(T_micro~T_macro+(1|Month),data=dt)))
  sl<-as.numeric(fixef(m)["T_macro"]); ic<-coef(m)$Month[,"(Intercept)"]
  c(slope=sl, equil=median(ic/(1-sl))) }
dtmax <- function(daily){ daily<-merge(daily, md, by="d"); daily[, dT:=Tmax-Tmac_max]
  c(dTmax_summer=mean(daily$dT,na.rm=TRUE), dTmax_hot=mean(daily$dT[daily$d%in%hotd],na.rm=TRUE)) }

# ---- features ----------------------------------------------------------------
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
tf <- c(CHM_std="std_res_10_m.tif", slope_t="slope_res_10_m.tif", aspect_sin="aspect_sin_res_10_m.tif",
        aspect_cos="aspect_cos_res_10_m.tif", northness="northness_slope_res_10_m.tif", twi="twi_res_10_m.tif")
stk <- rast(file.path(MET, tf)); names(stk) <- names(tf)
ex  <- as.data.table(terra::extract(stk, vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs=crs(stk))))[, -1]
df  <- cbind(df, ex); for (v in names(tf)) if (anyNA(df[[v]])) df[[v]][is.na(df[[v]])] <- mean(df[[v]], na.rm=TRUE)

# ---- observed -----------------------------------------------------------------
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[, time:=as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot%in%CFG_C3$ids_to_remove) & year(time)==YEAR & month(time)%in%MM,
         .(id_plot=as.character(id_plot), T_micro=t_hobo, hr=floor_date(time,"hour"))]
hb <- merge(hb, mh[,.(hr,T_macro,Month=month(d))], by="hr")
hbd <- hb[, .(Tmax=max(T_micro)), by=.(id_plot, d=as.Date(hr))]

# ---- sim reader ---------------------------------------------------------------
sim_metrics <- function(dir, id){ f<-file.path(dir, sprintf("musica_out_HOBO_%s.nc", id))
  if(!file.exists(f)) return(c(slope=NA,equil=NA,dTmax_summer=NA,dTmax_hot=NA))
  r<-as.data.table(get_tair_at_z(nc_open(f),1.0)); r<-r[year(time)==YEAR & month(time)%in%MM][, hr:=floor_date(time,"hour")]
  r<-merge(r, mh[,.(hr,T_macro,Month=month(d))], by="hr")
  g<-gril(r[,.(T_micro=Tair_sim,T_macro,Month)])
  rd<-r[,.(Tmax=max(Tair_sim)),by=.(d=as.Date(hr))]; c(g, dtmax(rd)) }

# ---- assemble all per-plot metrics -------------------------------------------
M <- rbindlist(lapply(df$id_plot, function(id){
  go<-gril(hb[id_plot==id,.(T_micro,T_macro,Month)]); do<-dtmax(hbd[id_plot==id,.(d,Tmax)])
  s0<-sim_metrics(NCDIRS["v3.2.0"],id); s3<-sim_metrics(NCDIRS["v3.2.3 iter"],id)
  data.table(id_plot=id,
    slope_obs=go["slope"], equil_obs=go["equil"], dTmaxS_obs=do["dTmax_summer"], dTmaxH_obs=do["dTmax_hot"],
    slope_v320=s0["slope"], equil_v320=s0["equil"], dTmaxS_v320=s0["dTmax_summer"], dTmaxH_v320=s0["dTmax_hot"],
    slope_v323=s3["slope"], equil_v323=s3["equil"], dTmaxS_v323=s3["dTmax_summer"], dTmaxH_v323=s3["dTmax_hot"]) }))
D <- merge(M, df[,.(id_plot,LAI_ALS,Hmax,fCover,VCI,LCV,CHM_std,northness,twi,slope_t)], by="id_plot")
fwrite(D, file.path(OUT_TAB,"recap_JJAS_plotdata.csv"))
cat(sprintf("n=%d | finite per metric:\n", nrow(D)))
for(c in grep("_obs$|_v320$|_v323$",names(D),value=TRUE)) cat(sprintf("  %-12s %d/53\n",c,sum(is.finite(D[[c]]))))

# ---- advanced-stats helpers --------------------------------------------------
r2<-function(o,p){k<-is.finite(o)&is.finite(p);1-sum((o[k]-p[k])^2)/sum((o[k]-mean(o[k]))^2)}
rmse<-function(o,p){k<-is.finite(o)&is.finite(p);sqrt(mean((o[k]-p[k])^2))}
mae<-function(o,p){k<-is.finite(o)&is.finite(p);mean(abs(o[k]-p[k]))}
boot_ci<-function(o,p,fn,B=2000){ k<-which(is.finite(o)&is.finite(p)); v<-replicate(B,{i<-sample(k,length(k),TRUE);fn(o[i],p[i])})
  quantile(v,c(.025,.975),na.rm=TRUE) }
do_loo<-function(form,data){ vars<-all.vars(form); dd<-data[complete.cases(data[,..vars])]
  p<-rep(NA_real_,nrow(dd)); for(i in seq_len(nrow(dd))){ m<-lm(form,dd[-i]); p[i]<-predict(m,dd[i]) }
  list(obs=dd[[vars[1]]],pred=p) }
# Williams test: compare r(x,y) vs r(x,z) with shared x (dependent correlations)
williams<-function(r_xy,r_xz,r_yz,n){ Rdet<-1-r_xy^2-r_xz^2-r_yz^2+2*r_xy*r_xz*r_yz
  t<-(r_xy-r_xz)*sqrt((n-1)*(1+r_yz)/(2*((n-1)/(n-3))*Rdet+((r_xy+r_xz)^2/4)*(1-r_yz)^3))
  c(t=t, df=n-3, p=2*pt(-abs(t),n-3)) }

# ---- head-to-head, advanced, per target --------------------------------------
TARG <- list(gril_slope=c("slope_obs","slope_v320","slope_v323"),
             dTmax_summer=c("dTmaxS_obs","dTmaxS_v320","dTmaxS_v323"),
             dTmax_hot=c("dTmaxH_obs","dTmaxH_v320","dTmaxH_v323"))
STAT_FORM <- list(gril_slope=slope_obs~LAI_ALS+fCover+Hmax,
                  dTmax_summer=dTmaxS_obs~LAI_ALS+fCover+Hmax,
                  dTmax_hot=dTmaxH_obs~LAI_ALS+fCover+Hmax)
STAT_FORMT<- list(gril_slope=slope_obs~LAI_ALS+fCover+Hmax+northness+twi,
                  dTmax_summer=dTmaxS_obs~LAI_ALS+fCover+Hmax+northness+twi,
                  dTmax_hot=dTmaxH_obs~LAI_ALS+fCover+Hmax+northness+twi)
RES<-list(); WIL<-list()
for(tg in names(TARG)){ o<-D[[TARG[[tg]][1]]]; v0<-D[[TARG[[tg]][2]]]; v3<-D[[TARG[[tg]][3]]]
  ec<-do_loo(STAT_FORM[[tg]],D); et<-do_loo(STAT_FORMT[[tg]],D)
  rows<-list(list("MuSICA v3.2.0",o,v0), list("MuSICA v3.2.3",o,v3),
             list("Gril LOO canopy",ec$obs,ec$pred), list("Gril LOO canopy+topo",et$obs,et$pred))
  RES[[tg]]<-rbindlist(lapply(rows,function(z){ ci_r<-boot_ci(z[[2]],z[[3]],function(a,b)cor(a,b)); ci_e<-boot_ci(z[[2]],z[[3]],rmse)
    data.table(target=tg, model=z[[1]], r=cor(z[[2]],z[[3]]), r_lo=ci_r[1], r_hi=ci_r[2],
               rho=cor(z[[2]],z[[3]],method="spearman"), R2=r2(z[[2]],z[[3]]),
               RMSE=rmse(z[[2]],z[[3]]), RMSE_lo=ci_e[1], RMSE_hi=ci_e[2], MAE=mae(z[[2]],z[[3]])) }))
  # Williams: v3.2.0 vs v3.2.3 (shared obs); and v3.2.0 vs statistical(+topo)
  nfin<-sum(is.finite(o)&is.finite(v0)&is.finite(v3))
  w1<-williams(cor(o,v0),cor(o,v3),cor(v0,v3),nfin)
  w2<-williams(cor(o,v0),cor(et$obs,et$pred),cor(v0,et$pred),length(et$pred))
  WIL[[tg]]<-rbind(
    data.table(target=tg, contrast="MuSICA v3.2.0 vs v3.2.3",     t=w1["t"], df=w1["df"], p=w1["p"]),
    data.table(target=tg, contrast="MuSICA v3.2.0 vs Gril+topo",  t=w2["t"], df=w2["df"], p=w2["p"])) }
RES<-rbindlist(RES); WIL<-rbindlist(WIL)
cat("\n=== HEAD-TO-HEAD (JJAS) — Pearson r [95% boot CI], Spearman rho, R2, RMSE [CI], MAE ===\n")
print(RES[,.(target,model,r=round(r,3),CI=sprintf("[%.2f,%.2f]",r_lo,r_hi),rho=round(rho,3),R2=round(R2,3),RMSE=round(RMSE,3),MAE=round(MAE,3))])
cat("\n=== WILLIAMS test for dependent correlations (shared = observed) ===\n")
print(WIL[,.(target,contrast,t=round(t,2),df=round(df),p=round(p,4))])
fwrite(RES, file.path(OUT_TAB,"recap_JJAS_headtohead.csv")); fwrite(WIL, file.path(OUT_TAB,"recap_JJAS_williams.csv"))

# ---- equilibrium (rank + robust), all 53 -------------------------------------
medae<-function(o,p){k<-is.finite(o)&is.finite(p);median(abs(o[k]-p[k]))}
cat("\n=== EQUILIBRIUM (all 53; Spearman + medAE) ===\n")
for(v in c("equil_v320","equil_v323")) cat(sprintf("  %-11s Spearman=%+.3f medAE=%.2f degC (Pearson=%+.3f)\n",
  sub("equil_","",v), cor(D$equil_obs,D[[v]],method="spearman"), medae(D$equil_obs,D[[v]]), cor(D$equil_obs,D[[v]])))
eqs<-do_loo(equil_obs~LAI_ALS+fCover+Hmax+northness+twi,D)
cat(sprintf("  Gril LOO    Spearman=%+.3f medAE=%.2f degC (Pearson=%+.3f)\n",
  cor(eqs$obs,eqs$pred,method="spearman"), medae(eqs$obs,eqs$pred), cor(eqs$obs,eqs$pred)))

# ---- attribution + variance partitioning (observed slope) --------------------
cat("\n=== ATTRIBUTION: observed Gril slope ~ canopy + topo (standardized) ===\n")
Z <- copy(D); for(v in c("LAI_ALS","fCover","Hmax","CHM_std","northness","twi","slope_t")) Z[[v]]<-as.numeric(scale(Z[[v]]))
ma <- lm(slope_obs ~ LAI_ALS+fCover+northness+twi, Z); ci<-confint(ma)
ATT <- data.table(term=names(coef(ma)), beta=coef(ma), lo=ci[,1], hi=ci[,2], p=summary(ma)$coefficients[,4])
print(ATT[,.(term,beta=round(beta,4),CI=sprintf("[%.3f,%.3f]",lo,hi),p=round(p,4))])
R2c <- summary(lm(slope_obs~LAI_ALS+fCover+Hmax,D))$r.squared
R2t <- summary(lm(slope_obs~northness+twi+slope_t,D))$r.squared
R2f <- summary(lm(slope_obs~LAI_ALS+fCover+Hmax+northness+twi+slope_t,D))$r.squared
cat(sprintf("\n  Variance partition (observed slope):  canopy-only R2=%.3f | topo-only R2=%.3f | full R2=%.3f\n", R2c,R2t,R2f))
cat(sprintf("  unique canopy=%.3f | unique topo=%.3f | shared=%.3f\n", R2f-R2t, R2f-R2c, (R2c+R2t)-R2f))
cat("\n  mechanistic residual ~ topo (per binary):\n")
for(v in c("v320","v323")){ D[, rr:=slope_obs-get(paste0("slope_",v))]; mr<-lm(rr~northness+twi+slope_t,D)
  cat(sprintf("    [%s] R2=%.3f p=%.4f | northness p=%.4f r=%+.3f\n", v, summary(mr)$r.squared,
    pf(summary(mr)$fstatistic[1],summary(mr)$fstatistic[2],summary(mr)$fstatistic[3],lower.tail=FALSE),
    summary(mr)$coefficients["northness",4], cor(D$rr,D$northness))) }; D[,rr:=NULL]

# ---- RECAP FIGURE ------------------------------------------------------------
RES[, model:=factor(model, levels=c("MuSICA v3.2.0","MuSICA v3.2.3","Gril LOO canopy","Gril LOO canopy+topo"))]
RES[, kind:=fifelse(grepl("MuSICA",model),"Mechanistic","Statistical")]
pA <- ggplot(RES[target!="dTmax_hot"], aes(r, model, colour=kind)) +
  geom_vline(xintercept=0, linetype="dotted", colour="grey70") +
  geom_vline(xintercept=1, linetype="dashed", colour="grey70") +
  geom_errorbarh(aes(xmin=r_lo, xmax=r_hi), height=.25, linewidth=.6) + geom_point(size=2.4) +
  facet_wrap(~target, ncol=1) + scale_colour_manual(values=c(Mechanistic="#0072B2",Statistical="#D55E00"), name=NULL) +
  labs(x="Pearson r vs observed [95% bootstrap CI]", y=NULL, subtitle="(a) Between-plot ranking skill (JJAS, n=53)") +
  theme_article(10) + theme(legend.position="top")
ec<-do_loo(STAT_FORMT[["gril_slope"]],D)
F11 <- rbind(data.table(model="MuSICA v3.2.0", obs=D$slope_obs, pred=D$slope_v320),
             data.table(model="MuSICA v3.2.3", obs=D$slope_obs, pred=D$slope_v323),
             data.table(model="Gril LOO +topo", obs=ec$obs, pred=ec$pred))
F11[, model:=factor(model, levels=c("MuSICA v3.2.0","MuSICA v3.2.3","Gril LOO +topo"))]
lim<-range(c(F11$obs,F11$pred),na.rm=TRUE)
pB <- ggplot(F11, aes(pred,obs)) + geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55") +
  geom_vline(xintercept=1,colour="grey85",linewidth=.3)+geom_hline(yintercept=1,colour="grey85",linewidth=.3)+
  geom_point(size=1.3,colour="#0072B2",alpha=.8) + facet_wrap(~model,nrow=1) + coord_fixed(xlim=lim,ylim=lim) +
  labs(x="Predicted Gril slope", y="Observed", subtitle="(b) Gril slope: sims compress toward 1") + theme_article(9)
pC <- ggplot(ATT[term!="(Intercept)"], aes(beta, reorder(term,beta))) +
  geom_vline(xintercept=0,linetype="dotted",colour="grey60") +
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=.2,colour="#1A9850")+geom_point(size=2.4,colour="#1A9850")+
  labs(x="Standardized effect on observed slope [95% CI]", y=NULL, subtitle="(c) Attribution: canopy drives buffering, topo ~0") + theme_article(10)
g <- (pA | (pB / pC)) + plot_layout(widths=c(1,1.15))
ggsave_article(file.path(OUT_FIG,"Fig_recap_JJAS_advanced"), g, 13, 7)
cat("\nDONE -> Fig_recap_JJAS_advanced + recap_JJAS_{plotdata,headtohead,williams}.csv\n")
