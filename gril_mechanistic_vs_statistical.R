# ==============================================================================
# Mechanistic (MuSICA v3.2.0 AND v3.2.3-iter) vs statistical (Gril et al.
# slope-and-equilibrium) comparison of forest microclimate buffering.
#
# Gril method (slope_equilibrium_MORFO_2022.Rmd), leaf-on, per plot:
#   lmer(T_micro ~ T_macro + (1|Month))
#   slope       = fixed T_macro coefficient (<1 buffering, >1 amplification)
#   equilibrium = intercept_month / (1 - slope)   [per month; T where micro==macro]
#
# Head-to-head on ONE target = OBSERVED Gril slope/equilibrium (per plot, n=53):
#   - MECHANISTIC: MuSICA STATIC_ALS, pure forward (never saw the loggers), BOTH
#                  v3.2.0 (article binary) and v3.2.3 + abl_flag=iter (ABL-coupled)
#   - STATISTICAL: Gril regression of obs metric ~ canopy(+topo), scored by LOO
# Caveats handled: equilibrium = i/(1-slope) blows up as slope->1 (where the sims
# pin, v3.2.3 worst); lead with slope, equilibrium robust (median, |1-slope|>=thr).
# T_sim~T_forcing is mechanically tighter than T_obs~T_forcing (sim derived from
# the forcing through the column) -> WHY sim slopes pin near 1, reported as finding.
#   Rscript gril_mechanistic_vs_statistical.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(lme4); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

YEAR <- 2021L; LEAFON <- 5:10
MET  <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Deciduous_Only"
OUT_TAB <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
OUT_FIG <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
EQ_THR <- 0.05
NCDIRS <- c("v3.2.0"=file.path(CFG_C3$out_dir,"nc","STATIC_ALS"),
            "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter","STATIC_ALS"))

gril_fit <- function(dt){
  dt <- dt[is.finite(T_micro) & is.finite(T_macro)]
  if (uniqueN(dt$Month) < 2 || nrow(dt) < 200) return(list(slope=NA, equil=NA))
  m <- suppressMessages(suppressWarnings(lmer(T_micro ~ T_macro + (1|Month), data=dt)))
  sl <- as.numeric(fixef(m)["T_macro"]); ic <- coef(m)$Month[,"(Intercept)"]
  # equilibrium computed for EVERY plot (Gril keeps all; extremes handled in the
  # comparison via rank stats), median over the monthly intercepts.
  list(slope=sl, equil=median(ic/(1-sl)))
}
fit_nc <- function(dir, id, mh){ f <- file.path(dir, sprintf("musica_out_HOBO_%s.nc", id))
  if (!file.exists(f)) return(list(slope=NA, equil=NA))
  r <- as.data.table(get_tair_at_z(nc_open(f),1.0))
  r <- r[year(time)==YEAR & month(time)%in%LEAFON][, hr:=floor_date(time,"hour")]
  r <- merge(r, mh[,.(hr,T_macro,Month=month(time))], by="hr")
  gril_fit(r[, .(T_micro=Tair_sim, T_macro, Month)]) }

# ---- macro hourly forcing (same Tair for both binaries) ----------------------
ncf <- nc_open(CFG_C3$forcing_file)
mh  <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                  T_macro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
mh  <- mh[year(time)==YEAR & month(time)%in%LEAFON][, hr:=floor_date(time,"hour")]

# ---- topo + structure features -----------------------------------------------
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
tf <- c(CHM_std="std_res_10_m.tif", slope_t="slope_res_10_m.tif", aspect_sin="aspect_sin_res_10_m.tif",
        aspect_cos="aspect_cos_res_10_m.tif", northness="northness_slope_res_10_m.tif", twi="twi_res_10_m.tif")
stk <- rast(file.path(MET, tf)); names(stk) <- names(tf)
ex  <- as.data.table(terra::extract(stk, vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs=crs(stk))))[, -1]
df  <- cbind(df, ex); for (v in names(tf)) if (anyNA(df[[v]])) df[[v]][is.na(df[[v]])] <- mean(df[[v]], na.rm=TRUE)

# ---- observed HOBO hourly ----------------------------------------------------
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[, time:=as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot%in%CFG_C3$ids_to_remove) & year(time)==YEAR & month(time)%in%LEAFON,
         .(id_plot=as.character(id_plot), T_micro=t_hobo, hr=floor_date(time,"hour"))]
hb <- merge(hb, mh[,.(hr,T_macro,Month=month(time))], by="hr")

# ---- per-plot Gril fits: obs + both binaries ---------------------------------
G <- rbindlist(lapply(df$id_plot, function(id){
  o  <- gril_fit(hb[id_plot==id, .(T_micro, T_macro, Month)])
  s0 <- fit_nc(NCDIRS["v3.2.0"], id, mh); s3 <- fit_nc(NCDIRS["v3.2.3 iter"], id, mh)
  data.table(id_plot=id, slope_obs=o$slope, equil_obs=o$equil,
             slope_v320=s0$slope, equil_v320=s0$equil, slope_v323=s3$slope, equil_v323=s3$equil) }))
D <- merge(G, df[,.(id_plot, LAI_ALS, Hmax, fCover, VCI, LCV, CHM_std, northness, twi, slope_t)], by="id_plot")
fwrite(D, file.path(OUT_TAB,"gril_slope_equil_plotdata.csv"))

cat("\n=== pole diagnostic: summary(1 - slope) ===\n")
for (nm in c("slope_obs","slope_v320","slope_v323")){ cat(sprintf("  %-11s", sub("slope_","",nm)))
  print(round(summary(1-D[[nm]]),3)) }
cat(sprintf("  plots in pole |1-slope|<%.2f : obs=%d  v3.2.0=%d  v3.2.3=%d\n", EQ_THR,
  sum(abs(1-D$slope_obs)<EQ_THR), sum(abs(1-D$slope_v320)<EQ_THR), sum(abs(1-D$slope_v323)<EQ_THR)))

# ---- metrics helpers ---------------------------------------------------------
r2<-function(o,p){k<-is.finite(o)&is.finite(p);1-sum((o[k]-p[k])^2)/sum((o[k]-mean(o[k]))^2)}
rmse<-function(o,p){k<-is.finite(o)&is.finite(p);sqrt(mean((o[k]-p[k])^2))}
do_loo <- function(form, data){ vars<-all.vars(form); dd<-data[complete.cases(data[,..vars])]
  p<-rep(NA_real_,nrow(dd)); for(i in seq_len(nrow(dd))){ m<-lm(form, dd[-i]); p[i]<-predict(m, dd[i]) }
  list(obs=dd[[vars[1]]], pred=p) }

# ---- head-to-head: SLOPE -----------------------------------------------------
cat("\n=== HEAD-TO-HEAD on OBSERVED Gril slope (target), n=53 ===\n")
ec <- do_loo(slope_obs ~ LAI_ALS + fCover + Hmax, D)
et <- do_loo(slope_obs ~ LAI_ALS + fCover + Hmax + northness + twi, D)
SL <- rbind(
  data.table(model="MuSICA v3.2.0 (forward)",  cor=cor(D$slope_obs,D$slope_v320,use="complete.obs"), R2=r2(D$slope_obs,D$slope_v320), RMSE=rmse(D$slope_obs,D$slope_v320)),
  data.table(model="MuSICA v3.2.3 (forward)",  cor=cor(D$slope_obs,D$slope_v323,use="complete.obs"), R2=r2(D$slope_obs,D$slope_v323), RMSE=rmse(D$slope_obs,D$slope_v323)),
  data.table(model="Gril LOO: canopy",         cor=cor(ec$obs,ec$pred), R2=r2(ec$obs,ec$pred), RMSE=rmse(ec$obs,ec$pred)),
  data.table(model="Gril LOO: canopy+topo",    cor=cor(et$obs,et$pred), R2=r2(et$obs,et$pred), RMSE=rmse(et$obs,et$pred)))
print(SL[, .(model, cor=round(cor,3), R2=round(R2,3), RMSE=round(RMSE,4))])

# ---- head-to-head: EQUILIBRIUM (all 53; rank + robust error) -----------------
# Equilibrium = i/(1-slope) is heavy-tailed near slope=1 (the sims), so the
# head-to-head is read on Spearman rank cor and median |error| (Pearson/RMSE are
# pole-dominated and shown only for reference).
cat("\n=== HEAD-TO-HEAD on OBSERVED equilibrium (all 53; Spearman + median|err|) ===\n")
medae <- function(o,p){k<-is.finite(o)&is.finite(p);median(abs(o[k]-p[k]))}
for (v in c("v320","v323")){ sc<-paste0("equil_",v); DE<-D[is.finite(equil_obs)&is.finite(get(sc))]
  cat(sprintf("  MuSICA %-6s n=%d/53  Spearman=%+.3f  medAE=%.2f degC  (Pearson=%+.3f RMSE=%.1f)\n", v, nrow(DE),
    cor(DE$equil_obs,DE[[sc]],method="spearman"), medae(DE$equil_obs,DE[[sc]]),
    cor(DE$equil_obs,DE[[sc]]), rmse(DE$equil_obs,DE[[sc]]))) }
DEs <- D[is.finite(equil_obs)]; eqs <- do_loo(equil_obs ~ LAI_ALS + fCover + Hmax + northness + twi, DEs)
cat(sprintf("  Gril LOO    n=%d/53  Spearman=%+.3f  medAE=%.2f degC  (Pearson=%+.3f RMSE=%.1f)\n",
  length(eqs$obs), cor(eqs$obs,eqs$pred,method="spearman"), medae(eqs$obs,eqs$pred),
  cor(eqs$obs,eqs$pred), rmse(eqs$obs,eqs$pred)))

# ---- attribution + mechanistic residual (both binaries) ----------------------
cat("\n=== STATISTICAL attribution: lm(slope_obs ~ LAI + fCover + northness + twi) ===\n")
ma <- lm(slope_obs ~ LAI_ALS + fCover + northness + twi, D); print(round(summary(ma)$coefficients,4))
cat(sprintf("  adj R2 = %.3f\n", summary(ma)$adj.r.squared))
cat("\n=== MECHANISTIC residual ~ topo, per binary ===\n")
for (v in c("v320","v323")){ D[, rr := slope_obs - get(paste0("slope_",v))]
  mr<-lm(rr ~ northness + twi + slope_t, D)
  cat(sprintf("  [%s] model R2=%.3f p=%.4f | northness: coef=%+.4f p=%.4f ; r(res,northness)=%+.3f p=%.4f\n",
    v, summary(mr)$r.squared, pf(summary(mr)$fstatistic[1],summary(mr)$fstatistic[2],summary(mr)$fstatistic[3],lower.tail=FALSE),
    coef(mr)["northness"], summary(mr)$coefficients["northness",4],
    cor.test(D$rr,D$northness)$estimate, cor.test(D$rr,D$northness)$p.value)) }
D[, rr:=NULL]

# ---- Figure: slope head-to-head 1:1, 4 panels --------------------------------
F <- rbind(
  data.table(model="MuSICA v3.2.0 (forward)", obs=D$slope_obs, pred=D$slope_v320),
  data.table(model="MuSICA v3.2.3 (forward)", obs=D$slope_obs, pred=D$slope_v323),
  data.table(model="Gril LOO (canopy)",       obs=ec$obs, pred=ec$pred),
  data.table(model="Gril LOO (canopy+topo)",  obs=et$obs, pred=et$pred))
F[, model:=factor(model, levels=c("MuSICA v3.2.0 (forward)","MuSICA v3.2.3 (forward)","Gril LOO (canopy)","Gril LOO (canopy+topo)"))]
ann <- F[, .(r=cor(obs,pred,use="complete.obs"), rmse=rmse(obs,pred)), by=model]
lim <- range(c(F$obs,F$pred), na.rm=TRUE)
g <- ggplot(F, aes(pred,obs)) + geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55") +
  geom_vline(xintercept=1,colour="grey80",linewidth=.3) + geom_hline(yintercept=1,colour="grey80",linewidth=.3) +
  geom_point(size=1.6, colour="#0072B2", alpha=.85) +
  geom_text(data=ann, aes(x=-Inf,y=Inf,label=sprintf("r=%.2f  RMSE=%.3f",r,rmse)), hjust=-0.06, vjust=1.5, size=2.9, inherit.aes=FALSE) +
  facet_wrap(~model, nrow=1) + coord_fixed(xlim=lim, ylim=lim) +
  labs(x="Predicted Gril slope", y="Observed Gril slope (HOBO)",
       subtitle=paste0("Mechanistic (MuSICA v3.2.0 & v3.2.3, pure forward) vs statistical (Gril LOO) on the Gril slope, leaf-on 2021. ",
                       "n=53. Dashed=1:1, grey=slope 1.\nSim slopes compress toward 1 (1 m air slaved to the forcing); v3.2.3 compresses hardest.")) +
  theme_article(10)
ggsave_article(file.path(OUT_FIG,"Fig_gril_mechanistic_vs_statistical_slope"), g, 12, 4.0)
cat("\nDONE -> Fig_gril_mechanistic_vs_statistical_slope + gril_slope_equil_plotdata.csv\n")
