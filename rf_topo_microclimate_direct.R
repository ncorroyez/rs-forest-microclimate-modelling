# ==============================================================================
# Direct empirical microclimate model + the powered topo test.
#
# Per plot (n=53 HOBO), summer (JJAS, 1 m, no shift — same extraction as the
# time-series figures):
#   obs_beta   = slope lm(Tobs ~ Tmacro) over all hours      (beta<1 = buffering)
#   obs_dTmaxS = mean(Tmax_obs - Tmax_macro) over summer
#   obs_dTmaxH = same over hot days (macro daily Tmax >= p90)
#   mus_*      = identical metrics from MuSICA STATIC_ALS v3.2.0 (physical baseline)
#
# Q1 (RF, gradient reproduction): does a parsimonious empirical model of the
#     observed gradient beat MuSICA's collapsed one? cor(obs, .) for MuSICA vs a
#     LOO lm(obs ~ LAI + fCover + Hmax). RF kept light; n=53 caveat acknowledged.
# Q2 (THE topo test, powered & interpretable): residual = obs - MuSICA, then
#     lm(residual ~ TWI + northness + slope) — 2-3 terms by physical rationale.
#     Partial r / signs say whether terrain captures what the 1-D column cannot.
#   Rscript rf_topo_microclimate_direct.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

WIN  <- as.Date(c("2021-06-01","2021-09-30"))   # JJAS plateau
MET  <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Deciduous_Only"
OUT_TAB <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
OUT_FIG <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
beta <- function(y,x){ ok<-is.finite(x)&is.finite(y); if(sum(ok)<24) return(NA_real_)
  as.numeric(coef(lm(y[ok]~x[ok]))[2]) }

# ---- macro hourly + daily Tmax + hot-day flag --------------------------------
ncf <- nc_open(CFG_C3$forcing_file)
mh  <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                  Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
mh  <- mh[as.Date(time)>=WIN[1] & as.Date(time)<=WIN[2]][, `:=`(hr=floor_date(time,"hour"), d=as.Date(time))]
md  <- mh[, .(Tmac_max=max(Tmacro,na.rm=TRUE)), by=d]
hotd <- md[Tmac_max >= quantile(Tmac_max, .90, na.rm=TRUE), d]
cat(sprintf("summer days=%d  hot days(>=p90)=%d\n", nrow(md), length(hotd)))

# ---- df + topo/structure features --------------------------------------------
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
tf <- c(CHM_std="std_res_10_m.tif", slope="slope_res_10_m.tif", aspect_sin="aspect_sin_res_10_m.tif",
        aspect_cos="aspect_cos_res_10_m.tif", northness="northness_slope_res_10_m.tif", twi="twi_res_10_m.tif")
stk <- rast(file.path(MET, tf)); names(stk) <- names(tf)
ex  <- as.data.table(terra::extract(stk, vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs=crs(stk))))[, -1]
df  <- cbind(df, ex)
for (v in names(tf)) if (anyNA(df[[v]])) df[[v]][is.na(df[[v]])] <- mean(df[[v]], na.rm=TRUE)

# ---- observed metrics per plot -----------------------------------------------
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[, time:=as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) &
         as.Date(time)>=WIN[1] & as.Date(time)<=WIN[2],
         .(id_plot=as.character(id_plot), Tobs=t_hobo, hr=floor_date(time,"hour"), d=as.Date(time))]
hb <- merge(hb, mh[,.(hr,Tmacro)], by="hr")
obs_b  <- hb[, .(obs_beta=beta(Tobs,Tmacro)), by=id_plot]
hbd <- hb[, .(Tmax_obs=max(Tobs,na.rm=TRUE)), by=.(id_plot,d)]
hbd <- merge(hbd, md, by="d")[, dT:=Tmax_obs-Tmac_max]
obs_d <- hbd[, .(obs_dTmaxS=mean(dT,na.rm=TRUE),
                 obs_dTmaxH=mean(dT[d %in% hotd],na.rm=TRUE)), by=id_plot]

# ---- MuSICA STATIC_ALS v3.2.0 metrics per plot -------------------------------
read_micro <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc)
  if(is.null(r)||!nrow(r))return(NULL); as.data.table(r)[as.Date(time)>=WIN[1] & as.Date(time)<=WIN[2]] }
ncdir <- file.path(CFG_C3$out_dir,"nc","STATIC_ALS")
mus <- rbindlist(lapply(df$id_plot, function(id){
  f<-file.path(ncdir, sprintf("musica_out_HOBO_%s.nc", id)); if(!file.exists(f))return(NULL)
  r<-read_micro(f); if(is.null(r))return(NULL); r[, hr:=floor_date(time,"hour")]; r<-merge(r, mh[,.(hr,Tmacro)], by="hr")
  rd<-r[, .(Tmax_sim=max(Tair_sim,na.rm=TRUE)), by=.(d=as.Date(time))]; rd<-merge(rd, md, by="d")[, dT:=Tmax_sim-Tmac_max]
  data.table(id_plot=id, mus_beta=beta(r$Tair_sim, r$Tmacro),
             mus_dTmaxS=mean(rd$dT,na.rm=TRUE), mus_dTmaxH=mean(rd$dT[rd$d %in% hotd],na.rm=TRUE)) }))

# ---- assemble ----------------------------------------------------------------
D <- Reduce(function(a,b) merge(a,b,by="id_plot"), list(obs_b, obs_d, mus,
            df[,.(id_plot, LAI_ALS, Hmax, fCover, VCI, LCV, CHM_std, slope, aspect_sin, aspect_cos, northness, twi)]))
cat(sprintf("assembled n=%d plots\n", nrow(D)))
fwrite(D, file.path(OUT_TAB,"microclim_direct_plotdata.csv"))

r2  <- function(o,p){ ok<-is.finite(o)&is.finite(p); 1 - sum((o[ok]-p[ok])^2)/sum((o[ok]-mean(o[ok]))^2) }
rmse<- function(o,p){ ok<-is.finite(o)&is.finite(p); sqrt(mean((o[ok]-p[ok])^2)) }
loo_lm <- function(form, data){ p<-rep(NA_real_,nrow(data)); for(i in seq_len(nrow(data))){
  m<-lm(form, data[-i]); p[i]<-predict(m, data[i]) }; p }

# ---- Q1: empirical vs MuSICA (gradient reproduction) -------------------------
TARG <- list(beta=c(obs="obs_beta",mus="mus_beta"),
             dTmaxS=c(obs="obs_dTmaxS",mus="mus_dTmaxS"),
             dTmaxH=c(obs="obs_dTmaxH",mus="mus_dTmaxH"))
q1 <- rbindlist(lapply(names(TARG), function(tg){ o<-D[[TARG[[tg]]["obs"]]]; mu<-D[[TARG[[tg]]["mus"]]]
  fe <- as.formula(paste(TARG[[tg]]["obs"], "~ LAI_ALS + fCover + Hmax"))
  pe <- loo_lm(fe, D)
  data.table(target=tg,
    MuSICA_cor=cor(o,mu,use="complete.obs"), MuSICA_R2=r2(o,mu), MuSICA_RMSE=rmse(o,mu),
    Emp_cor=cor(o,pe,use="complete.obs"),    Emp_R2_loo=r2(o,pe), Emp_RMSE_loo=rmse(o,pe)) }))
cat("\n=== Q1: between-plot gradient — MuSICA(physical) vs empirical LOO lm(obs~LAI+fCover+Hmax) ===\n")
print(q1[, lapply(.SD, function(x) if(is.numeric(x)) round(x,3) else x)])

# ---- Q2: residual ~ topo (THE powered test) ----------------------------------
D[, `:=`(res_beta = obs_beta-mus_beta, res_dTmaxS = obs_dTmaxS-mus_dTmaxS, res_dTmaxH = obs_dTmaxH-mus_dTmaxH)]
topo3 <- c("twi","northness","slope")
q2_fit <- list(); q2_pair <- list()
for (rv in c("res_beta","res_dTmaxS","res_dTmaxH")){
  f <- as.formula(paste(rv, "~ twi + northness + slope")); m <- lm(f, D); s <- summary(m)
  q2_fit[[rv]] <- data.table(residual=rv, term=rownames(s$coefficients),
    coef=s$coefficients[,1], t=s$coefficients[,3], p=s$coefficients[,4],
    model_R2=s$r.squared, model_p=pf(s$fstatistic[1],s$fstatistic[2],s$fstatistic[3],lower.tail=FALSE))
  q2_pair[[rv]] <- rbindlist(lapply(topo3, function(tv){ ct<-cor.test(D[[rv]], D[[tv]])
    data.table(residual=rv, topo=tv, pearson_r=ct$estimate, p=ct$p.value) }))
}
Q2FIT <- rbindlist(q2_fit); Q2PAIR <- rbindlist(q2_pair)
cat("\n=== Q2: lm(residual ~ TWI + northness + slope)  [residual = obs - MuSICA] ===\n")
print(Q2FIT[, .(residual, term, coef=round(coef,4), t=round(t,2), p=round(p,4),
                model_R2=round(model_R2,3), model_p=round(model_p,4))])
cat("\n--- pairwise Pearson r(residual, topo) ---\n")
print(Q2PAIR[, .(residual, topo, pearson_r=round(pearson_r,3), p=round(p,4))])
fwrite(Q2FIT, file.path(OUT_TAB,"microclim_direct_Q2_topo_lm.csv"))
fwrite(q1,    file.path(OUT_TAB,"microclim_direct_Q1_gradient.csv"))

# ---- Figure 1: Q1 obs vs predicted 1:1 (beta + dTmaxS) ------------------------
mk11 <- function(tg, lab){ o<-D[[TARG[[tg]]["obs"]]]; mu<-D[[TARG[[tg]]["mus"]]]
  fe<-as.formula(paste(TARG[[tg]]["obs"],"~ LAI_ALS + fCover + Hmax")); pe<-loo_lm(fe,D)
  rbind(data.table(target=lab, model="MuSICA (physical)", obs=o, pred=mu),
        data.table(target=lab, model="Empirical LOO (LAI+fCover+Hmax)", obs=o, pred=pe)) }
F1 <- rbind(mk11("beta","beta (micro/macro)"), mk11("dTmaxS","dTmax summer (degC)"))
ann <- F1[, .(r=cor(obs,pred,use="complete.obs"), rmse=rmse(obs,pred)), by=.(target,model)]
g1 <- ggplot(F1, aes(pred,obs)) + geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey55") +
  geom_point(size=1.6, colour="#0072B2", alpha=.8) +
  geom_text(data=ann, aes(x=-Inf,y=Inf,label=sprintf("r=%.2f  RMSE=%.2f",r,rmse)),
            hjust=-0.08, vjust=1.4, size=3, inherit.aes=FALSE) +
  facet_grid(target ~ model, scales="free") +
  labs(x="Predicted", y="Observed (HOBO)",
       subtitle="Q1 — between-plot gradient: physical (MuSICA STATIC_ALS v3.2.0) vs empirical LOO. n=53. Dashed = 1:1.") +
  theme_article(10)
ggsave_article(file.path(OUT_FIG,"Fig_microclim_Q1_gradient"), g1, 8, 6)

# ---- Figure 2: Q2 residual vs topo (the topo payload) ------------------------
F2 <- rbindlist(lapply(c("res_beta","res_dTmaxS"), function(rv)
  rbindlist(lapply(topo3, function(tv) data.table(residual=rv, topo=tv, x=D[[tv]], y=D[[rv]])))))
F2[, residual:=factor(residual, levels=c("res_beta","res_dTmaxS"),
                      labels=c("residual beta (obs-MuSICA)","residual dTmax summer (obs-MuSICA)"))]
lab2 <- Q2PAIR[residual %in% c("res_beta","res_dTmaxS")][, residual:=factor(residual,
          levels=c("res_beta","res_dTmaxS"), labels=levels(F2$residual))]
g2 <- ggplot(F2, aes(x,y)) + geom_hline(yintercept=0,linetype="dotted",colour="grey60") +
  geom_point(size=1.5, colour="#D55E00", alpha=.8) + geom_smooth(method="lm",se=TRUE,colour="#0072B2",linewidth=.6) +
  geom_text(data=lab2, aes(x=Inf,y=Inf,label=sprintf("r=%.2f  p=%.3f",pearson_r,p)),
            hjust=1.05, vjust=1.4, size=3, inherit.aes=FALSE) +
  facet_grid(residual ~ topo, scales="free") +
  labs(x="Topographic predictor", y="MuSICA residual (obs - sim)",
       subtitle="Q2 — does terrain explain what the 1-D column misses? residual = obs - MuSICA vs TWI / northness / slope. n=53.") +
  theme_article(10)
ggsave_article(file.path(OUT_FIG,"Fig_microclim_Q2_topo_residual"), g2, 9, 5.6)
cat("\nDONE -> Fig_microclim_Q1_gradient + Fig_microclim_Q2_topo_residual + 3 csv\n")
