# Chapter 3 @2×LAI — significance of the layered-hybrid scenarios.
# Per-plot RMSE + paired bootstrap ΔRMSE, and block-bootstrap ΔR², fullyear,
# vs DYN_S2_ANNUAL and vs CONST_ALS. Uses cached 2×LAI nc.
# Run from z_Example root:  Rscript c3_hybrid_signif.R
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot

# rebuild hybrid scenarios (same as c3_hybrid_musica.R)
key <- df %>% mutate(pid=sprintf("X%d_Y%d",round(x),round(y))) %>%
  transmute(pid, LAI_ALS, dopt=LAI_ALS_DOPT, below=pmax(LAI_ALS-LAI_ALS_DOPT,0))
stretch_frac <- function(doy,frac,s){approx(doy, frac, xout=200+(doy-200)/s, rule=2)$y}
build_hyb <- function(s){pids<-intersect(names(tb$annual),key$pid)
  setNames(lapply(pids,function(p){a<-tb$annual[[p]];k<-key[key$pid==p,]
    ft<-pmin(pmax(a$lai/k$LAI_ALS,0),1);fb<-pmin(pmax(stretch_frac(a$doy,ft,s),0),1)
    data.frame(doy=a$doy, lai=k$dopt*ft + k$below*fb)}),pids)}
.fn<-function(c)function(pr)as.numeric(pr[[c]])
all_sc <- make_all_scenarios_c3(tb, CFG_C3$list_year, CFG_C3$d_opt_m, CFG_C3$scenarios_mode)
sc <- list(DYN_S2_ANNUAL=all_sc[["DYN_S2_ANNUAL"]], CONST_ALS=all_sc[["CONST_ALS"]],
           STATIC_S2_ATBD=all_sc[["STATIC_S2_ATBD"]], STATIC_ALS=all_sc[["STATIC_ALS"]],
           DYN_RF=all_sc[["DYN_RF"]],
           DYN_HYBRID_s10=list(name="DYN_HYBRID_s10",lai_fn=.fn("LAI_ALS"),hmax_fn=.fn("Hmax"),fcover_fn=.fn("fCover"),lad_fn=make_lad_real,phenology_fn=make_phenology_fn_factory(build_hyb(1.0),CFG_C3$list_year)),
           DYN_HYBRID_s12=list(name="DYN_HYBRID_s12",lai_fn=.fn("LAI_ALS"),hmax_fn=.fn("Hmax"),fcover_fn=.fn("fCover"),lad_fn=make_lad_real,phenology_fn=make_phenology_fn_factory(build_hyb(1.2),CFG_C3$list_year)))
sc <- Filter(Negate(is.null), sc)
ds <- seq(as.Date("2021-01-01"),as.Date("2021-12-31"),by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds); hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
val <- validate_scenarios_at_hobos(df, hd, sc, file.path(CFG_C3$out_dir,"nc"), dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
d <- val$daily
pr <- d %>% group_by(scenario,id_plot) %>% summarise(rmse=sqrt(mean((Delta_sim-Delta_obs)^2,na.rm=T)),.groups="drop")
W <- pr %>% pivot_wider(names_from=scenario,values_from=rmse); W <- W[complete.cases(W),]
plots <- unique(d$id_plot); byp <- split(d, d$id_plot)
pooled_r2 <- function(x){x<-x[is.finite(x$Delta_obs)&is.finite(x$Delta_sim),];if(nrow(x)<3)NA else cor(x$Delta_obs,x$Delta_sim)^2}
set.seed(7); B<-3000
cat(sprintf("n=%d plots\n\n=== fullyear, vs each REFERENCE: ΔRMSE [paired boot] & ΔR² [block boot] ===\n", nrow(W)))
for (REF in c("DYN_S2_ANNUAL","CONST_ALS")) {
  cat(sprintf("\n--- vs %s ---\n", REF))
  for (s in setdiff(names(sc), REF)) {
    dl <- W[[s]]-W[[REF]]; ci <- quantile(replicate(B,mean(sample(dl,replace=T))),c(.025,.975))
    bd <- replicate(B,{sm<-sample(plots,replace=T);dd<-bind_rows(lapply(sm,function(p)byp[[p]]))
      pooled_r2(dd[dd$scenario==s,])-pooled_r2(dd[dd$scenario==REF,])}); cr<-quantile(bd,c(.025,.975),na.rm=T)
    cat(sprintf("  %-15s ΔRMSE=%+.3f[%+.3f,%+.3f]%s | ΔR²=%+.3f[%+.3f,%+.3f]%s | RMSE=%.3f R²=%.3f\n",
      s, mean(dl),ci[1],ci[2], if(ci[1]>0|ci[2]<0)"*"else" ", mean(bd,na.rm=T),cr[1],cr[2], if(cr[1]>0|cr[2]<0)"*"else" ",
      mean(W[[s]]), pooled_r2(d[d$scenario==s,]))) }
}
write.csv(W, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4d_hybrid_signif_2x.csv", row.names=FALSE)
cat("\nDONE\n")
