# ==============================================================================
# Chapter 3 — recompute every manuscript ΔTmax statistic under the OLD (own-max)
# and NEW (Ch1 time-matched) conventions, from perplot_dtmax_conventions.csv.
# Slope-of-coupling metrics are untouched by design.
# Writes: NC_Full/manuscripts/ch3/tables/timematched_reextraction.csv (+ sweep)
#   Rscript c3_tm_analysis.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table)})
source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
PP<-fread(file.path(CFG_C3$out_dir,"tables","perplot_dtmax_conventions.csv"))
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
CROSS<-3.86; HGATE<-17.8
df[,st:=ifelse(LAI_ALS<CROSS,"open","dense")]
OBS<-PP[scenario=="OBS",.(id_plot,do_old=d_old,do_new=d_new,no_old=n_old,no_new=n_new)]
key<-df[,.(id_plot,LAI_ALS,FORMS_H,Hmax,st)]
OBS<-merge(OBS,key,by="id_plot")
R2<-function(a,b)suppressWarnings(cor(a,b,use="complete.obs")^2)
RM<-function(a,b)sqrt(mean((a-b)^2,na.rm=TRUE))
rows<-list(); addrow<-function(metric,scenario,old,new)
  rows[[length(rows)+1]]<<-data.table(metric=metric,scenario=scenario,
    old=round(old,4),new=round(new,4))

# ---- per-scenario table (day) ------------------------------------------------
SCN<-setdiff(unique(PP$scenario),"OBS")
mk<-function(sv,ov,st){D<-st=="dense";O<-st=="open"
  c(dense=R2(sv[D],ov[D]),open=R2(sv[O],ov[O]),pooled=R2(sv,ov),
    bias=mean(sv-ov,na.rm=TRUE),RMSE=RM(sv,ov),SDrec=sd(sv,na.rm=TRUE)/sd(ov,na.rm=TRUE),
    maximin=min(R2(sv[D],ov[D]),R2(sv[O],ov[O])))}
allstrat<-rbindlist(lapply(SCN,function(s){
  M<-merge(OBS,PP[scenario==s,.(id_plot,ds_old=d_old,ds_new=d_new,ns_old=n_old,ns_new=n_new)],by="id_plot")
  o<-mk(M$ds_old,M$do_old,M$st); n<-mk(M$ds_new,M$do_new,M$st)
  cbind(data.table(scenario=s,conv=c("old","new")),rbindlist(list(as.list(o),as.list(n))))}))
fwrite(allstrat,file.path(TAB,"timematched_all_stratified.csv"))
for(s in SCN){a<-allstrat[scenario==s&conv=="old"];b<-allstrat[scenario==s&conv=="new"]
  for(m in c("dense","open","pooled","bias","RMSE","SDrec","maximin"))
    addrow(paste0("dT_",m),s,a[[m]],b[[m]])}

# ---- composites (switch / blend recomposition / oracle) -----------------------
L<-PP[scenario=="STATIC_ALS",.(id_plot,aL_old=d_old,aL_new=d_new)]
S<-PP[scenario=="STATIC_S2_ATBD",.(id_plot,aS_old=d_old,aS_new=d_new)]
M<-Reduce(function(x,y)merge(x,y,by="id_plot"),list(OBS,L,S))
comp<-function(sv,ov,st,lab,conv){D<-st=="dense";O<-st=="open"
  deb<-sv-ave(sv-ov,st)  # per-stratum de-biased
  data.table(metric=c("dT_dense","dT_open","dT_pooled","dT_bias","dT_RMSE","dT_RMSE_debias","dT_SDrec"),
    scenario=lab,conv=conv,
    value=round(c(R2(sv[D],ov[D]),R2(sv[O],ov[O]),R2(sv,ov),mean(sv-ov),RM(sv,ov),RM(deb,ov),sd(sv)/sd(ov)),4))}
CC<-list()
for(cv in c("old","new")){
  a<-if(cv=="old")M$aL_old else M$aL_new; s<-if(cv=="old")M$aS_old else M$aS_new
  o<-if(cv=="old")M$do_old else M$do_new
  CC[[paste0("swL_",cv)]] <-comp(ifelse(M$LAI_ALS>=CROSS,a,s),o,M$st,"Switch LAI-keyed",cv)
  CC[[paste0("swH_",cv)]] <-comp(ifelse(M$FORMS_H>=HGATE,a,s),o,M$st,"Switch FORMS-H",cv)
  CC[[paste0("orc_",cv)]] <-comp(ifelse(abs(a-o)<=abs(s-o),a,s),o,M$st,"Oracle per-plot",cv)}
CO<-rbindlist(CC)
CW<-dcast(CO,metric+scenario~conv,value.var="value")
for(i in seq_len(nrow(CW))) addrow(CW$metric[i],CW$scenario[i],CW$old[i],CW$new[i])

# ---- switch at k=0.5 (figQ) ------------------------------------------------------
K5<-PP[scenario=="STATIC_ALS_K05",.(id_plot,k_old=d_old,k_new=d_new)]
M5<-merge(M,K5,by="id_plot")
for(cv in c("old","new")){
  a<-if(cv=="old")M5$k_old else M5$k_new; s<-if(cv=="old")M5$aS_old else M5$aS_new
  o<-if(cv=="old")M5$do_old else M5$do_new
  CC[[paste0("sw5_",cv)]]<-comp(ifelse(M5$LAI_ALS>=CROSS,a,s),o,M5$st,"Switch k=0.5",cv)}
CO<-rbindlist(CC[grepl("sw5_",names(CC))])
CW5<-dcast(CO,metric+scenario~conv,value.var="value")
for(i in seq_len(nrow(CW5))) addrow(CW5$metric[i],CW5$scenario[i],CW5$old[i],CW5$new[i])

# ---- multiclass k-means selection (figT) ------------------------------------------
Xk<-scale(as.matrix(df[match(M$id_plot,df$id_plot),.(LAI_ALS,Hmax,fCover,VCI,LCV)]))
for(k in c(3,4)){set.seed(1);cl<-kmeans(Xk,centers=k,nstart=25)$cluster
  for(cv in c("old","new")){
    a<-if(cv=="old")M$aL_old else M$aL_new; s<-if(cv=="old")M$aS_old else M$aS_new
    o<-if(cv=="old")M$do_old else M$do_new
    pick<-sapply(1:k,function(c){ix<-cl==c;if(mean((a[ix]-o[ix])^2)<=mean((s[ix]-o[ix])^2))"L" else "S"})
    dv<-ifelse(pick[cl]=="L",a,s)
    if(cv=="old"){to_old<-c(R2(dv,o),R2(dv[M$st=="dense"],o[M$st=="dense"]),R2(dv[M$st=="open"],o[M$st=="open"]))}
    else{to_new<-c(R2(dv,o),R2(dv[M$st=="dense"],o[M$st=="dense"]),R2(dv[M$st=="open"],o[M$st=="open"]))}}
  addrow("multiclass_pooled",sprintf("kmeans_k%d",k),to_old[1],to_new[1])
  addrow("multiclass_dense",sprintf("kmeans_k%d",k),to_old[2],to_new[2])
  addrow("multiclass_open",sprintf("kmeans_k%d",k),to_old[3],to_new[3])}

# ---- threshold sweep (Fig 5) ---------------------------------------------------
sweep<-rbindlist(lapply(seq(2.0,5.5,0.25),function(t){
  de<-M$LAI_ALS>=t
  data.table(thr=t,n_dense=sum(de),n_open=sum(!de),
    LiDAR_dense_old=R2(M$aL_old[de],M$do_old[de]), LiDAR_dense_new=R2(M$aL_new[de],M$do_new[de]),
    S2_open_old=R2(M$aS_old[!de],M$do_old[!de]),   S2_open_new=R2(M$aS_new[!de],M$do_new[!de]),
    S2_dense_old=R2(M$aS_old[de],M$do_old[de]),    S2_dense_new=R2(M$aS_new[de],M$do_new[de]),
    LiDAR_open_old=R2(M$aL_old[!de],M$do_old[!de]),LiDAR_open_new=R2(M$aL_new[!de],M$do_new[!de]))}))
fwrite(sweep,file.path(TAB,"timematched_threshold_sweep.csv"))
for(i in seq_len(nrow(sweep))){t<-sweep$thr[i]
  addrow("sweep_LiDAR_dense_R2",sprintf("thr=%.2f",t),sweep$LiDAR_dense_old[i],sweep$LiDAR_dense_new[i])
  addrow("sweep_S2_open_R2",   sprintf("thr=%.2f",t),sweep$S2_open_old[i],sweep$S2_open_new[i])
  addrow("sweep_S2_dense_R2",  sprintf("thr=%.2f",t),sweep$S2_dense_old[i],sweep$S2_dense_new[i])
  addrow("sweep_LiDAR_open_R2",sprintf("thr=%.2f",t),sweep$LiDAR_open_old[i],sweep$LiDAR_open_new[i])}

# ---- segmented breakpoint of observed dTmax ~ LAI ------------------------------
bp<-function(y,x){ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok]
  cand<-seq(quantile(x,.2),quantile(x,.8),length.out=40)
  rss<-sapply(cand,function(c){xl<-pmin(x,c);xr<-pmax(x-c,0);sum(resid(lm(y~xl+xr))^2)})
  cand[which.min(rss)]}
addrow("breakpoint_obs_dTmax_vs_LAI","observed",bp(M$do_old,M$LAI_ALS),bp(M$do_new,M$LAI_ALS))

# ---- paired ΔR2 bootstrap A (dense): LiDAR - S2 ---------------------------------
D<-M[st=="dense"]; nD<-nrow(D)
bootA<-function(a,s,o){set.seed(42);pt<-R2(a,o)-R2(s,o)
  v<-replicate(2000,{i<-sample(nD,nD,replace=TRUE);R2(a[i],o[i])-R2(s[i],o[i])})
  c(pt,quantile(v,c(.025,.975),na.rm=TRUE))}
ao<-bootA(D$aL_old,D$aS_old,D$do_old); an<-bootA(D$aL_new,D$aS_new,D$do_new)
addrow("bootA_deltaR2_dense_LiDAR_minus_S2","point",ao[1],an[1])
addrow("bootA_deltaR2_dense_LiDAR_minus_S2","ci_lo",ao[2],an[2])
addrow("bootA_deltaR2_dense_LiDAR_minus_S2","ci_hi",ao[3],an[3])

# ---- per-stratum bootstrap CIs (LiDAR & S2, dense/open/pooled) --------------------
bootR2<-function(s,o,B=2000,seed=42){set.seed(seed);n<-length(s)
  v<-replicate(B,{i<-sample(n,n,replace=TRUE);R2(s[i],o[i])})
  quantile(v,c(.025,.975),na.rm=TRUE)}
for(sens in c("LiDAR","S2")){
  for(strat in c("dense","open","pooled")){
    ix<-if(strat=="pooled")rep(TRUE,nrow(M)) else M$st==strat
    so_<-if(sens=="LiDAR")list(M$aL_old[ix],M$aL_new[ix]) else list(M$aS_old[ix],M$aS_new[ix])
    co<-bootR2(so_[[1]],M$do_old[ix]); cn<-bootR2(so_[[2]],M$do_new[ix])
    addrow(paste0("bootCI_R2_",strat,"_lo"),sens,co[1],cn[1])
    addrow(paste0("bootCI_R2_",strat,"_hi"),sens,co[2],cn[2])}}

# ---- pooled bootstrap: switch − naive, switch − S2 ---------------------------------
NV<-PP[scenario=="NAIVE_S2_FORMSH",.(id_plot,nv_old=d_old,nv_new=d_new)]
MB<-merge(M,NV,by="id_plot")
bootD<-function(a,b,o,B=2000,seed=42){set.seed(seed);n<-length(a)
  pt<-R2(a,o)-R2(b,o)
  v<-replicate(B,{i<-sample(n,n,replace=TRUE);R2(a[i],o[i])-R2(b[i],o[i])})
  c(pt,quantile(v,c(.025,.975),na.rm=TRUE))}
for(cv in c("old","new")){
  aa<-if(cv=="old")MB$aL_old else MB$aL_new; ss<-if(cv=="old")MB$aS_old else MB$aS_new
  oo<-if(cv=="old")MB$do_old else MB$do_new; nv<-if(cv=="old")MB$nv_old else MB$nv_new
  sw<-ifelse(MB$LAI_ALS>=CROSS,aa,ss)
  r1<-bootD(sw,nv,oo); r2v<-bootD(sw,ss,oo)
  if(cv=="old"){sn_o<-r1;s2_o<-r2v}else{sn_n<-r1;s2_n<-r2v}}
addrow("bootPooled_switch_minus_naive","point",sn_o[1],sn_n[1])
addrow("bootPooled_switch_minus_naive","ci_lo",sn_o[2],sn_n[2])
addrow("bootPooled_switch_minus_naive","ci_hi",sn_o[3],sn_n[3])
addrow("bootPooled_switch_minus_S2","point",s2_o[1],s2_n[1])
addrow("bootPooled_switch_minus_S2","ci_lo",s2_o[2],s2_n[2])
addrow("bootPooled_switch_minus_S2","ci_hi",s2_o[3],s2_n[3])

# ---- LAI-space gradient correlations (obs & MuSICA vs LiDAR LAI) --------------------
addrow("gradient_r_obs_dTmax_vs_LAI","pooled",cor(M$do_old,M$LAI_ALS),cor(M$do_new,M$LAI_ALS))
addrow("gradient_r_obs_dTmax_vs_LAI","dense",cor(M$do_old[M$st=="dense"],M$LAI_ALS[M$st=="dense"]),
       cor(M$do_new[M$st=="dense"],M$LAI_ALS[M$st=="dense"]))
addrow("gradient_r_sim_dTmax_vs_LAI","STATIC_ALS",cor(M$aL_old,M$LAI_ALS),cor(M$aL_new,M$LAI_ALS))
addrow("gradient_r_FORMSH_vs_obs_dTmax","pooled",cor(M$FORMS_H,M$do_old),cor(M$FORMS_H,M$do_new))

# ---- night (ΔTmin) --------------------------------------------------------------
for(s in c("STATIC_ALS","STATIC_S2_ATBD","FUSION_H","BLEND_H","NAIVE_S2_FORMSH")){
  Mn<-merge(OBS,PP[scenario==s,.(id_plot,ns_old=n_old,ns_new=n_new)],by="id_plot")
  o<-mk(Mn$ns_old,Mn$no_old,Mn$st); n<-mk(Mn$ns_new,Mn$no_new,Mn$st)
  for(m in c("dense","open","pooled","bias"))
    addrow(paste0("night_dTmin_",m),s,o[[m]],n[[m]])}

# ---- opt in-domain (LAI_ALS>2 & Hmax>10, n=39) ----------------------------------
OBS[,inDom:=LAI_ALS>2&Hmax>10]
for(s in c("STATIC_S2_OPT","STATIC_S2_ATBD","STATIC_ALS","BLEND_H")){
  Md<-merge(OBS,PP[scenario==s,.(id_plot,ds_old=d_old,ds_new=d_new)],by="id_plot")[inDom==TRUE]
  addrow("dT_pooled_inDomain",s,R2(Md$ds_old,Md$do_old),R2(Md$ds_new,Md$do_new))
  addrow("dT_dense_inDomain",s,R2(Md$ds_old[Md$st=="dense"],Md$do_old[Md$st=="dense"]),
                               R2(Md$ds_new[Md$st=="dense"],Md$do_new[Md$st=="dense"]))
  addrow("dT_bias_inDomain",s,mean(Md$ds_old-Md$do_old),mean(Md$ds_new-Md$do_new))}

# ---- naive dense bias (Fig 7 talking point) --------------------------------------
Mn<-merge(OBS,PP[scenario=="NAIVE_S2_FORMSH",.(id_plot,d_old,d_new)],by="id_plot")
addrow("dT_bias_dense","NAIVE_S2_FORMSH",
  mean(Mn[st=="dense",d_old-do_old]),mean(Mn[st=="dense",d_new-do_new]))
addrow("dT_bias_open","NAIVE_S2_FORMSH",
  mean(Mn[st=="open",d_old-do_old]),mean(Mn[st=="open",d_new-do_new]))

# ---- obs description --------------------------------------------------------------
addrow("obs_mean_dTmax","OBS",mean(OBS$do_old),mean(OBS$do_new))
addrow("obs_sd_dTmax","OBS",sd(OBS$do_old),sd(OBS$do_new))

OUTT<-rbindlist(rows)
fwrite(OUTT,file.path(TAB,"timematched_reextraction.csv"))
cat(sprintf("wrote %s (%d rows)\n",file.path(TAB,"timematched_reextraction.csv"),nrow(OUTT)))

# ---- console digest ---------------------------------------------------------------
dig<-function(s)print(allstrat[scenario==s,.(conv,dense=round(dense,3),open=round(open,3),
  pooled=round(pooled,3),bias=round(bias,3),RMSE=round(RMSE,3),SDrec=round(SDrec,3))],row.names=FALSE)
cat("\n=== STATIC_ALS ===\n");dig("STATIC_ALS")
cat("\n=== STATIC_S2_ATBD ===\n");dig("STATIC_S2_ATBD")
cat("\n=== STATIC_ALS_DOPT ===\n");dig("STATIC_ALS_DOPT")
cat("\n=== NAIVE_S2_FORMSH ===\n");dig("NAIVE_S2_FORMSH")
cat("\n=== FUSION_H (dir) ===\n");dig("FUSION_H")
cat("\n=== BLEND_H ===\n");dig("BLEND_H")
cat("\ncomposites:\n");print(CW,row.names=FALSE)
cat(sprintf("\nbreakpoint obs dTmax~LAI: old=%.2f new=%.2f\n",
  bp(M$do_old,M$LAI_ALS),bp(M$do_new,M$LAI_ALS)))
cat("\nsweep (reversal check):\n");print(sweep[,.(thr,LiDAR_dense_new=round(LiDAR_dense_new,2),
  S2_dense_new=round(S2_dense_new,2),S2_open_new=round(S2_open_new,2),
  LiDAR_open_new=round(LiDAR_open_new,2))],row.names=FALSE)
