# ==============================================================================
# Append remaining old->new rows (RF, ablation, LAI-space, switch variants,
# blend debias, SDrec ceiling, wind-corr) to timematched_reextraction.csv.
#   Rscript c3_tm_assemble.R
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
T<-fread(file.path(TAB,"timematched_reextraction.csv"))
PP<-fread(file.path(CFG_C3$out_dir,"tables","perplot_dtmax_conventions.csv"))
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
CROSS<-3.86; df[,st:=ifelse(LAI_ALS<CROSS,"open","dense")]
OBS<-merge(PP[scenario=="OBS",.(id_plot,do_old=d_old,do_new=d_new)],
           df[,.(id_plot,LAI_ALS,st)],by="id_plot")
R2<-function(a,b)suppressWarnings(cor(a,b,use="complete.obs")^2)
rows<-list(); addrow<-function(metric,scenario,old,new)
  rows[[length(rows)+1]]<<-data.table(metric=metric,scenario=scenario,
    old=round(old,4),new=round(new,4))

# raw-bands RF
RB<-fread(file.path(TAB,"timematched_rawbands_rf.csv"))
addrow("rawbandsRF_pooled_R2","bands10",RB[conv=="old",pooled_r2],RB[conv=="new",pooled_r2])
addrow("rawbandsRF_dense_R2","bands10",RB[conv=="old",dense_r2],RB[conv=="new",dense_r2])
addrow("rawbandsRF_dense_ci_lo","bands10",RB[conv=="old",dense_ci_lo],RB[conv=="new",dense_ci_lo])
addrow("rawbandsRF_dense_ci_hi","bands10",RB[conv=="old",dense_ci_hi],RB[conv=="new",dense_ci_hi])
addrow("rawbandsRF_open_R2","bands10",RB[conv=="old",open_r2],RB[conv=="new",open_r2])
addrow("bootB_deltaR2_dense_rawbands_minus_S2","point",RB[conv=="old",B_deltaR2],RB[conv=="new",B_deltaR2])
addrow("bootB_deltaR2_dense_rawbands_minus_S2","ci_lo",RB[conv=="old",B_ci_lo],RB[conv=="new",B_ci_lo])
addrow("bootB_deltaR2_dense_rawbands_minus_S2","ci_hi",RB[conv=="old",B_ci_hi],RB[conv=="new",B_ci_hi])

# ablation
AB<-fread(file.path(TAB,"timematched_ablation.csv"))
for(i in seq_len(nrow(AB))){
  addrow("ablation_LOO_R2_LAI",AB$feature_set[i],AB$LOO_R2_LAI[i],AB$LOO_R2_LAI[i])
  addrow("ablation_cor_dTmax",AB$feature_set[i],AB$cor_dTmax_old[i],AB$cor_dTmax_new[i])}
addrow("ablation_deltaR2_S2FORMSH_minus_FORMSH","LAI target",
  AB[feature_set=="S2 + FORMS-H",LOO_R2_LAI]-AB[feature_set=="FORMS-H only",LOO_R2_LAI],
  AB[feature_set=="S2 + FORMS-H",LOO_R2_LAI]-AB[feature_set=="FORMS-H only",LOO_R2_LAI])

# LAI-space
LS<-fread(file.path(TAB,"timematched_laispace_cor.csv"))
for(i in seq_len(nrow(LS)))
  addrow("laispace_cor_dTmax",LS$product[i],LS$cor_dTmax_old[i],LS$cor_dTmax_new[i])

# switch variants (dense=LiDAR, open=X)
L<-PP[scenario=="STATIC_ALS",.(id_plot,aL_old=d_old,aL_new=d_new)]
for(sv in c("STATIC_S2_DOPT","STATIC_S2_RESCALED")){
  S<-PP[scenario==sv,.(id_plot,s_old=d_old,s_new=d_new)]
  M<-Reduce(function(x,y)merge(x,y,by="id_plot"),list(OBS,L,S))
  po<-mapply(function(cv){a<-M[[paste0("aL_",cv)]];s<-M[[paste0("s_",cv)]];o<-M[[paste0("do_",cv)]]
    sw<-ifelse(M$LAI_ALS>=CROSS,a,s);c(R2(sw,o),R2(sw[M$st=="open"],o[M$st=="open"]))},c("old","new"))
  addrow("switch_variant_pooled_R2",paste0("open=",sv),po[1,1],po[1,2])
  addrow("switch_variant_open_R2",paste0("open=",sv),po[2,1],po[2,2])}

# blend (BLEND_H dir) debiased RMSE
B<-merge(OBS,PP[scenario=="BLEND_H",.(id_plot,b_old=d_old,b_new=d_new)],by="id_plot")
deb<-function(sv,ov,st){d<-sv-ave(sv-ov,st);sqrt(mean((d-ov)^2))}
addrow("dT_RMSE_debias","BLEND_H",deb(B$b_old,B$do_old,B$st),deb(B$b_new,B$do_new,B$st))

# SDrec ceiling across scenario family (excl. k=0.5 rescale)
A<-fread(file.path(TAB,"timematched_all_stratified.csv"))
excl<-c("STATIC_ALS_K05")
addrow("max_SDrec_across_scenarios","excl K05",
  A[conv=="old"&!scenario%in%excl,max(SDrec)],A[conv=="new"&!scenario%in%excl,max(SDrec)])

# wind-corr robustness
WCR<-fread(file.path(TAB,"timematched_windcorr_robustness.csv"))
addrow("windcorr_bias_shift","STATIC_ALS",
  WCR[scenario=="STATIC_ALS",bias_old_wc-bias_old_nwc],WCR[scenario=="STATIC_ALS",bias_new_wc-bias_new_nwc])
addrow("windcorr_spearman_bias","across scenarios",
  cor(WCR$bias_old_nwc,WCR$bias_old_wc,method="spearman"),
  cor(WCR$bias_new_nwc,WCR$bias_new_wc,method="spearman"))
addrow("windcorr_spearman_SDrec","across scenarios",
  cor(WCR$SDrec_old_nwc,WCR$SDrec_old_wc,method="spearman"),
  cor(WCR$SDrec_new_nwc,WCR$SDrec_new_wc,method="spearman"))

OUT<-rbind(T,rbindlist(rows))
fwrite(OUT,file.path(TAB,"timematched_reextraction.csv"))
cat(sprintf("final timematched_reextraction.csv: %d rows\n",nrow(OUT)))
print(rbindlist(rows),nrows=60)
