# ==============================================================================
# FPC-typology PREVIEW, completed: add the LAD (real-vs-uniform) co-lead + the
# bootstrap of the LAI-LAD difference, under the NEW {LAI,Hmax,fCover,FPC1-3}
# labels, re-using EXISTING sims (zero new sims). Answers: does the dense co-lead
# of the vertical profile survive the FPC typology?
#   Rscript c1_fpc_preview_lad_bootstrap.R
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(dplyr); library(mclust); source(here::here("R/forest.R")) })
set.seed(1); NB <- 3000

# ---- new FPC labels per forest pixel (same recipe as c1_build_fpc_typology) --
fc <- readRDS("outputs/figures_pipeline_z05/data/forest_fpca.rds")
fo <- as.data.frame(fc$forest); sco <- fc$fpca$fpca$scores
fo$FPC1<-sco[,1]; fo$FPC2<-sco[,2]; fo$FPC3<-sco[,3]
cl <- label_clusters(fo, k=4, vars=c("LAI","Hmax","fCover","FPC1","FPC2","FPC3"))$df
ord <- as.data.table(cl)[, .(LAI=mean(LAI)), by=Cluster][order(LAI)]
remap <- setNames(paste0("P",seq_len(nrow(ord))), ord$Cluster)
cl$P <- remap[as.character(cl$Cluster)]
fonew <- as.data.table(cl)[, .(x,y,Pnew=P)]

# ---- frozen 400 (have sims) -> new label, + VCI -------------------------------
froz <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
froz <- froz[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]; froz[, pid:=sprintf("S%04d",.I)]
froz <- merge(froz, fonew, by=c("x","y"), all.x=TRUE)

M <- rbindlist(lapply(list.files("out_files/Chapter1/tables/metrics6_v320","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
M <- M[metric=="Tmax_all"]
M <- merge(M, froz[, .(pid, Pnew, VCI)], by="pid", all.x=TRUE)
M <- M[!is.na(Pnew)]
M[, fCov_add:=ifelse(fCover>=0.95,0,fCov_add)][, fCov_rem:=ifelse(fCover<=0.5001,NA_real_,fCov_rem)]

# ---- signed per +1 within-archetype SD (incl. LAD per SD of (1-VCI)) ----------
sens <- M[, {
  sLAI<-sd(LAI/2,na.rm=T); sFC<-sd(fCover,na.rm=T)/0.1; sHM<-sd(Hmax,na.rm=T)/5; sVCI<-sd(VCI,na.rm=T)
  .(LAI   = round(median((LAI_add-LAI_rem)/2 * sLAI, na.rm=T),3),
    fCover= round(median((fCov_add-fCov_rem)/2 * sFC, na.rm=T),3),
    Hmax  = round(median((Hmax_add-Hmax_rem)/2 * sHM, na.rm=T),3),
    LAD   = round(median(dT_LAD/pmax(1-VCI,0.05) * sVCI, na.rm=T),3), n=.N)
}, by=Pnew][order(Pnew)]
cat("=== NEW FPC typology — ΔTmax sensitivity per +1 SD (signed; LAD = real-vs-uniform) ===\n"); print(sens)

# ---- bootstrap LAI - LAD importance difference per new archetype --------------
imp_of <- function(d){
  sL <- median(abs((d$LAI_add-d$LAI_rem)/2),na.rm=T)*sd(d$LAI,na.rm=T)/2
  c(LAI=sL, LAD=abs(median(d$dT_LAD,na.rm=T)))
}
ci <- function(v) c(med=median(v), lo=quantile(v,.025,names=F), hi=quantile(v,.975,names=F))
res <- rbindlist(lapply(c("P1","P2","P3","P4"), function(p){
  d <- M[Pnew==p]; n <- nrow(d)
  B <- replicate(NB, { i<-sample(n,n,TRUE); v<-imp_of(d[i]); unname(v["LAI"]-v["LAD"]) })
  c0 <- ci(B); pt <- imp_of(d)
  data.table(P=p, n=n, impLAI=round(pt["LAI"],3), impLAD=round(pt["LAD"],3),
             LAI_minus_LAD=round(c0[1],3), lo=round(c0[2],3), hi=round(c0[3],3),
             straddle0 = c0[2]<=0 & c0[3]>=0)
}))
cat("\n=== bootstrap co-lead (LAI-LAD, 95% CI) under NEW labels ===\n"); print(res)
cat("\n→ Frozen: LAI-LAD straddled 0 ONLY in P4 (co-lead). Holds under FPC typology?\n")
fwrite(sens, "out_files/Chapter1/tables/tab_fpc_preview_sensitivity.csv")
fwrite(res,  "out_files/Chapter1/tables/tab_fpc_preview_bootstrap.csv")
