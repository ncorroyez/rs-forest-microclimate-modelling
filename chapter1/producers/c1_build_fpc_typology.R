# ==============================================================================
# NEW typology with {LAI, Hmax, fCover, FPC1-3} (user request, 2026-06-30).
# ISOLATED — writes only NEW files, overwrites nothing. STAGE 1: typology + cLHS
# + a zero-sim PREVIEW (re-aggregate existing metrics6 sensitivities under the
# new labels) to see if the headline survives before launching 3200 new sims.
#   Rscript c1_build_fpc_typology.R
# Out (NEW): out_files/Sensitivity_Analysis/clhs_sample_fpc_v1.rds
#            out_files/Chapter1/tables/tab_fpc_typology_{centroids,preview}.csv
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(dplyr); library(purrr); library(cli) })
source(here::here("R/config.R")); source(here::here("R/forest.R"))

# ---- 1. forest + cached FPCA scores (same FPCA the article uses) -------------
fc  <- readRDS("outputs/figures_pipeline_z05/data/forest_fpca.rds")
fo  <- as.data.frame(fc$forest)
sco <- fc$fpca$fpca$scores
stopifnot(nrow(sco) == nrow(fo))
fo$FPC1 <- sco[,1]; fo$FPC2 <- sco[,2]; fo$FPC3 <- sco[,3]

CV <- c("LAI","Hmax","fCover","FPC1","FPC2","FPC3")
cat("\n== New typology: k-means on LAI+Hmax+fCover+FPC1-3 ==\n")
cl <- label_clusters(fo, k = 4, vars = CV)          # seed 42, scale, kmeans nstart=25
fo_cl <- cl$df
# relabel raw k-means code -> P1..P4 by ASCENDING LAI (keep density convention)
ord <- fo_cl %>% group_by(Cluster) %>% summarise(LAI=mean(LAI), .groups="drop") %>% arrange(LAI)
remap <- setNames(paste0("P", seq_len(nrow(ord))), ord$Cluster)
fo_cl$P <- remap[as.character(fo_cl$Cluster)]

cat("\n=== nouveaux archétypes FPC (centroïdes, ordre densité) ===\n")
cen <- as.data.table(fo_cl)[, .(n=.N, LAI=round(mean(LAI),2), Hmax=round(mean(Hmax),1),
        fCover=round(mean(fCover),2), VCI=round(mean(VCI),3),
        FPC1=round(mean(FPC1),2), FPC2=round(mean(FPC2),2)), by=P][order(P)]
print(cen)
fwrite(cen, "out_files/Chapter1/tables/tab_fpc_typology_centroids.csv")

# ---- 2. NEW cLHS sample (100/cluster), isolated -------------------------------
# sample_clhs_per_cluster keys on $Cluster; feed the new (raw) code, carry P + FPC
df_s <- sample_clhs_per_cluster(fo_cl, n_per_cluster = 100, vars = CV, iter = 10000)
df_s$P <- remap[as.character(df_s$Cluster)]
df_sf <- df_s; df_sf$fCover[df_sf$fCover < 0.5] <- 0.5
saveRDS(df_sf, "out_files/Sensitivity_Analysis/clhs_sample_fpc_v1.rds")
cli_alert_success("NEW sample: {nrow(df_sf)} plots -> clhs_sample_fpc_v1.rds (rien d'écrasé)")
cat("tailles par archétype:\n"); print(table(df_sf$P))

# ---- 3. ZERO-SIM PREVIEW: re-aggregate EXISTING sensitivities by NEW labels ---
# map the FROZEN 400 pixels (which already have sims) to their NEW FPC label,
# then recompute the per-archetype ΔTmax sensitivity (deltadelta NORM=sd logic).
froz <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
froz <- froz[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
froz[, pid := sprintf("S%04d", .I)]
fonew <- as.data.table(fo_cl)[, .(x, y, Pnew = P)]
froz <- merge(froz, fonew, by = c("x","y"), all.x = TRUE)        # new label per frozen pixel
cat(sprintf("\nfrozen pixels remappés au nouveau clustering: %d/%d\n", sum(!is.na(froz$Pnew)), nrow(froz)))

M <- rbindlist(lapply(list.files("out_files/Chapter1/tables/metrics6_v320","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
M <- M[metric=="Tmax_all"]
M <- merge(M, froz[, .(pid, Pnew)], by="pid", all.x=TRUE)
M[, `:=`(sl_LAI=(LAI_add-LAI_rem)/2, sl_fCover=(fCov_add-fCov_rem)/2, sl_Hmax=(Hmax_add-Hmax_rem)/2, sl_LAD=dT_LAD)]
M[, fCov_add:=ifelse(fCover>=0.95,0,fCov_add)][, fCov_rem:=ifelse(fCover<=0.5001,NA_real_,fCov_rem)]
# per +1 within-NEW-archetype SD
M[, `:=`(sLAI=sd(LAI/2,na.rm=T), sFC=sd(fCover,na.rm=T)/0.1, sHM=sd(Hmax,na.rm=T)/5), by=Pnew]
prev <- M[!is.na(Pnew), .(
  LAI    = round(median(sl_LAI*sLAI, na.rm=T),3),
  fCover = round(median(sl_fCover*sFC, na.rm=T),3),
  Hmax   = round(median(sl_Hmax*sHM, na.rm=T),3),
  n=.N), by=Pnew][order(Pnew)]
cat("\n=== PREVIEW (sims réutilisées, ré-agrégées par NOUVEAUX labels) — ΔTmax sensi par +1 SD ===\n")
print(prev)
fwrite(prev, "out_files/Chapter1/tables/tab_fpc_typology_preview.csv")
cat("\n→ compare au gelé (P1 LAI≈0.57 dominant → P4 LAI≈0.12). Le headline (LAI domine, érosion densité) tient-il ?\n")
cat("NB: preview = approximation (mêmes 400 pixels regroupés) ; le full re-run fera un vrai cLHS+sims sur clhs_sample_fpc_v1.rds.\n")
