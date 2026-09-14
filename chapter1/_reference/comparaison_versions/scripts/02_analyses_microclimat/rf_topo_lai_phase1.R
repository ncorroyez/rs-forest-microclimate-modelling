# ==============================================================================
# Phase 1 (decisive, cheap): does adding topography + full structure improve the
# RF prediction of LAI_ALS over the existing DYN_RF feature set? Compare
#   RF_base = {LAI_S2_ATBD, Hmax, fCover, VCI, LCV}        (= existing DYN_RF)
#   RF_topo = base + {CHM_std, slope, aspect_sin, aspect_cos, northness, TWI}
# under leave-one-out CV (n=53 plots), reporting LOO R²/RMSE and RF_topo variable
# importance. Only if topo materially lifts LOO R² is it worth running the 53
# MuSICA sims for a DYN_RF_TOPO scenario (terrain->canopy density is a weak link).
#   Rscript rf_topo_lai_phase1.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(data.table); library(randomForest); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

MET <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Deciduous_Only"
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
stopifnot(all(c("LAI_ALS","LAI_S2_ATBD","Hmax","fCover","VCI","LCV","x","y") %in% names(df)))

# ---- extract CHM_std + topo at plot locations --------------------------------
topo_files <- c(CHM_std="std_res_10_m.tif", slope="slope_res_10_m.tif",
                aspect_sin="aspect_sin_res_10_m.tif", aspect_cos="aspect_cos_res_10_m.tif",
                northness="northness_slope_res_10_m.tif", twi="twi_res_10_m.tif")
stk <- rast(file.path(MET, topo_files)); names(stk) <- names(topo_files)
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs=crs(stk))
ex  <- as.data.table(terra::extract(stk, pts))[, -1]
df  <- cbind(df, ex)
cat("Extracted topo (NA per var):\n"); print(colSums(is.na(ex)))

base  <- c("LAI_S2_ATBD","Hmax","fCover","VCI","LCV")
topo  <- c("CHM_std","slope","aspect_sin","aspect_cos","northness","twi")
feats_base <- base
feats_topo <- c(base, topo)

# impute residual NA by column mean (few cells off-raster)
for (v in unique(c(feats_topo))) if (anyNA(df[[v]])) df[[v]][is.na(df[[v]])] <- mean(df[[v]], na.rm=TRUE)
D <- as.data.frame(df[, c("LAI_ALS", feats_topo), with=FALSE])

# ---- leave-one-out CV --------------------------------------------------------
loo_pred <- function(feats, ntree=500){
  p <- numeric(nrow(D))
  for (i in seq_len(nrow(D))){
    set.seed(42)
    m <- randomForest(x=D[-i, feats, drop=FALSE], y=D$LAI_ALS[-i], ntree=ntree,
                      mtry=max(1, floor(length(feats)/3)))
    p[i] <- predict(m, newdata=D[i, feats, drop=FALSE])
  }
  pmax(p, 0)
}
r2  <- function(o,p) 1 - sum((o-p)^2)/sum((o-mean(o))^2)
rmse<- function(o,p) sqrt(mean((o-p)^2))
mae <- function(o,p) mean(abs(o-p))

p_base <- loo_pred(feats_base); p_topo <- loo_pred(feats_topo)
o <- D$LAI_ALS
res <- data.table(
  model = c("RF_base (DYN_RF features)","RF_topo (+CHM_std +topo)"),
  n_feat= c(length(feats_base), length(feats_topo)),
  LOO_R2  = c(r2(o,p_base),  r2(o,p_topo)),
  LOO_RMSE= c(rmse(o,p_base),rmse(o,p_topo)),
  LOO_MAE = c(mae(o,p_base), mae(o,p_topo)))
cat("\n=== LOO-CV: predicting LAI_ALS (n=53) ===\n")
print(res[, .(model, n_feat, LOO_R2=round(LOO_R2,3), LOO_RMSE=round(LOO_RMSE,3), LOO_MAE=round(LOO_MAE,3))])
cat(sprintf("\n  delta LOO_R2 (topo - base) = %+.3f | delta RMSE = %+.3f m2/m2\n",
            r2(o,p_topo)-r2(o,p_base), rmse(o,p_topo)-rmse(o,p_base)))

# also: S2 raw (no RF) and the in-pipeline LAI_RF as anchors
cat(sprintf("  anchor: raw S2_ATBD vs LAI_ALS  R2=%.3f RMSE=%.3f\n", r2(o,df$LAI_S2_ATBD), rmse(o,df$LAI_S2_ATBD)))

# ---- full-sample importance for RF_topo --------------------------------------
set.seed(42)
m_full <- randomForest(x=D[,feats_topo,drop=FALSE], y=D$LAI_ALS, ntree=1000, importance=TRUE,
                       mtry=max(1, floor(length(feats_topo)/3)))
imp <- as.data.table(importance(m_full), keep.rownames="var")[order(-`%IncMSE`)]
imp[, is_topo := var %in% topo]
cat("\n=== RF_topo variable importance (%IncMSE, desc) ===\n")
print(imp[, .(var, pct_IncMSE=round(`%IncMSE`,2), IncNodePurity=round(IncNodePurity,2), is_topo)])
cat(sprintf("\n  topo block summed %%IncMSE = %.1f / total %.1f  (%.0f%%)\n",
            sum(imp[is_topo==TRUE]$`%IncMSE`), sum(imp$`%IncMSE`),
            100*sum(imp[is_topo==TRUE]$`%IncMSE`)/sum(imp$`%IncMSE`)))

# save predictions for a possible Phase 2 (MuSICA scenario)
out <- df[, .(plot_id, x, y, LAI_ALS, LAI_S2_ATBD, LAI_RF_base_loo=p_base, LAI_RF_topo_loo=p_topo)]
fwrite(out, "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/tables/rf_topo_lai_loo.csv")
fwrite(res, "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/tables/rf_topo_lai_phase1_metrics.csv")
cat("\nDONE -> rf_topo_lai_loo.csv + rf_topo_lai_phase1_metrics.csv\n")
