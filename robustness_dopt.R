# Robustness of the d_opt out-of-sample gain (cLHS-train -> HOBO-test).
# 50 RF seeds x 4 models incl. a random-noise negative control, + bootstrap CI
# on the structure-vs-(+d_opt) RMSE gain. No MuSICA.
suppressPackageStartupMessages({ library(terra); library(randomForest) })
NCF <- path.expand("~/Documents/NC_Full")
r_s2   <- rast(file.path(NCF,"output/intermediate/sm6/Blois/s2lai_summer_atbd_T_res_10_m.tif"))
r_chm  <- rast(file.path(NCF,"output/intermediate/sm6/Blois/chm_sd_block_20_m.tif"))
r_dopt <- rast(file.path(NCF,"output/intermediate/lai_als_dopt/Blois/LAI_ALS_dopt_common.tif"))
ex <- function(r,x,y) as.numeric(terra::extract(r, vect(cbind(x,y),crs="EPSG:32631"))[,2])

cl <- as.data.frame(readRDS("out_files/Sensitivity_Analysis/clhs_sample.rds"))
tr <- data.frame(LAI_ALS=cl$LAI, Hmax=cl$Hmax, fCover=cl$fCover, VCI=cl$VCI,
                 LAI_S2_ATBD=ex(r_s2,cl$x,cl$y), CHM_std=ex(r_chm,cl$x,cl$y),
                 LAI_ALS_DOPT=ex(r_dopt,cl$x,cl$y))
tr <- tr[complete.cases(tr),]
hb <- readRDS("out_files/Chapter3/lai_prep/df_plots_lai.rds")
te <- data.frame(LAI_ALS=hb$LAI_ALS, Hmax=hb$Hmax, fCover=hb$fCover, VCI=hb$VCI,
                 LAI_S2_ATBD=hb$LAI_S2_ATBD, LAI_ALS_DOPT=hb$LAI_ALS_DOPT,
                 CHM_std=ex(r_chm,hb$x,hb$y))
te <- te[complete.cases(te),]
cat(sprintf("train=%d test=%d\n", nrow(tr), nrow(te)))

base <- c("LAI_S2_ATBD","Hmax","fCover","VCI")
feats <- list(structure=base, dopt=c(base,"LAI_ALS_DOPT"),
              het=c(base,"CHM_std"), noise=c(base,"NOISE"))
rmse <- function(p,o) sqrt(mean((p-o)^2))

S <- 50
M <- matrix(NA, S, length(feats), dimnames=list(NULL, names(feats)))
predmat <- array(NA, c(S, nrow(te), length(feats)))
for (s in 1:S) {
  set.seed(s)
  tr$NOISE <- rnorm(nrow(tr)); te$NOISE <- rnorm(nrow(te))
  for (j in seq_along(feats)) {
    f <- feats[[j]]
    set.seed(1000+s)
    m <- randomForest(x=tr[,f,drop=FALSE], y=tr$LAI_ALS, ntree=800)
    p <- as.numeric(predict(m, te[,f,drop=FALSE]))
    predmat[s,,j] <- p; M[s,j] <- rmse(p, te$LAI_ALS)
  }
}
cat("\n=== RMSE over 50 seeds (mean +/- sd) ===\n")
for (j in seq_along(feats))
  cat(sprintf("  %-10s %.3f +/- %.3f\n", names(feats)[j], mean(M[,j]), sd(M[,j])))

dg <- M[,"structure"] - M[,"dopt"]      # >0 means d_opt helps
dn <- M[,"structure"] - M[,"noise"]     # negative control
cat(sprintf("\nDelta RMSE (structure - +d_opt): mean=%+.3f sd=%.3f | %% seeds d_opt helps=%.0f%%\n",
            mean(dg), sd(dg), 100*mean(dg>0)))
cat(sprintf("Delta RMSE (structure - +noise) [control]: mean=%+.3f (should be ~0)\n", mean(dn)))

# bootstrap test set: paired structure vs +d_opt on mean predictions across seeds
ps <- apply(predmat[,,1],2,mean); pd <- apply(predmat[,,2],2,mean); o <- te$LAI_ALS
set.seed(7)
bd <- replicate(2000, { i<-sample(length(o),replace=TRUE); rmse(ps[i],o[i])-rmse(pd[i],o[i]) })
cat(sprintf("\nBootstrap (test plots, 2000x) Delta RMSE struct-dopt: mean=%+.3f  95%%CI=[%+.3f, %+.3f]  P(>0)=%.0f%%\n",
            mean(bd), quantile(bd,.025), quantile(bd,.975), 100*mean(bd>0)))
