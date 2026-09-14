# Train S2->full-LAI correction on the cLHS LiDAR sample (NOT at HOBO plots),
# test on HOBO plots. Isolates the out-of-sample value of (a) d_opt and (b)
# heterogeneity (CHM_std) as predictors, then stratifies HOBO error by Ch2
# heterogeneity class. No MuSICA (LAI recovery — the diagnostic before forcing).
suppressPackageStartupMessages({ library(terra); library(randomForest) })
NCF <- path.expand("~/Documents/NC_Full")
r_s2   <- rast(file.path(NCF, "output/intermediate/sm6/Blois/s2lai_summer_atbd_T_res_10_m.tif"))
r_chm  <- rast(file.path(NCF, "output/intermediate/sm6/Blois/chm_sd_block_20_m.tif"))
r_dopt <- rast(file.path(NCF, "output/intermediate/lai_als_dopt/Blois/LAI_ALS_dopt_common.tif"))

ex <- function(r, x, y) {
  pts <- vect(cbind(x, y), type = "points", crs = "EPSG:32631")
  as.numeric(terra::extract(r, pts)[, 2])
}

# ---- TRAIN: cLHS sample (full LiDAR LAI in 'LAI'); extract S2/CHM/d_opt ----
cl <- as.data.frame(readRDS("out_files/Sensitivity_Analysis/clhs_sample.rds"))
tr <- data.frame(
  LAI_ALS      = cl$LAI, Hmax = cl$Hmax, fCover = cl$fCover, VCI = cl$VCI,
  LAI_S2_ATBD  = ex(r_s2,  cl$x, cl$y),
  CHM_std      = ex(r_chm, cl$x, cl$y),
  LAI_ALS_DOPT = ex(r_dopt, cl$x, cl$y)
)
tr <- tr[complete.cases(tr), ]

# ---- TEST: HOBO plots (df already has LAI_ALS, S2, structure, d_opt) ----
hb <- readRDS("out_files/Chapter3/lai_prep/df_plots_lai.rds")
te <- data.frame(
  LAI_ALS = hb$LAI_ALS, Hmax = hb$Hmax, fCover = hb$fCover, VCI = hb$VCI,
  LAI_S2_ATBD = hb$LAI_S2_ATBD, LAI_ALS_DOPT = hb$LAI_ALS_DOPT,
  CHM_std = ex(r_chm, hb$x, hb$y)
)
te <- te[complete.cases(te), ]
cat(sprintf("train(cLHS)=%d  test(HOBO)=%d | LAI_ALS mean train=%.2f test=%.2f\n",
            nrow(tr), nrow(te), mean(tr$LAI_ALS), mean(te$LAI_ALS)))

feats <- list(
  s2only        = c("LAI_S2_ATBD"),
  structure     = c("LAI_S2_ATBD","Hmax","fCover","VCI"),
  struct_dopt   = c("LAI_S2_ATBD","Hmax","fCover","VCI","LAI_ALS_DOPT"),
  struct_het    = c("LAI_S2_ATBD","Hmax","fCover","VCI","CHM_std"),
  struct_all    = c("LAI_S2_ATBD","Hmax","fCover","VCI","LAI_ALS_DOPT","CHM_std")
)
set.seed(42)
metr <- function(pr, ob) c(rmse=sqrt(mean((pr-ob)^2)), r2=cor(pr,ob)^2, bias=mean(pr-ob))
res <- do.call(rbind, lapply(names(feats), function(fn) {
  f <- feats[[fn]]
  m <- randomForest(x = tr[, f, drop=FALSE], y = tr$LAI_ALS, ntree = 800)
  pr <- as.numeric(predict(m, te[, f, drop=FALSE]))
  data.frame(model = fn, t(round(metr(pr, te$LAI_ALS), 3)))
}))
cat("\n=== cLHS-train -> HOBO-test recovery of full LAI_ALS ===\n")
print(res[order(res$rmse), ], row.names = FALSE)

# ---- Heterogeneity stratification (Ch2): error of best model by CHM_std tercile ----
best <- res$model[which.min(res$rmse)]
f <- feats[[best]]
m <- randomForest(x = tr[, f, drop=FALSE], y = tr$LAI_ALS, ntree = 800)
te$pred <- as.numeric(predict(m, te[, f, drop=FALSE]))
te$abserr <- abs(te$pred - te$LAI_ALS)
q <- quantile(te$CHM_std, c(0, 1/3, 2/3, 1), na.rm = TRUE)
te$het <- cut(te$CHM_std, q, labels = c("Low","Med","High"), include.lowest = TRUE)
cat(sprintf("\n=== HOBO error by Ch2 heterogeneity class (best model: %s) ===\n", best))
strat <- aggregate(cbind(MAE = abserr) ~ het, data = te, mean)
strat$n <- as.numeric(table(te$het)[strat$het])
strat$CHM_std_mean <- aggregate(CHM_std ~ het, data = te, mean)$CHM_std
print(format(strat, digits = 3), row.names = FALSE)
cat(sprintf("\nCHM_std vs |error| correlation: r=%.2f\n",
            cor(te$CHM_std, te$abserr, use = "complete.obs")))
