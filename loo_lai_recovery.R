# Leave-one-plot-out recovery of full LiDAR LAI (LAI_ALS) from various predictors.
# Tests whether d_opt physics / S2-only methods can recover full canopy LAI out
# of sample (proxy for "beyond LiDAR"), vs structure-based RF (needs LiDAR feats).
# No MuSICA — pure LAI cross-validation. SAFE sourcing not needed.
suppressPackageStartupMessages(library(randomForest))
p <- readRDS("out_files/Chapter3/lai_prep/df_plots_lai.rds")
d_opt <- 7
need <- c("LAI_ALS","LAI_S2_ATBD","LAI_ALS_DOPT","Hmax","fCover","VCI","LCV")
p <- p[complete.cases(p[, need]), need]
n <- nrow(p); cat(sprintf("LOO over %d complete plots\n", n))

# predictor builders: each returns predicted LAI_ALS for the held-out row,
# trained ONLY on the train rows.
methods <- list(
  raw_s2          = function(tr, te) te$LAI_S2_ATBD,
  global_rescale  = function(tr, te) te$LAI_S2_ATBD * mean(tr$LAI_ALS / tr$LAI_S2_ATBD),
  dopt_physics    = function(tr, te) {
    # S2 sees top d_opt of canopy height -> amplify by H/d_opt, calibrated on train
    raw <- te$LAI_S2_ATBD * (te$Hmax / d_opt)
    cal <- mean(tr$LAI_ALS) / mean(tr$LAI_S2_ATBD * (tr$Hmax / d_opt))
    raw * cal
  },
  dopt_lai_rescale= function(tr, te) {
    # uses per-plot LAI_ALS_DOPT (needs LiDAR) — only valid where LiDAR exists
    te$LAI_S2_ATBD * mean(tr$LAI_ALS / tr$LAI_ALS_DOPT)
  },
  rf_s2only       = function(tr, te) {
    m <- randomForest(LAI_ALS ~ LAI_S2_ATBD, data = tr, ntree = 500)
    as.numeric(predict(m, te))
  },
  rf_struct_only  = function(tr, te) {
    m <- randomForest(LAI_ALS ~ Hmax + fCover + VCI + LCV, data = tr, ntree = 500)
    as.numeric(predict(m, te))
  },
  rf_full         = function(tr, te) {
    m <- randomForest(LAI_ALS ~ LAI_S2_ATBD + Hmax + fCover + VCI + LCV, data = tr, ntree = 500)
    as.numeric(predict(m, te))
  }
)

set.seed(42)
preds <- sapply(names(methods), function(mn) {
  sapply(seq_len(n), function(i) methods[[mn]](p[-i, ], p[i, ]))
})
obs <- p$LAI_ALS
res <- do.call(rbind, lapply(colnames(preds), function(mn) {
  pr <- preds[, mn]
  data.frame(method = mn,
             rmse = sqrt(mean((pr - obs)^2)),
             r2   = cor(pr, obs)^2,
             bias = mean(pr - obs))
}))
res <- res[order(res$rmse), ]
cat("\n=== LOO recovery of full LAI_ALS (lower RMSE = better) ===\n")
print(format(res, digits = 3), row.names = FALSE)
cat("\nNote: methods using only LAI_S2_ATBD (raw_s2, global_rescale, rf_s2only)\n")
cat("are the only ones applicable BEYOND LiDAR. dopt_physics adds the d_opt depth.\n")
cat("rf_struct_only / rf_full / dopt_lai_rescale require LiDAR-derived inputs.\n")

# --- persist for the chapter (durable source for plan §3.2) ---
res$applicable_beyond_lidar <- res$method %in% c("raw_s2","global_rescale","rf_s2only")
write.csv(res, "/home/corroyez/Documents/NC_Full/output/tables/Table6_LAIrecovery_LOO.csv", row.names = FALSE)
write.csv(res, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table6_LAIrecovery_LOO.csv", row.names = FALSE)
cat("\nWritten Table6_LAIrecovery_LOO.csv\n")
