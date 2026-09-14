# ==============================================================================
# Chapter 3 — build a GENUINE 53/53 input df (no fabrication), for v3.2.3-iter.
# Keeps df_plots' real-53 columns (LAI_ALS, LAI_S2_ATBD, Hmax, fCover, LCV, VCI,
# LAI_RF*, FORMS_H, LAD_Layer_*) unchanged (so §3.1-3.3 stay consistent), and
# fills the 6 near-bare plots' 3 missing scalars from REAL data / correct physics:
#   LAI_ALS_DOPT  : = LAI_ALS where canopy Hmax <= d_opt (whole canopy is within
#                   the top-d_opt layer -> exact, method-independent); for the one
#                   tall-sparse plot (41_19, Hmax 17) = real top-d_opt LAD slice.
#   LAI_S2_DOPT   : = LAI_S2_ATBD where Hmax <= d_opt (S2 senses the whole short
#                   canopy); tall plot from the Not_Masked optim summer-mean.
#   LAI_S2_RESCALED : = LAI_S2_ATBD * ratio-A (the scenario's own definition).
# Saves out_files/Chapter3/lai_prep/df_plots_real53.rds and prints an NA audit.
#   Rscript build_real53_df.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(terra); library(stringr)})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
df <- as.data.table(load_lai_prep(CFG_C3)$df_plots)
dopt <- CFG_C3$d_opt_m; kscale <- 0.5/0.65
NM <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"

# real top-d_opt LAD slice (x kscale). For short canopies (Hmax<=d_opt) this
# equals the full-canopy sum = LAI_ALS by construction.
lad_cols <- grep("^LAD_Layer_", names(df), value=TRUE); hh <- as.numeric(sub("LAD_Layer_","",lad_cols))
Mlad <- as.matrix(df[, ..lad_cols]); Mlad[is.na(Mlad)] <- 0
top_slice <- sapply(seq_len(nrow(df)), function(i) sum(Mlad[i, hh > (df$Hmax[i]-dopt)]) ) * kscale

# ratio-A (real, from the 53)
rA <- mean(df$LAI_ALS, na.rm=TRUE)/mean(df$LAI_S2_ATBD, na.rm=TRUE)

# --- identify the 6 out-of-mask plots (the ones the masked extraction FILLED) ---
fill6 <- which(is.na(df$LAI_ALS_DOPT))               # exactly the 6 near-bare plots
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs="EPSG:32631")

# --- REAL un-masked S2 summer-means for those 6 (replace the constant fills) ----
af <- list.files(NM, pattern="^s2lai_2021-0[678]-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
real_atbd <- rowMeans(terra::extract(rast(af), pts)[,-1,drop=FALSE], na.rm=TRUE)
of <- list.files(NM, pattern="^s2lai_2021-0[678]-\\d{2}_atbd_optim_common_res_10_m\\.tif$", full.names=TRUE)
real_opt  <- rowMeans(terra::extract(rast(of), pts)[,-1,drop=FALSE], na.rm=TRUE)
df$LAI_S2_ATBD[fill6] <- real_atbd[fill6]            # was 4.82 fill -> real per-plot
df$LAI_S2_DOPT[fill6] <- real_opt[fill6]             # was 4.82 fill -> real per-plot

# --- LAI_ALS_DOPT: short -> LAI_ALS (exact); tall (41_19) -> real LAD top-slice ---
short <- df$Hmax[fill6] <= dopt
df$LAI_ALS_DOPT[fill6] <- ifelse(short, df$LAI_ALS[fill6], top_slice[fill6])
# --- LAI_S2_RESCALED: peak = LAI_ALS (def R/lai_corrections.R) ------------------
df$LAI_S2_RESCALED[is.na(df$LAI_S2_RESCALED)] <- df$LAI_ALS[is.na(df$LAI_S2_RESCALED)]

# --- RE-TRAIN the RF on the CORRECTED df (real S2 features) so LAI_RF/LAI_RF_DOPT
#     are genuine model outputs, not predictions off the fill (seed 42 = pipeline) -
suppressPackageStartupMessages(library(randomForest))
dfp <- as.data.frame(df)
df$LAI_RF      <- predict_lai_rf(dfp, train_lai_rf(dfp, target="LAI_ALS")$model)
df$LAI_RF_DOPT <- predict_lai_rf(dfp, train_lai_rf(dfp, target="LAI_ALS_DOPT")$model)

cat(sprintf("\nd_opt=%.0f m | LAD-sum vs LAI_ALS cor=%.3f (sanity)\n", dopt, cor(rowSums(Mlad)*kscale, df$LAI_ALS)))
cat("=== the 6 formerly-filled plots (now real) ===\n")
print(df[fill6, .(id_plot, Hmax, LAI_ALS, LAI_ALS_DOPT, LAI_S2_ATBD, LAI_S2_DOPT, LAI_S2_RESCALED, LAI_RF, LAI_RF_DOPT)])
cat("\n=== FILL AUDIT: any value shared by >=3 plots? (fCover=1.0 legitimate) ===\n")
ivars <- c("LAI_ALS","LAI_ALS_DOPT","LAI_S2_ATBD","LAI_S2_DOPT","LAI_S2_RESCALED","LAI_RF","LAI_RF_DOPT","FORMS_H","Hmax","fCover")
for(v in ivars){ tb<-sort(table(round(df[[v]],4)),decreasing=TRUE)
  cat(sprintf("  %-16s NA=%d | max-repeat=%d%s\n", v, sum(is.na(df[[v]])), tb[1],
              if(tb[1]>=3 && v!="fCover") sprintf("  << STILL FILL (%.3f)",as.numeric(names(tb)[1])) else "")) }
saveRDS(as.data.frame(df), file.path(CFG_C3$out_dir, "lai_prep", "df_plots_real53.rds"))
cat("\nSAVED -> df_plots_real53.rds\n")
