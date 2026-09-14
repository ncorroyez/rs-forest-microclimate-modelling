# ==============================================================================
# Chapter 3 — two paired bootstraps that fix the wording of two load-bearing
# sentences in §3.4:
#  (1) open-canopy parity: is S2 ATBD |ΔTmax error| different from LiDAR full's
#      over the open-canopy plots? Paired per-plot, 95% CI on the mean difference.
#  (2) the LAI threshold where S2 saturation overtakes LiDAR: bootstrap CI on the
#      abs-error LOESS crossing, plus sensitivity to the LOESS span.
# Reads Table4n_perplot_error_structure.csv. Run from z_Example root.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table) })
set.seed(42)
pp <- fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4n_perplot_error_structure.csv")
w <- dcast(pp, id_plot + LAI_ALS ~ scenario, value.var="aerr")
setnames(w, c("LiDAR full","S2 ATBD (dyn)"), c("aL","aS"))
med <- median(unique(pp[,.(id_plot,LAI_ALS)])$LAI_ALS)
cat(sprintf("median LiDAR LAI (open/dense split) = %.2f\n", med))

boot_mean <- function(x, B=5000){ n<-length(x); m<-replicate(B, mean(x[sample.int(n,n,TRUE)])); quantile(m,c(.025,.5,.975)) }
report <- function(lbl, d){ ci<-boot_mean(d)
  cat(sprintf("%-26s n=%2d  mean Δ|err|(S2−LiDAR)=%+.3f °C  95%% CI [%+.3f, %+.3f]  -> %s\n",
              lbl, length(d), mean(d), ci[1], ci[3],
              ifelse(ci[1]<=0 & ci[3]>=0, "INDISTINGUISHABLE (CI spans 0)",
                     ifelse(ci[1]>0,"S2 worse","S2 better")))) }
cat("\n=== (1) paired per-plot |ΔTmax error|: S2 ATBD minus LiDAR full ===\n")
report("open canopy (LAI<med)", w[LAI_ALS<med, aS-aL])
report("dense canopy (LAI>=med)", w[LAI_ALS>=med, aS-aL])
report("all 53 plots", w[, aS-aL])

# (2) bootstrap CI on the abs-error LOESS crossing LAI
g <- pp[scenario %in% c("LiDAR full","S2 ATBD (dyn)"), .(id_plot,LAI_ALS,scenario,aerr)]
gx <- seq(min(g$LAI_ALS), max(g$LAI_ALS), length.out=300)
crossing <- function(d, span){
  fL<-tryCatch(predict(loess(aerr~LAI_ALS, d[scenario=="LiDAR full"], span=span), gx), error=function(e)rep(NA,length(gx)))
  fS<-tryCatch(predict(loess(aerr~LAI_ALS, d[scenario=="S2 ATBD (dyn)"], span=span), gx), error=function(e)rep(NA,length(gx)))
  diff<-fS-fL; ok<-is.finite(diff)
  # last LAI where S2 goes from <= LiDAR to > LiDAR (saturation onset)
  idx<-which(ok & diff[ -length(diff)]<=0 & c(diff[-1],NA)>0)
  if(!length(idx)) return(NA_real_); gx[max(idx)] }
cat("\n=== (2) LAI crossing where S2 |error| overtakes LiDAR ===\n")
for(sp in c(0.75,1.0,1.25)) cat(sprintf("  full-data crossing, LOESS span=%.2f : LAI = %.2f\n", sp, crossing(g, sp)))
ids <- unique(g$id_plot)
bc <- replicate(2000, { s<-sample(ids, length(ids), TRUE)
  d<-rbindlist(lapply(seq_along(s), function(i) g[id_plot==s[i]][, id_plot:=paste0(id_plot,"_",i)]))
  crossing(d, 1.0) })
bc<-bc[is.finite(bc)]; q<-quantile(bc, c(.025,.5,.975))
cat(sprintf("  bootstrap crossing (span=1, %d/%d valid): median LAI=%.2f, 95%% CI [%.2f, %.2f]\n",
            length(bc), 2000, q[2], q[1], q[3]))
cat("\nDONE\n")
