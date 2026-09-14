# ==============================================================================
# Robustesse au décalage temporel -2 h (degré de liberté du chercheur).
# Recalcule, sur les 400 sims BASE (nc_shapley2x *_1111), ΔTmax (max journalier)
# pour SHIFT ∈ {0,1,2,3} h + la pente micro~ERA5 (qui N'utilise PAS le shift, l.39
# de c1_metrics6) → confirme l'insensibilité. NO new sims.
#   Rscript c1_shift_sensitivity_check.R
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(lubridate); library(data.table); library(dplyr); source("R/cluster_relabel.R") })
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_|cluster_relabel)", src)]
invisible(lapply(src, source))                                  # CFG, build_era5_hourly, extract_macro_daily
Z_FIX <- 1.0; SHIFTS <- 0:3
ds   <- CFG$date_seq
dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, ds))
era5 <- as.data.table(build_era5_hourly(CFG$forcing_file, ds))

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N,min(.N,100))], by=Cluster]; sub[, pid:=sprintf("S%04d",.I)]
sub[, P := relabel_cluster(Cluster)]

one <- function(path) {   # ΔTmax pour chaque SHIFT + slope (sans shift)
  out <- c(setNames(rep(NA_real_, length(SHIFTS)), paste0("dT_s", SHIFTS)), slope=NA_real_)
  if (!file.exists(path) || file.size(path) < 1e5) return(out)
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(out); on.exit(nc_close(nc))
  if (!all(c("Tair_z","relative_height","veget_height_top") %in% names(nc$var))) return(out)
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu), tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE)
  zl <- rh*vh; Z <- Z_FIX
  if (Z<=zl[1]){ilo<-1L;ihi<-1L;w<-0}else if(Z>=zl[length(zl)]){ilo<-length(zl);ihi<-ilo;w<-0}else{ilo<-max(which(zl<=Z));ihi<-ilo+1L;w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo])}
  Tc <- ((1-w)*Tk[ilo,] + w*Tk[ihi,]) - 273.15
  for (s in SHIFTS) {
    tvec <- t0 + dhours(th) - lubridate::hours(s)
    dd <- data.table(date=as.Date(floor_date(tvec,"hour")), Tc=Tc)[date %in% ds, .(Tmax=max(Tc)), by=date]
    m1 <- merge(dd, dm, by="date")
    out[paste0("dT_s", s)] <- mean(m1$Tmax - m1$Tmax_macro, na.rm=TRUE)
  }
  mic <- data.table(time=floor_date(t0+dhours(th),"hour"), Tmic=Tc); mm <- merge(mic[as.Date(time)%in%ds], era5, by="time")
  if (nrow(mm)>10) out["slope"] <- as.numeric(coef(lm(Tmic~Tair_era5, mm))[2])
  out
}

R <- rbindlist(lapply(seq_len(nrow(sub)), function(i){
  v <- one(file.path("out_files/Chapter1/nc_shapley2x", sprintf("%s_1111.nc", sub$pid[i])))
  as.data.table(c(list(pid=sub$pid[i], P=as.character(sub$P[i])), as.list(v)))
}))
R <- R[is.finite(dT_s2)]
cat(sprintf("n plots = %d\n\n", nrow(R)))

# 1) ΔTmax par archétype, par shift
cat("=== ΔTmax moyen (°C) par archétype × shift ===\n")
tab <- R[, lapply(.SD, function(x) round(mean(x,na.rm=TRUE),3)), by=P, .SDcols=paste0("dT_s",SHIFTS)][order(P)]
print(tab)

# 2) écart max au choix canonique (s2) + corrélation plot-à-plot
cat("\n=== écart à SHIFT=2 (référence) ===\n")
for (s in setdiff(SHIFTS,2)) {
  d <- R[[paste0("dT_s",s)]] - R$dT_s2
  cat(sprintf("  SHIFT=%d : mean Δ=%+.4f °C  max|Δ|=%.4f °C  r(vs s2)=%.4f\n",
              s, mean(d,na.rm=TRUE), max(abs(d),na.rm=TRUE), cor(R[[paste0("dT_s",s)]], R$dT_s2, use="complete.obs")))
}
cat("\n=== pente micro~ERA5 : identique quel que soit le shift (n'utilise pas SHIFT, c1_metrics6 l.39) ===\n")
cat(sprintf("  slope range = [%.4f, %.4f], median = %.4f (un seul calcul, shift-indépendant)\n",
            min(R$slope,na.rm=TRUE), max(R$slope,na.rm=TRUE), median(R$slope,na.rm=TRUE)))
fwrite(R, "out_files/Chapter1/tables/tab_shift_sensitivity_check.csv")
cat("\nDONE -> tab_shift_sensitivity_check.csv\n")
