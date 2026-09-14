# ==============================================================================
# ANNEX — Shapley vs leave-one-out (LOO) vs forward (add-one), on the GLOBAL
# baseline, ΔTmax (all period + hot days), from cached nc (nc_shapley2x_global).
# Demonstrates that under LAI↔fCover collinearity, LOO under-credits both
# correlated traits (masking → ~0) and forward over-credits both (double count),
# while Shapley splits the shared credit and brackets the two. NO MuSICA re-run.
#   Rscript c1_loo_forward_chunk.R <chunk> <K>
# Writes out_files/Chapter1/tables/loofwd_parts/part_<chunk>.csv
#   (long: pid, Cluster, metric, trait, shapley, loo, forward)
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(lubridate); library(data.table); library(dplyr) })
args <- commandArgs(trailingOnly=TRUE); chunk <- as.integer(args[1]); K <- as.integer(args[2])
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
SHIFT <- 2L; Z_FIX <- 1.0; HOTQ <- 0.90; ds <- CFG$date_seq
ncdir <- "out_files/Chapter1/nc_shapley2x_global"

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N, min(.N,100))], by=Cluster]; sub[, pid:=sprintf("S%04d",.I)]
my <- which(((seq_len(nrow(sub))-1) %% K)==(chunk-1))
dm <- as.data.table(extract_macro_daily(CFG$forcing_file, ds))
hot <- dm[Tmax_macro >= quantile(Tmax_macro, HOTQ, na.rm=TRUE), date]
Fv <- c("LAI","Hmax","fCover","LAD"); n <- 4
coal <- as.matrix(expand.grid(LAI=0:1,Hmax=0:1,fCover=0:1,LAD=0:1)); bits_all <- apply(coal,1,paste,collapse="")
wsh <- function(s) factorial(s)*factorial(n-s-1)/factorial(n)

dtmax_one <- function(path) {   # returns c(all, hot) ΔTmax for one coalition nc
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(c(NA,NA))
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu),tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE)
  zl <- rh*vh; Z<-Z_FIX
  if (Z<=zl[1]){ilo<-1L;ihi<-1L;w<-0} else if (Z>=zl[length(zl)]){ilo<-length(zl);ihi<-ilo;w<-0} else {ilo<-max(which(zl<=Z));ihi<-ilo+1L;w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo])}
  tvec <- t0 + dhours(th) - hours(SHIFT)
  Tc <- ((1-w)*Tk[ilo,]+w*Tk[ihi,]) - 273.15
  dd <- data.table(date=as.Date(floor_date(tvec,"hour")), Tc=Tc)[date %in% ds, .(Tmax=max(Tc)), by=date]
  m1 <- merge(dd, dm, by="date")
  c(all = mean(m1$Tmax - m1$Tmax_macro, na.rm=TRUE), hot = mean(m1[date %in% hot, Tmax - Tmax_macro], na.rm=TRUE))
}
shap <- function(vals){ phi<-setNames(numeric(n),Fv)
  for (v in Fv){ others<-setdiff(Fv,v); acc<-0
    for (m in 0:length(others)) for (S in (if(m==0) list(character(0)) else combn(others,m,simplify=FALSE))){
      b0<-setNames(integer(n),Fv); b0[S]<-1L; b1<-b0; b1[v]<-1L
      f1<-vals[paste(b1[Fv],collapse="")]; f0<-vals[paste(b0[Fv],collapse="")]
      if (is.finite(f1)&&is.finite(f0)) acc<-acc+wsh(m)*(f1-f0) }
    phi[v]<-acc }; phi }
sgl <- function(tr){ b<-setNames(integer(n),Fv); b[tr]<-1L; paste(b[Fv],collapse="") }   # add-one coalition
lo  <- function(tr){ b<-setNames(rep(1L,n),Fv); b[tr]<-0L; paste(b[Fv],collapse="") }     # leave-one-out coalition

one_plot <- function(i){
  prow <- as.data.frame(sub[i])
  Va <- setNames(rep(NA_real_,16), bits_all); Vh <- Va
  for (b in bits_all){ d<-dtmax_one(file.path(ncdir, sprintf("%s_%s.nc", prow$pid, b))); Va[b]<-d["all"]; Vh[b]<-d["hot"] }
  mk <- function(V, lab){
    ph <- shap(V)
    rbindlist(lapply(Fv, function(tr) data.table(pid=prow$pid, Cluster=prow$Cluster, metric=lab, trait=tr,
      shapley=ph[tr], forward=V[sgl(tr)]-V["0000"], loo=V["1111"]-V[lo(tr)])))
  }
  rbind(mk(Va,"Tmax_all"), mk(Vh,"Tmax_hot"))
}
pdir <- "out_files/Chapter1/tables/loofwd_parts"; dir.create(pdir,recursive=TRUE,showWarnings=FALSE)
pf <- file.path(pdir, sprintf("part_%02d.csv", chunk)); if (file.exists(pf)) file.remove(pf); nok<-0L
for (j in seq_along(my)){ r<-tryCatch(one_plot(my[j]),error=function(e) NULL)
  if(!is.null(r)){ fwrite(r,pf,append=file.exists(pf)); nok<-nok+1L }
  if(j%%5==0) cat(sprintf("  chunk %d: %d/%d\n",chunk,j,length(my))) }
cat(sprintf("chunk %d DONE (%d/%d)\n", chunk, nok, length(my)))
