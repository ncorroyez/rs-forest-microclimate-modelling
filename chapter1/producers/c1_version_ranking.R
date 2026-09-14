# ==============================================================================
# Trait-ranking robustness to MuSICA version (legacy v3.2.0 vs v3.2.3), with NO
# new sims: both versions already have the full 2^4 coalition lattice at the 53
# HOBO plots (out_files/musica_hobo_z05[/_v323]/<coalition>/musica_out_HOBO_*.nc).
# For each version we read ΔTmax (fixed 1 m, -2 h, vs macro) for all 16 coalitions
# x 53 plots, and form each trait's MEAN MARGINAL EFFECT (average over the 8 pairs
# that toggle that trait real-vs-baseline). We then compare the trait ranking and
# the per-plot, per-trait effects between versions. Bit order: LAI,Hmax,fCover,LAD.
#   Rscript c1_version_ranking.R
# Out: tab_version_ranking.csv (+ console)
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
})
CFG_forcing <- "in_files/musica_in_Blois.nc"
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
dm <- as.data.table(extract_macro_daily(CFG_forcing, ds))
Z_FIX <- 1.0; SHIFT <- 2L
dTmax <- function(path) {
  if (!file.exists(path) || file.size(path) < 1000) return(NA_real_)
  nc <- try(nc_open(path), silent = TRUE); if (inherits(nc, "try-error")) return(NA_real_)
  on.exit(nc_close(nc))
  if (!all(c("Tair_z","relative_height","veget_height_top") %in% names(nc$var))) return(NA_real_)
  tu <- ncatt_get(nc, "time", "units")$value; t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time"); Tk <- ncvar_get(nc, "Tair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh; Z <- Z_FIX
  if (Z <= zl[1]) { ilo <- 1L; ihi <- 1L; w <- 0 }
  else if (Z >= zl[length(zl)]) { ilo <- length(zl); ihi <- ilo; w <- 0 }
  else { ilo <- max(which(zl <= Z)); ihi <- ilo + 1L; w <- (Z - zl[ilo]) / (zl[ihi] - zl[ilo]) }
  tvec <- t0 + dhours(th) - lubridate::hours(SHIFT); Tc <- ((1 - w) * Tk[ilo, ] + w * Tk[ihi, ]) - 273.15
  dd <- data.table(date = as.Date(floor_date(tvec, "hour")), Tc = Tc)[date %in% ds, .(Tmax = max(Tc)), by = date]
  m1 <- merge(dd, dm, by = "date"); mean(m1$Tmax - m1$Tmax_macro, na.rm = TRUE)
}

coalitions <- do.call(paste0, expand.grid(0:1,0:1,0:1,0:1)[,4:1])  # "0000".."1111"
traits <- c(LAI=1L, Hmax=2L, fCover=3L, LAD=4L)
VDIR <- c("v3.2.0"="out_files/musica_hobo_z05", "v3.2.3"="out_files/musica_hobo_z05_v323")

# ids present in both versions' REF
ids <- sub("musica_out_HOBO_(.+)\\.nc$","\\1",
           list.files(file.path(VDIR[1],"1111"), "\\.nc$"))
cat(sprintf("plots: %d | coalitions: %d | versions: %d\n", length(ids), length(coalitions), length(VDIR)))

read_version <- function(vdir) {
  M <- matrix(NA_real_, nrow=length(ids), ncol=length(coalitions),
              dimnames=list(ids, coalitions))
  for (co in coalitions) for (id in ids)
    M[id, co] <- dTmax(file.path(vdir, co, sprintf("musica_out_HOBO_%s.nc", id)))
  M
}
marg <- function(M, p) {  # mean marginal effect of trait at bit position p, per plot
  cs <- strsplit(coalitions, "")
  off <- coalitions[sapply(cs, function(b) b[p]=="0")]
  sapply(seq_len(nrow(M)), function(i) {
    on <- sapply(off, function(c0){ b<-strsplit(c0,"")[[1]]; b[p]<-"1"; paste(b,collapse="") })
    mean(M[i, on] - M[i, off], na.rm=TRUE)
  })
}

res <- list(); perplot <- list()
for (v in names(VDIR)) {
  M <- read_version(VDIR[[v]])
  pe <- sapply(names(traits), function(t) marg(M, traits[[t]]))   # plots x traits
  rownames(pe) <- ids; perplot[[v]] <- pe
  imp <- colMeans(abs(pe), na.rm=TRUE)
  res[[v]] <- data.table(version=v, trait=names(imp), importance=round(imp,3),
                         mean_signed=round(colMeans(pe, na.rm=TRUE),3))
  res[[v]] <- res[[v]][order(-importance)][, rank := .I]
  cat(sprintf("\n=== %s : trait ranking (mean |marginal ΔTmax|, 53 HOBO) ===\n", v)); print(res[[v]])
}

OUT <- rbindlist(res); fwrite(OUT, "out_files/Chapter1/tables/tab_version_ranking.csv")
cat("\n=== cross-version agreement (per-plot marginal effect, Spearman) ===\n")
for (t in names(traits)) {
  a <- perplot[["v3.2.0"]][,t]; b <- perplot[["v3.2.3"]][,t]
  cat(sprintf("  %-7s rho = %+.2f  (v320 imp %.3f, v323 imp %.3f)\n", t,
              cor(a,b,method="spearman",use="complete.obs"),
              mean(abs(a),na.rm=TRUE), mean(abs(b),na.rm=TRUE)))
}
o1 <- res[["v3.2.0"]][order(-importance),trait]; o2 <- res[["v3.2.3"]][order(-importance),trait]
cat(sprintf("\nRANK ORDER  v3.2.0: %s\n            v3.2.3: %s\n            IDENTICAL: %s\n",
            paste(o1,collapse=" > "), paste(o2,collapse=" > "), identical(o1,o2)))
cat("DONE -> tab_version_ranking.csv\n")
