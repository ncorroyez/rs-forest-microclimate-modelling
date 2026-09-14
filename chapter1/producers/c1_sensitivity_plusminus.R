# ==============================================================================
# Marginal ΔTmax of ADDING (+1) vs REMOVING (-1) one native unit of each variable,
# at each cluster's mean operating point (supervisor: "une unité en + ou en -
# également ; Hmax 1 m c'est rien"). Reveals the ASYMMETRY that saturation implies:
# near the plateau (dense clusters) adding leaf area does little but removing it
# still warms. Steps: LAI ±1 one-sided unit (±2 stored); Hmax ±5 m (1 m is
# negligible); fCover ±0.1 (capped at [0.5,1], per-0.1 normalized by actual step).
# Legacy binary, cluster-mean LAD shape, 1 m metric. 24 new runs (+ cached baseline).
#   Rscript c1_sensitivity_plusminus.R
# Out: out_files/Chapter1/figures/FigSh_sensitivity_plusminus.png + tab_sensitivity_plusminus.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(ggplot2)
  library(rmusica); library(musica.tools); source("R/cluster_relabel.R")
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
NCDIR <- "out_files/Chapter1/nc_sensitivity"; dir.create(NCDIR, recursive=TRUE, showWarnings=FALSE)
ds <- CFG$date_seq
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
cl_mean <- samp[, .(LAI=mean(LAI), Hmax=mean(Hmax), fCover=mean(fCover)), by=Cluster][order(Cluster)]
LAD_CLUSTER <- make_lad_cluster_type_factory(as.data.frame(samp))

esat_hpa <- function(Tc) 6.108*exp(17.27*Tc/(Tc+237.3)); Z_FIX<-1.0; SHIFT<-2L; P_HPA<-1013
dm <- as.data.table(extract_macro_daily(CFG$forcing_file, ds))
metrics_one <- function(path) {
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(NA_real_); on.exit(nc_close(nc))
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu), tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE)
  zl <- rh*vh; Z <- Z_FIX
  if (Z<=zl[1]) {ilo<-1L;ihi<-1L;w<-0} else if (Z>=zl[length(zl)]) {ilo<-length(zl);ihi<-ilo;w<-0} else {
    ilo<-max(which(zl<=Z));ihi<-ilo+1L;w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo])}
  tvec <- t0+dhours(th)-lubridate::hours(SHIFT); Tc <- ((1-w)*Tk[ilo,]+w*Tk[ihi,])-273.15
  dd <- data.table(date=as.Date(floor_date(tvec,"hour")), Tc=Tc)[date %in% ds, .(Tmax=max(Tc)), by=date]
  m1 <- merge(dd, dm, by="date"); mean(m1$Tmax-m1$Tmax_macro, na.rm=TRUE)
}
run_get <- function(cid, tag, lai, hmax, fcov) {
  nc <- file.path(NCDIR, sprintf("C%d_%s.nc", cid, tag)); prow <- data.frame(Cluster=cid, x=0, y=0)
  sc <- list(lai_fn=function(p) lai, hmax_fn=function(p) hmax, fcover_fn=function(p) fcov,
             lad_fn=LAD_CLUSTER, phenology_fn=NULL)
  run_musica_one(prow, sc, nc, CFG$forcing_file, CFG$musica_cmd)
  if (!file.exists(nc)) return(NA_real_); metrics_one(nc)
}
DLAI_STORED <- 2   # +1 one-sided LAI unit = +2 stored ; DHMAX <- 5 m ; DFC <- 0.1
res <- rbindlist(lapply(seq_len(nrow(cl_mean)), function(i) {
  cid <- cl_mean$Cluster[i]; b <- cl_mean[i]
  base <- run_get(cid, "baseline", b$LAI, b$Hmax, b$fCover)
  # endpoints (clamp fCover to [0.5,1])
  fcp <- min(b$fCover+0.1, 1); fcm <- max(b$fCover-0.1, 0.5)
  laip <- run_get(cid,"LAIp1", b$LAI+DLAI_STORED, b$Hmax, b$fCover); laim <- run_get(cid,"LAIm1", b$LAI-DLAI_STORED, b$Hmax, b$fCover)
  hmp  <- run_get(cid,"Hmaxp5", b$LAI, b$Hmax+5, b$fCover);          hmm  <- run_get(cid,"Hmaxm5", b$LAI, b$Hmax-5, b$fCover)
  fcpr <- run_get(cid,"fCovp", b$LAI, b$Hmax, fcp);                  fcmr <- run_get(cid,"fCovm", b$LAI, b$Hmax, fcm)
  data.table(Cluster=cid,
    # per-unit effects: ADD (+) and REMOVE (-) one unit. fCover normalized to per-0.1 by actual step.
    rbind(
      data.table(variable="LAI~(per~1~unit)",      dir="+1 unit", dT=laip-base),
      data.table(variable="LAI~(per~1~unit)",      dir="-1 unit", dT=laim-base),
      data.table(variable="Hmax~(per~5~m)",        dir="+5 m",    dT=hmp-base),
      data.table(variable="Hmax~(per~5~m)",        dir="-5 m",    dT=hmm-base),
      data.table(variable="fCover~(per~0.1)",      dir="+0.1",    dT=(fcpr-base)/((fcp-b$fCover)/0.1)),
      data.table(variable="fCover~(per~0.1)",      dir="-0.1",    dT=(fcmr-base)/((b$fCover-fcm)/0.1))
    ))
}))
res[, P := relabel_cluster(Cluster)]
fwrite(res, "out_files/Chapter1/tables/tab_sensitivity_plusminus.csv")
cat("=== ΔTmax of +unit vs -unit, per cluster (°C) ===\n")
print(dcast(res, P + variable ~ dir, value.var="dT")[, lapply(.SD, function(x) if (is.numeric(x)) round(x,3) else x)])

res[, variable := factor(variable, levels=c("LAI~(per~1~unit)","fCover~(per~0.1)","Hmax~(per~5~m)"))]
res[, sign := ifelse(grepl("^\\+", dir), "add (+)", "remove (-)")]
p <- ggplot(res, aes(P, dT, fill=sign)) +
  geom_hline(yintercept=0, colour="grey55", linewidth=0.4) +
  geom_col(position=position_dodge(width=0.78), width=0.7, colour="grey25", linewidth=0.2) +
  facet_wrap(~variable, scales="free_y", nrow=1, labeller=label_parsed) +
  scale_fill_manual(values=c("add (+)"="#2C7BB6","remove (-)"="#D7191C"), name=NULL) +
  labs(x=NULL, y=expression(Delta*T[max]~change~("°C")),
       title="Effect of ADDING vs REMOVING one unit, at each cluster operating point") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey92"), strip.text=element_text(face="bold"),
        legend.position="bottom", plot.title=element_text(size=11, face="bold"))
ggsave("out_files/Chapter1/figures/FigSh_sensitivity_plusminus.png", p, width=11, height=4.6, dpi=200, bg="white")
cat("DONE -> FigSh_sensitivity_plusminus.png + tab_sensitivity_plusminus.csv\n")
