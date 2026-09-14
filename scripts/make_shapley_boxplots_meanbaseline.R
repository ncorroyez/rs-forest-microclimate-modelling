# ==============================================================================
# EXACT SHAPLEY φ boxplots on ΔTmax & ΔVPDmax at 1 m, mean-fCover baseline.
# Replaces the LOO boxplots. 16-coalition set assembled with mean fCover baseline:
#   file resolution = "v10_fcovmean if present, else v9"  →  fCover-baseline (f0)
#   coalitions come from v10 (fcov_b = mean = 0.869); fCover-real (f1) from v9
#   (identical to v10). Both extracted at FIXED 1 m (interp) with the legacy
#   -2 h shift, to match extract_deltatmax_one.
#
#   φ_v = trait v's share of metric(REF=1111) − metric(Baseline=0000).
#   φ_v < 0 = trait buffers (lowers daily max). Σφ = REF − Baseline (additivity).
# Periods: whole season & 10% hottest days (macro Tmax).
# Outputs: fig_shapley_boxplot_{tmax,vpd}_by_cluster.{png,pdf}
#          fig_shapley_boxplot_pooled.{png,pdf} + tab_shapley_perplot_meanbaseline.csv
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4); library(sf)
  library(tidyverse); library(lubridate); library(patchwork)
})
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/musica.R")); source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT   <- here::here("outputs/figs_MEB2026_final")
P_HPA <- 1013; Z_FIX <- 1.0
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))
bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
          "1000","1001","1010","1011","1100","1101","1110","1111")

df_macro <- as.data.table(extract_macro_daily(CFG$forcing_file, CFG$date_seq))
thr_hot  <- quantile(df_macro$Tmax_macro, 0.90, na.rm = TRUE)
hot_days <- df_macro[Tmax_macro >= thr_hot, date]
cli_alert("Hot threshold {round(thr_hot,1)}°C ; {length(hot_days)} days")

# ---- file resolution: mean-baseline set --------------------------------------
nc_path <- function(bit, id) {
  v10 <- here::here("out_files/musica_hobo_v10_fcovmean", bit, sprintf("musica_out_HOBO_%s.nc", id))
  v9  <- here::here("out_files/musica_hobo_v9", bit, sprintf("musica_out_HOBO_%s.nc", id))
  if (file.exists(v10)) v10 else v9
}

# daily Tmax & VPDmax at fixed 1 m (interp), legacy -2 h shift
daily_1m <- function(path) {
  nc <- try(nc_open(path), silent = TRUE); if (inherits(nc, "try-error")) return(NULL)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value
  t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  tvec <- t0 + dhours(ncvar_get(nc, "time")) - lubridate::hours(2)   # legacy shift
  Tk <- ncvar_get(nc, "Tair_z"); wmr <- ncvar_get(nc, "wair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm=TRUE)
  zl <- rh * vh
  if (Z_FIX <= zl[1]) { ilo<-1L; ihi<-1L; w<-0 }
  else if (Z_FIX >= zl[length(zl)]) { ilo<-length(zl); ihi<-ilo; w<-0 }
  else { ilo <- max(which(zl <= Z_FIX)); ihi <- ilo+1L; w <- (Z_FIX-zl[ilo])/(zl[ihi]-zl[ilo]) }
  Tc <- ((1-w)*Tk[ilo,] + w*Tk[ihi,]) - 273.15
  wv <- (1-w)*wmr[ilo,] + w*wmr[ihi,]
  vpd <- pmax(esat_hpa(Tc) - (wv/(1+wv))*P_HPA, 0)/10              # kPa
  data.table(date = as.Date(floor_date(tvec,"hour")), Tc = Tc, vpd = vpd)[
    date %in% CFG$date_seq, .(Tmax = max(Tc), VPDmax = max(vpd)), by = date]
}

# ---- collect 16 coalitions × 53 plots ---------------------------------------
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
ids <- hobo_pts$id_plot
cli_h1("Extracting Tmax/VPDmax @1m, 16 coalitions × {length(ids)} plots")
rows <- list()
for (bit in bits) for (id in ids) {
  d <- daily_1m(nc_path(bit, id)); if (is.null(d) || nrow(d)==0) next
  rows[[length(rows)+1L]] <- data.table(
    id_plot=id, bit=bit,
    Tmax_all=mean(d$Tmax), Tmax_hot=mean(d[date %in% hot_days, Tmax]),
    VPD_all=mean(d$VPDmax), VPD_hot=mean(d[date %in% hot_days, VPDmax]))
}
A <- rbindlist(rows)
miss <- A[, .N, by=id_plot][N < 16]
if (nrow(miss)) cli_alert_warning("plots with <16 coalitions: {nrow(miss)}")

# ---- exact Shapley per plot --------------------------------------------------
shap_one <- function(vec) {
  v0 <- as.numeric(vec["0000"]); v <- function(b) as.numeric(vec[b]) - v0
  mb <- function(bo,n=4L){ b<-rep("0",n); b[bo]<-"1"; paste(b,collapse="") }
  out <- numeric(4L)
  for (vb in 1:4){ others<-setdiff(1:4,vb); tot<-0
    for (s in 0:3){ subs <- if(s==0) list(integer(0)) else combn(others,s,simplify=FALSE)
      w <- factorial(s)*factorial(4-s-1)/factorial(4)
      for (S in subs) tot <- tot + w*(v(mb(c(S,vb)))-v(mb(S))) }
    out[vb] <- tot }
  setNames(out, c("LAI","Hmax","fCover","LAD"))
}
phi_for <- function(col, period_lbl) {
  W <- dcast(A, id_plot ~ bit, value.var = col)
  rbindlist(lapply(seq_len(nrow(W)), function(i){
    vec <- as.numeric(unlist(W[i, ..bits])); names(vec) <- bits
    if (any(is.na(vec))) return(NULL)
    p <- shap_one(vec)
    data.table(id_plot=W$id_plot[i], trait=names(p), phi=unname(p), period=period_lbl) }))
}
phi <- rbindlist(list(
  cbind(metric="Delta*T[max]~(degree*C)", phi_for("Tmax_all","All period")),
  cbind(metric="Delta*T[max]~(degree*C)", phi_for("Tmax_hot","10% hottest")),
  cbind(metric="Delta*VPD[max]~(kPa)",    phi_for("VPD_all","All period")),
  cbind(metric="Delta*VPD[max]~(kPa)",    phi_for("VPD_hot","10% hottest"))))

# ---- cluster labels ----------------------------------------------------------
df_forest <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))[, .(x,y,Cluster=as.character(Cluster))]
hxy <- sf::st_coordinates(hobo_pts)
clu <- data.table(id_plot=ids, x=hxy[,"X"], y=hxy[,"Y"])
clu[, Cluster := relabel_cluster(sapply(seq_len(.N), function(i)
  df_forest$Cluster[which.min((df_forest$x-x[i])^2+(df_forest$y-y[i])^2)]))]
phi <- merge(phi, clu[,.(id_plot,Cluster)], by="id_plot")
phi[, period := factor(period, levels=c("All period","10% hottest"))]
phi[, Cluster := factor(Cluster, levels=paste0("P",1:4))]
fwrite(phi, file.path(OUT, "tab_shapley_perplot_meanbaseline.csv"))

cli_h2("Median Shapley φ per metric × period × trait")
print(phi[, .(median=round(median(phi,na.rm=TRUE),3)), by=.(metric,period,trait)][order(metric,period,median)])

# additivity check
chk <- phi[metric=="Delta*T[max]~(degree*C)" & period=="All period", .(sphi=sum(phi)), by=id_plot]
W <- dcast(A, id_plot~bit, value.var="Tmax_all")[, .(id_plot, tot=`1111`-`0000`)]
cli_alert("Additivity ΔTmax: median Σφ={round(median(merge(chk,W,by='id_plot')$sphi),3)} vs REF-Base={round(median(W$tot),3)}")

# ---- plotting helpers --------------------------------------------------------
trait_lab <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD")
ordv <- c("LAI","fCover","LAD","Hmax")
phi[, trait_lab := factor(trait_lab[trait], levels=trait_lab[ordv])]
fill <- c(LAI="#4575B4", Hmax="#91BFDB", fCover="#FC8D59", LAD="#D73027")
yl_of <- function(d,q=c(.04,.96)){ z<-quantile(d,q,na.rm=TRUE); p<-diff(z)*.10; c(z[1]-p,z[2]+p) }

by_cluster <- function(mk, ylab, fname){
  d <- phi[metric==mk]; yl <- yl_of(d$phi); nc <- sum(d$phi<yl[1]|d$phi>yl[2],na.rm=TRUE)
  p <- ggplot(d, aes(trait_lab, phi, fill=trait)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    geom_boxplot(width=0.62, outlier.shape=NA, alpha=0.85) +
    geom_jitter(width=0.12, size=0.9, alpha=0.4, colour="grey20") +
    facet_grid(period ~ Cluster) + coord_cartesian(ylim=yl) +
    scale_x_discrete(labels=function(x) parse(text=x)) + scale_fill_manual(values=fill) +
    labs(x=NULL, y=ylab,
         caption=sprintf("Exact Shapley, 1 m, fCover baseline=mean. φ<0 = buffers. P1 sparse → P4 dense.%s",
                         if(nc>0) sprintf(" %d clipped.",nc) else "")) +
    theme_bw(base_size=18) +
    theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
          strip.text=element_text(face="bold", size=16), axis.text.x=element_text(size=14,angle=30,hjust=1),
          plot.caption=element_text(size=11,colour="grey40"), legend.position="none")
  ggsave(file.path(OUT,paste0(fname,".png")), p, width=16, height=8.5, dpi=300, bg="white")
  ggsave(file.path(OUT,paste0(fname,".pdf")), p, width=16, height=8.5, device=cairo_pdf)
  cli_alert_success("Saved {fname}")
}
by_cluster("Delta*T[max]~(degree*C)", expression(varphi[v]^"Shapley"~"on "*Delta*T[max]~"("*degree*"C)"), "fig_shapley_boxplot_tmax_by_cluster")
by_cluster("Delta*VPD[max]~(kPa)",    expression(varphi[v]^"Shapley"~"on "*Delta*VPD[max]~"(kPa)"),       "fig_shapley_boxplot_vpd_by_cluster")

# pooled (metric × period)
pooled_panel <- function(mk, ylab){
  d <- phi[metric==mk]; yl <- yl_of(d$phi, c(.02,.98)); nc <- sum(d$phi<yl[1]|d$phi>yl[2],na.rm=TRUE)
  ggplot(d, aes(trait_lab, phi, fill=trait)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.85) +
    geom_jitter(width=0.12, size=1.0, alpha=0.35, colour="grey20") +
    facet_wrap(~period) + coord_cartesian(ylim=yl) +
    scale_x_discrete(labels=function(x) parse(text=x)) + scale_fill_manual(values=fill) +
    labs(x=NULL, y=ylab, subtitle=if(nc>0) sprintf("(%d points clipped)",nc) else NULL) +
    theme_bw(base_size=17) +
    theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
          strip.text=element_text(face="bold"), axis.text.x=element_text(face="bold",size=15),
          plot.subtitle=element_text(size=10,colour="grey45"), legend.position="none")
}
p_pool <- (pooled_panel("Delta*T[max]~(degree*C)", expression(varphi[v]^"Shapley"~"on "*Delta*T[max]~"("*degree*"C)")) /
           pooled_panel("Delta*VPD[max]~(kPa)",    expression(varphi[v]^"Shapley"~"on "*Delta*VPD[max]~"(kPa)"))) +
  plot_annotation(caption="Exact Shapley, 1 m, fCover baseline=mean. φ<0 = trait buffers. Σφ = REF − baseline.",
                  theme=theme(plot.caption=element_text(size=11,colour="grey40")))
ggsave(file.path(OUT,"fig_shapley_boxplot_pooled.png"), p_pool, width=13, height=10, dpi=300, bg="white")
ggsave(file.path(OUT,"fig_shapley_boxplot_pooled.pdf"), p_pool, width=13, height=10, device=cairo_pdf)
cli_alert_success("Saved fig_shapley_boxplot_pooled")
