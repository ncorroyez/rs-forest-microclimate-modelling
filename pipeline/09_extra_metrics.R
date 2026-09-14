# ==============================================================================
# PIPELINE STAGE 09 (ANNEX) — extra response metrics: ΔTmin, diurnal amplitude
# (Tmax−Tmin), and T-stability (SD of hourly Tair). Exact Shapley per plot,
# all-period, at 1 m. Kept OUT of the main spine (ΔTmax/ΔVPDmax). Temperature-
# only metrics → HOBO-validatable in principle (annex shows attribution only).
# Sign of "buffering" differs per metric (see caption); φ is on the raw metric.
# Outputs (OUT_FIG/annex): fig_annex_shapley_extra_by_cluster.{png,pdf}
#                          tab_annex_shapley_extra_ci.csv
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
ANN <- file.path(PIPE$OUT_FIG, "annex"); dir.create(ANN, recursive=TRUE, showWarnings=FALSE)
cli_h1("STAGE 09 (annex) — ΔTmin / amplitude / T-stability")

bits <- PIPE$BITS
hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet=TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
ids <- hobo_pts$id_plot
nc_of <- function(bit, id) file.path(PIPE$MUSICA_DIR_HOBO, bit, sprintf("musica_out_HOBO_%s.nc", id))

# daily Tmin/Tmax @1m (interp, -2h) + hourly Tc for stability
extra_1m <- function(path) {
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(NULL); on.exit(nc_close(nc))
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu), tz="UTC")
  tvec <- t0 + dhours(ncvar_get(nc,"time")) - lubridate::hours(PIPE$TMAX_SHIFT_HR)
  Tk <- ncvar_get(nc,"Tair_z"); rh <- ncvar_get(nc,"relative_height")
  vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE); zl <- rh*vh; Z <- PIPE$Z_FIX
  if (Z<=zl[1]) {ilo<-1L;ihi<-1L;w<-0} else if (Z>=zl[length(zl)]) {ilo<-length(zl);ihi<-ilo;w<-0} else {
    ilo<-max(which(zl<=Z)); ihi<-ilo+1L; w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo]) }
  Tc <- ((1-w)*Tk[ilo,] + w*Tk[ihi,]) - 273.15
  dt <- data.table(date=as.Date(floor_date(tvec,"hour")), Tc=Tc)[date %in% CFG$date_seq]
  dd <- dt[, .(Tmin=min(Tc), Tmax=max(Tc)), by=date]
  data.table(Tmin=mean(dd$Tmin), Amp=mean(dd$Tmax-dd$Tmin), Stab=sd(dt$Tc))
}

cli_alert("Extracting ΔTmin/Amp/Stab @1m, 16 × {length(ids)} plots ...")
rows <- list()
for (bit in bits) for (id in ids) {
  e <- extra_1m(nc_of(bit,id)); if (is.null(e)) next
  rows[[length(rows)+1L]] <- data.table(id_plot=id, bit=bit, Tmin=e$Tmin, Amp=e$Amp, Stab=e$Stab)
}
A <- rbindlist(rows)

# exact Shapley per plot per metric
shap_one <- function(vec){ v0<-as.numeric(vec[PIPE$BASE_BIT]); v<-function(b) as.numeric(vec[b])-v0
  mb<-function(bo,n=4L){b<-rep("0",n);b[bo]<-"1";paste(b,collapse="")}; out<-numeric(4L)
  for (vb in 1:4){others<-setdiff(1:4,vb);tot<-0
    for (s in 0:3){subs<-if(s==0) list(integer(0)) else combn(others,s,simplify=FALSE)
      w<-factorial(s)*factorial(4-s-1)/factorial(4); for (S in subs) tot<-tot+w*(v(mb(c(S,vb)))-v(mb(S)))}
    out[vb]<-tot}; setNames(out,c("LAI","Hmax","fCover","LAD")) }
phi_for <- function(col, mlab){ W<-dcast(A, id_plot~bit, value.var=col)
  rbindlist(lapply(seq_len(nrow(W)), function(i){ vec<-as.numeric(unlist(W[i,..bits])); names(vec)<-bits
    if (any(is.na(vec))) return(NULL); p<-shap_one(vec)
    data.table(id_plot=W$id_plot[i], trait=names(p), phi=unname(p), metric=mlab) })) }
phi <- rbindlist(list(phi_for("Tmin","Delta*T[min]~(degree*C)"),
                      phi_for("Amp","Diurnal~amplitude~(degree*C)"),
                      phi_for("Stab","T-stability~(SD*','~degree*C)")))

# clusters
df_forest <- as.data.table(readRDS(PIPE$CLUSTER_SAMPLE))[, .(x,y,Cluster=as.character(Cluster))]
hxy <- sf::st_coordinates(hobo_pts); clu <- data.table(id_plot=ids, x=hxy[,"X"], y=hxy[,"Y"])
clu[, Cluster := relabel_cluster(sapply(seq_len(.N), function(i)
  df_forest$Cluster[which.min((df_forest$x-x[i])^2+(df_forest$y-y[i])^2)]))]
phi <- merge(phi, clu[,.(id_plot,Cluster)], by="id_plot")
phiC <- rbind(phi, copy(phi)[, Cluster:="All"]); phiC[, Cluster:=factor(Cluster, levels=c(paste0("P",1:4),"All"))]
mlevs <- c("Delta*T[min]~(degree*C)","Diurnal~amplitude~(degree*C)","T-stability~(SD*','~degree*C)")
phiC[, metric := factor(metric, levels=mlevs)]

# bootstrap CI of median + rank (table)
boot_ci <- function(x){x<-x[is.finite(x)];n<-length(x); if(n<2) return(c(median=if(n)median(x) else NA,NA,NA))
  set.seed(42L); m<-vapply(1:2000,function(b) median(x[sample.int(n,n,TRUE)]),numeric(1))
  c(median=median(x), quantile(m,c(.025,.975),names=FALSE))}
ci <- phiC[, as.list(boot_ci(phi)), by=.(metric,Cluster,trait)]
setnames(ci, c("metric","Cluster","trait","median","ci_lo","ci_hi"))
ci[, excludes_zero := is.finite(ci_lo) & (ci_lo>0 | ci_hi<0)]
ci[, rank := frank(-abs(median), ties.method="min"), by=.(metric,Cluster)]
fwrite(ci[order(metric,Cluster,rank)], file.path(PIPE$OUT_TAB, "tab_annex_shapley_extra_ci.csv"))
cli_h2("Annex extra-metric median φ (All)"); print(ci[Cluster=="All"][order(metric,rank), .(metric,rank,trait,median=round(median,3))])

# by-cluster boxplot (3 metric rows × 5 cluster cols), free-y per row
tl <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD"); ordv<-c("LAI","fCover","LAD","Hmax")
phiC[, trait_lab := factor(tl[trait], levels=tl[ordv])]
fill <- c(LAI="#4575B4", Hmax="#91BFDB", fCover="#FC8D59", LAD="#D73027")
yl_of <- function(d){z<-quantile(d,c(.03,.97),na.rm=TRUE);p<-diff(z)*.1;c(z[1]-p,z[2]+p)}
rrng <- phiC[, {yl<-yl_of(phi); .(ymin=yl[1],ymax=yl[2])}, by=metric]
blankC <- rbindlist(lapply(levels(phiC$Cluster), function(cc)
  data.table(metric=rep(rrng$metric,2), phi=c(rrng$ymin,rrng$ymax), Cluster=cc)))
blankC[, `:=`(trait_lab=factor(tl["LAI"],levels=tl[ordv]), Cluster=factor(Cluster,levels=levels(phiC$Cluster)),
              metric=factor(metric,levels=mlevs))]
jitC <- merge(phiC, rrng, by="metric")[phi>=ymin & phi<=ymax]

p <- ggplot(phiC, aes(trait_lab, phi, fill=trait)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_blank(data=blankC, aes(trait_lab, phi), inherit.aes=FALSE) +
  geom_boxplot(width=0.62, outlier.shape=NA, alpha=0.85) +
  geom_jitter(data=jitC, width=0.12, size=0.7, alpha=0.35, colour="grey20") +
  facet_grid(metric ~ Cluster, scales="free_y", labeller=labeller(metric=label_parsed)) +
  scale_x_discrete(labels=function(x) parse(text=x)) + scale_fill_manual(values=fill) +
  labs(x=NULL, y=expression(varphi[v]^"Shapley"),
       title="ANNEX — Shapley attribution of extra metrics (all summer, 1 m)",
       caption="φ>0 raises the metric, φ<0 lowers it. Buffering: ΔTmin → φ>0 (warmer nights); amplitude & T-stability → φ<0 (damped/stable). P1 sparse → P4 dense.") +
  theme_bw(base_size=15) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold", size=12), axis.text.x=element_text(size=12,angle=30,hjust=1),
        plot.caption=element_text(size=10,colour="grey40"), legend.position="none")
ggsave(file.path(ANN,"fig_annex_shapley_extra_by_cluster.png"), p, width=16, height=9, dpi=300, bg="white")
ggsave(file.path(ANN,"fig_annex_shapley_extra_by_cluster.pdf"), p, width=16, height=9, device=cairo_pdf)
cli_alert_success("Saved annex extra-metric figure + CI table")
