# ==============================================================================
# Step-by-step forward validation on the OTHER field metric: the micro-macro
# SLOPE (buffering/amplifying), the second HOBO-validatable quantity. Same design
# as the ΔTmax forward (c1_forward_percluster.R): per cluster, add LiDAR variables
# cumulatively in that cluster's own ΔTmax-importance order; at each step compare
# the simulated per-plot slope (1 m, hourly regression vs ERA5, no shift) with the
# observed slope across the cluster's loggers (r, RMSE). Reads the cached 848 z05
# HOBO nc (no MuSICA re-run); slopes cached to CSV for fast re-runs.
#   PIPE_BRANCH=z05 Rscript c1_forward_slope.R
# Out: out_files/Chapter1/figures/FigSh_forward_slope.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork); library(ncdf4); library(lubridate) })
if (!exists("PIPE")) { Sys.setenv(PIPE_BRANCH = "z05"); source("pipeline/00_config.R") }

# 1 m interpolated hourly slope vs ERA5 (no shift) — same interpolation as
# c1_metrics_chunk.R (works on the z05 nc, unlike get_tair_at_z).
slope_at_1m <- function(f, era5) {
  nc <- try(nc_open(f), silent = TRUE); if (inherits(nc, "try-error")) return(NA_real_); on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value; t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time"); Tk <- ncvar_get(nc, "Tair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh; Z <- PIPE$Z_FIX
  if (Z <= zl[1]) { ilo<-1L; ihi<-1L; w<-0 } else if (Z >= zl[length(zl)]) { ilo<-length(zl); ihi<-ilo; w<-0 } else {
    ilo <- max(which(zl <= Z)); ihi <- ilo + 1L; w <- (Z - zl[ilo]) / (zl[ihi] - zl[ilo]) }
  Tc <- ((1 - w) * Tk[ilo, ] + w * Tk[ihi, ]) - 273.15
  mic <- data.table(time = floor_date(t0 + dhours(th), "hour"), Tmic = Tc)
  mm <- merge(mic[as.Date(time) %in% CFG$date_seq], era5, by = "time")
  if (nrow(mm) < 24) return(NA_real_)
  as.numeric(coef(lm(Tmic ~ Tair_era5, mm))[2])
}

Fv <- c("LAI","Hmax","fCover","LAD"); pos <- c(LAI=1L, Hmax=2L, fCover=3L, LAD=4L)
clu <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "clusters.rds")))     # id_plot, Cluster (P1..P4)
V   <- readRDS(file.path(PIPE$OUT_DATA, "ref_validation.rds"))
obs <- as.data.table(V$obs_slope)[, .(id_plot, slope_obs)]

# ---- simulated slope per (plot, bit), cached --------------------------------
slope_csv <- file.path(PIPE$OUT_TAB, "tab_hobo_slope_by_coalition.csv")
if (file.exists(slope_csv)) {
  sl <- fread(slope_csv)
} else {
  era5 <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)
  nc_of <- function(bit, id) file.path(PIPE$MUSICA_DIR_HOBO, bit, sprintf("musica_out_HOBO_%s.nc", id))
  rows <- list()
  for (bit in PIPE$BITS) for (id in clu$id_plot) {
    f <- nc_of(bit, id); if (!file.exists(f)) next
    s <- tryCatch(slope_at_1m(f, era5), error = function(e) NA_real_)
    if (is.finite(s)) rows[[length(rows)+1L]] <- data.table(id_plot = id, bit = bit, slope_sim = s)
  }
  sl <- rbindlist(rows); fwrite(sl, slope_csv)
  cat(sprintf("extracted slope for %d (plot,coalition) pairs\n", nrow(sl)))
}

# ---- per-cluster addition order = ΔTmax mean|phi| order (same as ΔTmax forward)
M <- rbindlist(lapply(list.files("out_files/Chapter1/tables/metrics_parts","part_.*csv$",full.names=TRUE), fread))
M <- M[metric == "Tmax_all"]; M[, Cluster := relabel_cluster(Cluster)]
ord_of <- function(dt) names(sort(sapply(Fv, function(t) mean(abs(dt[[t]]), na.rm=TRUE)), decreasing=TRUE))
orders <- lapply(setNames(paste0("P",1:4), paste0("P",1:4)), function(cl) ord_of(M[Cluster==cl]))
orders[["All"]] <- ord_of(M)
mk <- function(a) { b <- rep("0",4L); b[a] <- "1"; paste(b, collapse="") }

build_row <- function(rowname, order, ids) {
  active <- integer(0); steps <- list(list(bit="0000", add="Baseline", step=0L))
  for (k in seq_along(order)) { active <- c(active, pos[[order[k]]])
    steps[[k+1L]] <- list(bit=mk(active), step=k,
                          add=if (k==length(order)) sprintf("+%s (REF)", order[k]) else sprintf("+%s", order[k])) }
  rbindlist(lapply(steps, function(s) {
    d <- merge(sl[bit==s$bit & id_plot %in% ids, .(id_plot, slope_sim)], obs, by="id_plot")
    d[, `:=`(rowlab=rowname, step=s$step, add=s$add)]; d }))
}
rows <- rbindlist(c(lapply(paste0("P",1:4), function(cl) build_row(cl, orders[[cl]], clu[Cluster==cl, id_plot])),
                    list(build_row("All", orders[["All"]], clu$id_plot))))
rows[, rowlab := factor(rowlab, levels=c(paste0("P",1:4),"All"))]

ann <- rows[, { ok <- is.finite(slope_sim) & is.finite(slope_obs)
  .(r = if (sum(ok)>2) cor(slope_sim[ok], slope_obs[ok]) else NA_real_,
    RMSE = sqrt(mean((slope_sim[ok]-slope_obs[ok])^2, na.rm=TRUE)), n = sum(ok)) }, by=.(rowlab, step)]
ann <- merge(ann, unique(rows[,.(rowlab,step,add)]), by=c("rowlab","step"))
ann[, trait := sub(" \\(REF\\)","", sub("^\\+","", add))][add=="Baseline", trait := ""]
ann[, is_lai := trait=="LAI"]
fwrite(ann[order(rowlab,step)], file.path(PIPE$OUT_TAB, "tab_forward_slope.csv"))
cat("\n=== r (sim vs obs SLOPE) per row × step ===\n"); print(dcast(ann, rowlab~step, value.var="r")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,2) else x)])

pal <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850", All="grey30")
mk_panel <- function(metric, ylab, ylim=NULL, drop_baseline=FALSE) {
  dd <- copy(ann); dd[, val := get(metric)]; if (drop_baseline) dd <- dd[is.finite(val)]
  set.seed(1)
  g <- ggplot(dd, aes(step, val, colour=rowlab, group=rowlab)) +
    geom_line(linewidth=0.9) +
    geom_point(data=dd[!(is_lai)], size=2) + geom_point(data=dd[(is_lai)], shape=18, size=4.2) +
    geom_text_repel(aes(label=trait), size=2.8, fontface="bold", show.legend=FALSE,
                    max.overlaps=Inf, min.segment.length=0, segment.size=0.2, box.padding=0.25) +
    scale_colour_manual(values=pal, name=NULL) +
    scale_x_continuous(breaks=0:4, labels=c("Baseline","step 1","step 2","step 3","step 4"),
                       expand=expansion(mult=c(0.04,0.04))) +
    labs(x=NULL, y=ylab) + theme_bw(base_size=12) +
    theme(panel.grid.minor=element_blank(), legend.position="right")
  g + coord_cartesian(xlim=c(-0.15,4.15), ylim=ylim)
}
pr  <- mk_panel("r", "r (sim vs obs slope)", ylim=c(0.5,1.0), drop_baseline=TRUE) + scale_y_continuous(breaks=seq(0.5,1,0.1)) +
  annotate("text", x=1.0, y=0.52, hjust=0, size=2.7, colour="#1A9850",
           label="P4: +LAD alone r = 0.11 (off scale); LAI rescues it")
prm <- mk_panel("RMSE", "RMSE (slope, dimensionless)")
p <- (pr / prm) + plot_layout(guides="collect") +
  plot_annotation(caption="Micro-macro slope (1 m, hourly vs ERA5). Variables added in each cluster's ΔTmax-importance order. Diamond = LAI step.") &
  theme(legend.position="right", plot.caption=element_text(size=9, colour="grey40"))
ggsave("out_files/Chapter1/figures/FigSh_forward_slope.png", p, width=9.5, height=7.2, dpi=200, bg="white")
cat("DONE -> out_files/Chapter1/figures/FigSh_forward_slope.png\n")
