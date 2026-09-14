# ==============================================================================
# Chapter 3 — compare the monthly between-plot slope R² (slope-and-equilibrium
# thermal coupling, lm(T_micro~T_macro) per plot/month) between MuSICA v3.2.0
# (existing nc) and v3.2.3 + ABL_flag='iter' (nc_v323iter), for the 8 scenarios.
# IDENTICAL extraction pipeline for both binaries; the only difference is which
# nc directory is read. Reads the April->Nov window; the SOLID anchor is the
# summer plateau (Jun-Sep), April being n~30 (directional). Outputs a tidy table
# and a faceted figure.
#   Rscript cmp_v320_v323iter_ch3.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr); library(data.table)
  library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

ROOT_V320 <- file.path(CFG_C3$out_dir, "nc")           # existing v3.2.0 article binary
ROOT_V323 <- file.path(CFG_C3$out_dir, "nc_v323iter")  # v3.2.3 + iter
SCN <- c("STATIC_ALS","STATIC_ALS_DOPT","DYN_S2_ATBD_NM","DYN_S2_OPT_NM",
         "DYN_ATBD_rfull_NM","DYN_OPT_rdopt_NM","DYN_ALS_S2TIMING_NM","DYN_ALSoOPT_S2TIMING_NM",
         "FUSION_H","DYN_RF","DYN_S2_ANNUAL","CONST_ALS","NAIVE_S2_FORMSH")
LAB <- c(STATIC_ALS="LiDAR full", STATIC_ALS_DOPT="LiDAR d_opt",
         DYN_S2_ATBD_NM="S2 ATBD", DYN_S2_OPT_NM="S2 opt",
         DYN_ATBD_rfull_NM="S2 ATBD x->full", DYN_OPT_rdopt_NM="S2 opt x->d_opt",
         DYN_ALS_S2TIMING_NM="S2t.ALS/ATBD", DYN_ALSoOPT_S2TIMING_NM="S2t.ALS/opt",
         FUSION_H="Fusion-H (layered)", DYN_RF="S2 RF (dyn)", DYN_S2_ANNUAL="S2 annual (dyn)",
         CONST_ALS="LiDAR const", NAIVE_S2_FORMSH="S2 naive")

mlab <- function(d) format(d, "%Y-%m")
WIN  <- as.Date(c("2021-04-15","2021-11-15")); ds_all <- seq(WIN[1], WIN[2], by="day")
slp  <- function(y,x) if (sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_

# ---- macro hourly + observed per-plot/month slope (version-independent) -------
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all][, `:=`(hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))]
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), Tobs=t_hobo, hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))]
hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
OBS <- hb[, .(so=slp(Tobs,Tmacro)), by=.(id_plot,mo)]

read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  as.data.table(res)[as.Date(time) %in% ds_all, .(Tsim=Tair_sim, hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))] }

extract_version <- function(root, version){
  out <- list()
  for (scn in SCN){
    d <- file.path(root, scn); fs <- list.files(d, pattern="\\.nc$", full.names=TRUE)
    if (!length(fs)) { cat(sprintf("  [%s] %s: no nc\n", version, scn)); next }
    sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)")
      s<-read_sim(f); if(is.null(s))return(NULL); s[,id_plot:=id][] }), fill=TRUE)
    if (!nrow(sim)) next
    smh <- merge(sim, macro_h[,.(hr,Tmacro)], by="hr")
    SIM <- smh[, .(ss=slp(Tsim,Tmacro)), by=.(id_plot,mo)]
    M <- merge(SIM, OBS, by=c("id_plot","mo"))
    r <- M[, .(version=version, scenario=scn, n=.N,
               sl_r2 = if (.N>2) cor(ss,so,use="complete.obs")^2 else NA_real_), by=mo][order(mo)]
    out[[scn]] <- r; cat(sprintf("  [%s] %s: %d months, %d nc\n", version, scn, nrow(r), length(fs)))
  }
  rbindlist(out)
}
cat("=== extracting v3.2.0 ===\n"); A <- extract_version(ROOT_V320, "v3.2.0")
cat("=== extracting v3.2.3 iter ===\n"); B <- extract_version(ROOT_V323, "v3.2.3 iter")
D <- rbind(A, B)
D[, month := factor(mlab2 <- substr(mo,6,7))]
D[, label := factor(LAB[scenario], levels=unname(LAB))]
D[, version := factor(version, levels=c("v3.2.0","v3.2.3 iter"))]
fwrite(D, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table_v320_v323iter_ch3.csv")

# ---- plateau (Jun-Sep) vs April summary --------------------------------------
D[, fam := fifelse(grepl("STATIC|CONST_ALS", scenario), "LiDAR",
            fifelse(grepl("S2TIMING|FUSION_H", scenario), "Fusion", "S2"))]
plateau <- D[substr(mo,6,7) %in% c("06","07","08","09"),
             .(plateau_R2=mean(sl_r2,na.rm=TRUE)), by=.(version,scenario,label,fam)]
april   <- D[substr(mo,6,7)=="04", .(april_R2=sl_r2, april_n=n), by=.(version,scenario)]
S <- merge(plateau, april, by=c("version","scenario"))[order(scenario,version)]
cat("\n=== plateau (Jun-Sep mean R²) + April R² ===\n"); print(S[,.(version,label,fam,plateau_R2=round(plateau_R2,3),april_R2=round(april_R2,3),april_n)])
fwrite(S, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table_v320_v323iter_ch3_summary.csv")

# ---- figure: monthly slope R², version x scenario ----------------------------
D[, mo_short := factor(month.abb[as.integer(substr(mo,6,7))], levels=month.abb)]
g <- ggplot(D[!is.na(sl_r2)], aes(mo_short, sl_r2, colour=version, group=version)) +
  geom_line(linewidth=0.6) + geom_point(size=1.4) +
  facet_wrap(~label, ncol=4) +
  scale_colour_manual(values=c("v3.2.0"="#0072B2","v3.2.3 iter"="#D55E00"), name=NULL) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="2021", y="Between-plot slope R² (per month)",
       subtitle="Monthly between-plot ranking skill (slope-and-equilibrium) per scenario. Solid anchor = summer plateau (Jun-Sep); April n~30 directional. v3.2.3 over-smooths the field (c1: ΔTmax r 0.93->0.63).") +
  theme_article(10) + theme(legend.position="top")
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_v320_v323iter_monthly_slopeR2"), g, 10, 5.4)
cat("\nDONE -> Table_v320_v323iter_ch3{,_summary}.csv + Fig_v320_v323iter_monthly_slopeR2\n")
