# ==============================================================================
# Chapter 3 — two free diagnostics (no new sim) gating the leaf-out reframe:
#  (1) DATE-vs-TRAJECTORY mechanism: is the per-plot leaf-out advantage of DYN
#      over the LiDAR-magnitude parametric curve (STATIC_ALS) concentrated in
#      plots whose S2 leaf-out date departs most from the parametric DOY 115?
#  (2) d_opt confound: partial Spearman(frac_below_dopt, ΔRMSE | total LAI).
# Run:  Rscript c3_freechecks.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
# partial Spearman of (x,y) controlling z: Pearson on rank-residuals.
pcor_spearman <- function(x, y, z) {
  rx<-rank(x); ry<-rank(y); rz<-rank(z)
  ex<-residuals(lm(rx~rz)); ey<-residuals(lm(ry~rz))
  ct<-cor.test(ex, ey)
  list(estimate=unname(ct$estimate), p.value=ct$p.value)
}
src <- list.files("R", pattern="\\.R$", full.names=TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df_hobo <- prep$df_plots; ts_list <- prep$ts_by_plot

REF <- "DYN_S2_ANNUAL"; PAR <- "STATIC_ALS"   # both LiDAR magnitude; differ in pheno
want <- c(REF, PAR, "CONST_ALS", "STATIC_S2_ATBD")
all_sc <- make_all_scenarios_c3(ts_list, CFG_C3$list_year, CFG_C3$d_opt_m, CFG_C3$scenarios_mode)
nc_parent <- file.path(CFG_C3$out_dir, "nc")
sc <- Filter(function(s) s$name %in% want &&
       length(list.files(file.path(nc_parent, s$name), "\\.nc$"))>0, all_sc)

ds <- seq(as.Date("2021-04-01"), as.Date("2021-05-31"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
val <- validate_scenarios_at_hobos(df_hobo, hd, sc, nc_parent, dm, ds,
                                   CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
pr <- val$daily %>% group_by(scenario,id_plot) %>%
  summarise(rmse=sqrt(mean((Delta_sim-Delta_obs)^2,na.rm=TRUE)), .groups="drop") %>%
  pivot_wider(names_from=scenario, values_from=rmse)

# per-plot S2 leaf-out DOY (half-max up) from annual series
ann <- ts_list$annual
lo_doy <- function(s){d<-s$doy;l<-s$lai;mid<-min(l)+0.5*(max(l)-min(l));d[which(l>=mid)[1]]}
key <- df_hobo %>% mutate(pid=sprintf("X%d_Y%d",round(x),round(y))) %>%
  dplyr::select(id_plot, pid, LAI, LAI_ALS, Hmax)
s2lo <- data.frame(pid=names(ann), s2_leafout=sapply(ann, lo_doy))
key <- key %>% left_join(s2lo, by="pid")

# below-d_opt LAD fraction
lad_cols <- grep("^LAD_Layer_", names(df_hobo), value=TRUE)
lad_h <- as.numeric(sub("^LAD_Layer_","",lad_cols)); below <- lad_h < CFG_C3$d_opt_m
LADm <- as.matrix(df_hobo[,lad_cols]); LADm[!is.finite(LADm)]<-0
key$frac_below <- rowSums(LADm[,below,drop=FALSE])/pmax(rowSums(LADm),1e-9)
key$total_lai  <- ifelse(is.finite(key$LAI_ALS), key$LAI_ALS, key$LAI)

D <- pr %>% inner_join(key, by="id_plot") %>%
  mutate(adv_vs_param = STATIC_ALS - DYN_S2_ANNUAL,      # >0 = DYN better at leaf-out
         dRMSE_dyn_param = DYN_S2_ANNUAL - pmin(CONST_ALS,STATIC_ALS,STATIC_S2_ATBD),
         mismatch = abs(s2_leafout - 115)) %>%
  filter(is.finite(adv_vs_param), is.finite(mismatch))

cat(sprintf("n=%d plots\n", nrow(D)))

# ---- CHECK 1: date-vs-trajectory ----
cat("\n=== CHECK 1: leaf-out advantage (DYN - STATIC_ALS) ~ |S2 leaf-out - DOY115| ===\n")
cat(sprintf("DYN better than STATIC_ALS at leaf-out in %d/%d plots; mean adv=%+.3f °C\n",
    sum(D$adv_vs_param>0), nrow(D), mean(D$adv_vs_param)))
rho1 <- cor(D$mismatch, D$adv_vs_param, method="spearman")
set.seed(41); b1 <- replicate(5000,{i<-sample(nrow(D),replace=TRUE)
  cor(D$mismatch[i], D$adv_vs_param[i], method="spearman")})
ci1 <- quantile(b1,c(.025,.975),na.rm=TRUE)
cat(sprintf("Spearman(mismatch, advantage) = %.3f  CI95[%.3f,%.3f]  (%s)\n",
    rho1, ci1[1], ci1[2], if(ci1[1]>0|ci1[2]<0)"SIG -> DATE effect" else "ns"))
# split: advantage in high vs low mismatch
md <- median(D$mismatch)
cat(sprintf("  mean advantage: low-mismatch(<=%.0f) = %+.3f | high-mismatch(>%.0f) = %+.3f\n",
    md, mean(D$adv_vs_param[D$mismatch<=md]), md, mean(D$adv_vs_param[D$mismatch>md])))

# ---- CHECK 2: d_opt confound ----
cat("\n=== CHECK 2: d_opt result, raw vs partialling out total LAI ===\n")
raw <- cor(D$frac_below, D$dRMSE_dyn_param, method="spearman")
cat(sprintf("raw   Spearman(frac_below_dopt, ΔRMSE_dyn-param) = %.3f\n", raw))
cat(sprintf("cor(frac_below_dopt, total_LAI) = %.3f\n",
    cor(D$frac_below, D$total_lai, method="spearman")))
pc <- pcor_spearman(D$frac_below, D$dRMSE_dyn_param, D$total_lai)
cat(sprintf("partial Spearman | total_LAI = %.3f  (p=%.3f)  -> %s\n",
    pc$estimate, pc$p.value, if(pc$p.value<0.05)"survives" else "DIES (confounded by LAI)"))
pcH <- pcor_spearman(D$frac_below, D$dRMSE_dyn_param, D$Hmax)
cat(sprintf("partial Spearman | Hmax     = %.3f  (p=%.3f)\n", pcH$estimate, pcH$p.value))

write.csv(D %>% dplyr::select(id_plot,DYN_S2_ANNUAL,STATIC_ALS,adv_vs_param,
          s2_leafout,mismatch,frac_below,total_lai,Hmax,dRMSE_dyn_param),
          "/home/corroyez/Documents/NC_Full/output/tables/Table5b_freechecks.csv", row.names=FALSE)
cat("\nDONE\n")
