# ==============================================================================
# Chapter 3 — regime test: does the S2 dynamic beat the parametric curve at
# LEAF-OUT (where 2021 timing departs from the generic curve)?
#   (A) leaf-out & autumn: paired-bootstrap ΔRMSE + block-bootstrap ΔR²
#       of DYN_S2_ANNUAL vs {CONST_ALS, STATIC_ALS, STATIC_S2_ATBD}.
#   (B) PRIMARY pre-registered test: per-plot ΔRMSE(DYN − best parametric) in
#       leaf-out ~ below-d_opt LAD fraction (structural, S2-independent).
# NO re-simulation (cached nc). Run:  Rscript c3_leafout_regime.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
source("Chapter3_config.R")

prep <- load_lai_prep(CFG_C3); df_hobo <- prep$df_plots; ts_list <- prep$ts_by_plot
REF  <- "DYN_S2_ANNUAL"
cmp  <- c("CONST_ALS","STATIC_ALS","STATIC_S2_ATBD")
want <- c(REF, cmp)
all_sc <- make_all_scenarios_c3(ts_list, list_year = CFG_C3$list_year,
                                d_opt = CFG_C3$d_opt_m, mode = CFG_C3$scenarios_mode)
nc_parent <- file.path(CFG_C3$out_dir, "nc")
sc <- Filter(function(s) s$name %in% want &&
       dir.exists(file.path(nc_parent, s$name)) &&
       length(list.files(file.path(nc_parent, s$name), "\\.nc$")) > 0, all_sc)

WIN <- list(leafout = seq(as.Date("2021-04-01"), as.Date("2021-05-31"), by="day"),
            autumn  = seq(as.Date("2021-10-01"), as.Date("2021-10-31"), by="day"))

pooled_r2 <- function(df){df<-df[is.finite(df$Delta_obs)&is.finite(df$Delta_sim),]
  if(nrow(df)<3) NA_real_ else cor(df$Delta_obs,df$Delta_sim)^2}
set.seed(21); B <- 5000

daily_by_win <- list()
for (wn in names(WIN)) {
  ds <- WIN[[wn]]
  dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
  hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
  val <- validate_scenarios_at_hobos(df_hobo, hd, sc, nc_parent, dm, ds,
                                     CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
  d <- val$daily; daily_by_win[[wn]] <- d
  pr <- d %>% group_by(scenario,id_plot) %>%
    summarise(rmse=sqrt(mean((Delta_sim-Delta_obs)^2,na.rm=TRUE)), .groups="drop")
  W <- pr %>% pivot_wider(names_from=scenario, values_from=rmse)
  W <- W[stats::complete.cases(W),]; np <- nrow(W)
  plots <- unique(d$id_plot); by_plot <- split(d, d$id_plot)
  cat(sprintf("\n===== %s (n=%d plots) — ΔRMSE (paired) & ΔR² (block) vs %s =====\n",
              toupper(wn), np, REF))
  cat(sprintf("  %-16s RMSE=%.3f  R2=%.3f  (reference)\n",
              REF, mean(W[[REF]]), pooled_r2(d[d$scenario==REF,])))
  for (s in cmp) {
    dlt <- W[[s]] - W[[REF]]
    ci  <- quantile(replicate(B, mean(sample(dlt,replace=TRUE))), c(.025,.975))
    bd  <- replicate(B, {samp<-sample(plots,replace=TRUE)
      dd<-bind_rows(lapply(samp,function(p)by_plot[[p]]))
      pooled_r2(dd[dd$scenario==s,])-pooled_r2(dd[dd$scenario==REF,])})
    cir <- quantile(bd, c(.025,.975), na.rm=TRUE)
    cat(sprintf("  %-16s RMSE=%.3f ΔRMSE=%+.3f[%+.3f,%+.3f]%s | R2=%.3f ΔR2=%+.3f[%+.3f,%+.3f]%s\n",
        s, mean(W[[s]]), mean(dlt), ci[1], ci[2], if(ci[1]>0|ci[2]<0)"*"else"ns",
        pooled_r2(d[d$scenario==s,]), mean(bd,na.rm=TRUE), cir[1], cir[2],
        if(cir[1]>0|cir[2]<0)"*"else"ns"))
  }
}

# ---- (B) PRIMARY: leaf-out per-plot ΔRMSE(DYN − best parametric) ~ below-dopt frac
d_lo <- daily_by_win[["leafout"]]
pr_lo <- d_lo %>% group_by(scenario,id_plot) %>%
  summarise(rmse=sqrt(mean((Delta_sim-Delta_obs)^2,na.rm=TRUE)), .groups="drop") %>%
  pivot_wider(names_from=scenario, values_from=rmse)
# best parametric per plot = min RMSE among the 3 non-S2-dynamic comparators
pr_lo$param_best <- pmin(pr_lo$CONST_ALS, pr_lo$STATIC_ALS, pr_lo$STATIC_S2_ATBD, na.rm=TRUE)
pr_lo$dRMSE_dyn_minus_param <- pr_lo$DYN_S2_ANNUAL - pr_lo$param_best  # <0 = DYN better

# below-d_opt LAD fraction per plot (structural, S2-independent)
lad_cols <- grep("^LAD_Layer_", names(df_hobo), value=TRUE)
lad_h    <- as.numeric(sub("^LAD_Layer_", "", lad_cols))      # layer mid-height (m)
dopt     <- CFG_C3$d_opt_m
below    <- lad_h < dopt
LADm <- as.matrix(df_hobo[, lad_cols]); LADm[!is.finite(LADm)] <- 0
frac_below <- rowSums(LADm[, below, drop=FALSE]) / pmax(rowSums(LADm), 1e-9)
strat <- data.frame(id_plot = df_hobo$id_plot, frac_below_dopt = frac_below)

reg <- pr_lo %>% dplyr::select(id_plot, DYN_S2_ANNUAL, param_best,
                               dRMSE_dyn_minus_param) %>%
  inner_join(strat, by="id_plot") %>% filter(is.finite(dRMSE_dyn_minus_param))

cat(sprintf("\n===== PRIMARY TEST (leaf-out, n=%d) =====\n", nrow(reg)))
cat(sprintf("below-d_opt(%dm) LAD fraction: median %.2f IQR[%.2f,%.2f]\n",
    dopt, median(reg$frac_below_dopt), quantile(reg$frac_below_dopt,.25),
    quantile(reg$frac_below_dopt,.75)))
cat(sprintf("DYN beats best-parametric (ΔRMSE<0) in %d/%d plots (%.0f%%); mean ΔRMSE=%+.3f\n",
    sum(reg$dRMSE_dyn_minus_param<0), nrow(reg),
    100*mean(reg$dRMSE_dyn_minus_param<0), mean(reg$dRMSE_dyn_minus_param)))
rho <- suppressWarnings(cor(reg$frac_below_dopt, reg$dRMSE_dyn_minus_param, method="spearman"))
set.seed(31)
rb <- replicate(5000,{i<-sample(nrow(reg),replace=TRUE)
  suppressWarnings(cor(reg$frac_below_dopt[i], reg$dRMSE_dyn_minus_param[i], method="spearman"))})
cir <- quantile(rb, c(.025,.975), na.rm=TRUE)
cat(sprintf("Spearman(frac_below_dopt, ΔRMSE_dyn-param) = %.3f  boot CI95[%.3f,%.3f]  (%s)\n",
    rho, cir[1], cir[2], if(cir[1]>0|cir[2]<0)"SIG" else "ns"))
cat("  (>0 expected if more below-d_opt canopy -> S2-top phenology less representative\n")
cat("   -> DYN less advantageous; <0 if understory-leads-out dominates. Exploratory.)\n")

write.csv(reg, "/home/corroyez/Documents/NC_Full/output/tables/Table5_leafout_regime_dopt.csv", row.names=FALSE)
write.csv(reg, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table5_leafout_regime_dopt.csv", row.names=FALSE)
cat("\nDONE\n")
