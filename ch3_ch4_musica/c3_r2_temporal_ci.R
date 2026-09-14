# ==============================================================================
# Chapter 3 — follow-up: (A) block-bootstrap ΔR² (by plot) fullyear, so the R²
# ranking is tested not just ranked; (B) block-bootstrap (by plot) CI of the
# two-way-demeaned temporal LAI<->buffering r, plus the pure-temporal-anomaly
# variance share that reconciles the annex with the forcing null.
# NO re-simulation (cached nc). Run:  Rscript c3_r2_temporal_ci.R
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
want <- c("CONST_ALS","DYN_RF","STATIC_S2_ATBD","STATIC_ALS","DYN_S2_ATBD",
          "DYN_S2_ANNUAL","NAIVE_S2_FORMSH")
all_sc <- make_all_scenarios_c3(ts_list, list_year = CFG_C3$list_year,
                                d_opt = CFG_C3$d_opt_m, mode = CFG_C3$scenarios_mode)
nc_parent <- file.path(CFG_C3$out_dir, "nc")
sc <- Filter(function(s) s$name %in% want &&
       dir.exists(file.path(nc_parent, s$name)) &&
       length(list.files(file.path(nc_parent, s$name), "\\.nc$")) > 0, all_sc)

ds <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
val <- validate_scenarios_at_hobos(df_hobo, hd, sc, nc_parent, dm, ds,
                                   CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
d <- val$daily

# ---- (A) pooled R² per scenario + block-bootstrap ΔR² (resample plots) -------
pooled_r2 <- function(df) {
  df <- df[is.finite(df$Delta_obs) & is.finite(df$Delta_sim), ]
  if (nrow(df) < 3) return(NA_real_)
  cor(df$Delta_obs, df$Delta_sim)^2
}
plots <- unique(d$id_plot)
by_plot <- split(d, d$id_plot)            # list of per-plot daily frames

r2_obs <- sapply(want, function(s) pooled_r2(d[d$scenario == s, ]))

set.seed(11); B <- 5000
# bootstrap distribution of pooled R² per scenario, resampling PLOTS (blocks)
boot_r2 <- function(s) {
  replicate(B, {
    samp <- sample(plots, replace = TRUE)
    dd <- bind_rows(lapply(samp, function(p) by_plot[[p]]))
    pooled_r2(dd[dd$scenario == s, ])
  })
}
# Paired ΔR² vs REF: use the SAME resampled plot set each iteration.
boot_dR2 <- function(s) {
  replicate(B, {
    samp <- sample(plots, replace = TRUE)
    dd <- bind_rows(lapply(samp, function(p) by_plot[[p]]))
    pooled_r2(dd[dd$scenario == s, ]) - pooled_r2(dd[dd$scenario == REF, ])
  })
}
cat("=== fullyear pooled R² + ΔR² vs DYN_S2_ANNUAL (block-bootstrap by plot, B=5000) ===\n")
r2_tab <- data.frame()
for (s in want) {
  if (s == REF) { cat(sprintf("  %-16s R2=%.3f  (reference)\n", s, r2_obs[s]))
    r2_tab <- rbind(r2_tab, data.frame(scenario=s, R2=round(r2_obs[s],3),
      dR2=0, ci_lo=0, ci_hi=0, sig="ref")); next }
  bd <- boot_dR2(s); ci <- quantile(bd, c(.025,.975), na.rm = TRUE)
  sig <- if (ci[1] > 0 | ci[2] < 0) "SIG" else "ns"
  cat(sprintf("  %-16s R2=%.3f  ΔR2 vs REF=%+.3f CI95[%+.3f,%+.3f] %s\n",
              s, r2_obs[s], mean(bd, na.rm=TRUE), ci[1], ci[2], sig))
  r2_tab <- rbind(r2_tab, data.frame(scenario=s, R2=round(r2_obs[s],3),
    dR2=round(mean(bd,na.rm=TRUE),4), ci_lo=round(ci[1],4),
    ci_hi=round(ci[2],4), sig=sig))
}
write.csv(r2_tab, "/home/corroyez/Documents/NC_Full/output/tables/Table4b_dR2_bootstrap.csv", row.names=FALSE)
write.csv(r2_tab, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4b_dR2_bootstrap.csv", row.names=FALSE)

# ---- (B) temporal annex: two-way r + block-bootstrap by plot + var share -----
ann <- ts_list$annual
key <- df_hobo %>% mutate(pid = sprintf("X%d_Y%d", round(x), round(y))) %>%
  dplyr::select(id_plot, pid)
obs <- d %>% filter(scenario == REF) %>%
  dplyr::select(id_plot, date, Delta_obs) %>% mutate(doy = yday(as.Date(date))) %>%
  left_join(key, by = "id_plot") %>% filter(pid %in% names(ann))
lai_long <- bind_rows(lapply(unique(obs$pid), function(p)
  data.frame(pid = p, doy = ann[[p]]$doy, lai = ann[[p]]$lai)))
dat <- obs %>% inner_join(lai_long, by = c("pid","doy")) %>%
  filter(is.finite(lai), is.finite(Delta_obs))

twoway <- function(df) {
  gl <- mean(df$lai); gd <- mean(df$Delta_obs)
  dm_ <- df %>% group_by(doy) %>% summarise(ld=mean(lai), dd=mean(Delta_obs), .groups="drop")
  pm_ <- df %>% group_by(pid) %>% summarise(lp=mean(lai), dp=mean(Delta_obs), .groups="drop")
  x <- df %>% left_join(dm_, "doy") %>% left_join(pm_, "pid") %>%
    mutate(lai_tw = lai-ld-lp+gl, del_tw = Delta_obs-dd-dp+gd)
  list(r = cor(x$lai_tw, x$del_tw, use="complete.obs"),
       vshare = var(x$lai_tw)/var(df$lai))
}
base <- twoway(dat)
pids <- unique(dat$pid)
by_pid <- split(dat, dat$pid)
set.seed(13)
br <- replicate(2000, {
  samp <- sample(pids, replace = TRUE)
  dd <- bind_rows(lapply(seq_along(samp), function(i){
    z <- by_pid[[samp[i]]]; z$pid <- paste0(z$pid,"_",i); z }))  # unique block ids
  twoway(dd)$r
})
ci <- quantile(br, c(.025,.975), na.rm = TRUE)
cat(sprintf("\n=== temporal two-way-demeaned r (Blois, %d plots) ===\n", length(pids)))
cat(sprintf("  r = %.3f  block-bootstrap CI95[%.3f, %.3f]  (%s)\n",
            base$r, ci[1], ci[2], if (ci[1]>0|ci[2]<0) "SIG" else "ns"))
cat(sprintf("  pure-temporal-anomaly LAI variance share = %.3f  (i.e. %.1f%% of LAI\n",
            base$vshare, 100*base$vshare))
cat("    temporal variance survives plot+day demeaning -> the channel that can act\n")
cat("    on the forcing is small, reconciling the annex with the RMSE/R² nulls.)\n")
annex2 <- data.frame(
  quantity = c("two_way_r","two_way_r_ci_lo","two_way_r_ci_hi",
               "two_way_LAI_var_share","n_plots_blois"),
  value = c(round(base$r,3), round(ci[1],3), round(ci[2],3),
            round(base$vshare,3), length(pids)))
write.csv(annex2, "/home/corroyez/Documents/NC_Full/output/tables/TableA1_temporal_r_annex.csv", row.names=FALSE)
write.csv(annex2, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/TableA1_temporal_r_annex.csv", row.names=FALSE)
cat("\nDONE\n")
