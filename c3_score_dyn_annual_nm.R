# ==============================================================================
# c3_score_dyn_annual_nm.R
#
# Score the re-run dynamic scenario against the observations, in the same frame
# as Table_archetype_scores_CHS41.csv, and report what changes.
#
# Extraction follows c3_tm_extract_CHS41.R: sub-canopy air at 1 m, read at the
# hour of the macroclimatic daily maximum, averaged over June-September 2021.
#
# Input : out_files/Chapter3_CHS41/nc_dyn_annual_nm/DYN_S2_ANNUAL/*.nc
# Output: out_files/Chapter3_CHS41/tables/perplot_dtmax_dyn_annual_nm.csv
#   Rscript c3_score_dyn_annual_nm.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr); library(data.table)
  library(parallel); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
source("Chapter3_config_CHS41.R")

ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
FORC <- CFG_C3$forcing_file
MREF <- macro_ref(FORC, ds)
NCD  <- file.path(CFG_C3$out_dir, "nc_dyn_annual_nm", "DYN_S2_ANNUAL")
OUT  <- file.path(CFG_C3$out_dir, "tables", "perplot_dtmax_dyn_annual_nm.csv")

one <- function(f) {
  id <- str_extract(basename(f), "(?<=HOBO_).*(?=\\.nc)")
  nc <- tryCatch(nc_open(f), error = function(e) NULL)
  if (is.null(nc)) return(NULL)
  tr <- tryCatch(get_tair_at_z(nc, 1.0), error = function(e) NULL)
  nc_close(nc)
  if (is.null(tr) || !nrow(tr)) return(NULL)
  m <- as.data.table(tr)[, .(time = floor_date(time, "hour"), Tmic = Tair_sim)]
  m <- m[, .(Tmic = mean(Tmic, na.rm = TRUE)), by = time][as.Date(time) %in% ds]
  data.table(id_plot = id, d_new = delta_tmax_mean(m, MREF, ds))
}
fs  <- list.files(NCD, "\\.nc$", full.names = TRUE)
sim <- rbindlist(mclapply(fs, one, mc.cores = 4), fill = TRUE)
sim <- sim[is.finite(d_new)]
cat(sprintf("placettes extraites : %d / %d\n", nrow(sim), length(fs)))
fwrite(sim, OUT)

# ---- score by archetype, same recipe as c3_archetype_scores_CHS41.R ----------
cl <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")
D  <- merge(sim, cl[, .(id_plot, obs = dTmax_obs, P)], by = "id_plot")

N_BOOT <- 2000L; BOOT_SEED <- 123L
score <- function(obs, s) {
  ok <- is.finite(obs) & is.finite(s); o <- obs[ok]; v <- s[ok]; n <- length(o)
  if (n < 3) return(list(n = n, R2 = NA_real_, R2_lo = NA_real_, R2_hi = NA_real_))
  set.seed(BOOT_SEED)
  bs <- replicate(N_BOOT, { i <- sample.int(n, n, replace = TRUE)
    if (sd(o[i]) == 0 || sd(v[i]) == 0) NA_real_ else cor(o[i], v[i])^2 })
  ci <- unname(quantile(bs, c(.025, .975), na.rm = TRUE))
  list(n = n, R2 = cor(o, v)^2, R2_lo = ci[1], R2_hi = ci[2])
}
new <- rbind(D[, score(obs, d_new), by = .(group = P)],
             D[, score(obs, d_new)][, group := "All (53)"], fill = TRUE)
setorder(new, group)

old <- fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table_archetype_scores_CHS41.csv")
old <- old[scenario == "DYN_S2_ANNUAL", .(group, n_old = n, R2_old = R2)]
cmp <- merge(new, old, by = "group", all.x = TRUE)
cat("\n=== Combinaison dynamique : avant / apres ===\n")
print(cmp[, .(archetype = group, n_avant = n_old, n_apres = n,
              R2_avant = round(R2_old, 3), R2_apres = round(R2, 3),
              IC = sprintf("[%.2f, %.2f]", R2_lo, R2_hi))])

# ---- write the new DYN_S2_ANNUAL rows into the chapter's score table --------
# The other four scenarios are untouched: only the dynamic combination was
# re-run. The shipped table is backed up before being rewritten.
score_full <- function(obs, s) {
  ok <- is.finite(obs) & is.finite(s); o <- obs[ok]; v <- s[ok]; n <- length(o)
  if (n < 3) return(list(n = n, R2 = NA_real_, R2_lo = NA_real_, R2_hi = NA_real_,
                         bias = NA_real_, RMSE = NA_real_))
  set.seed(BOOT_SEED)
  bs <- replicate(N_BOOT, { i <- sample.int(n, n, replace = TRUE)
    if (sd(o[i]) == 0 || sd(v[i]) == 0) NA_real_ else cor(o[i], v[i])^2 })
  ci <- unname(quantile(bs, c(.025, .975), na.rm = TRUE))
  list(n = n, R2 = cor(o, v)^2, R2_lo = ci[1], R2_hi = ci[2],
       bias = mean(v - o), RMSE = sqrt(mean((v - o)^2)))
}
rows <- rbind(D[, score_full(obs, d_new), by = .(group = P)],
              D[, score_full(obs, d_new)][, group := "All (53)"], fill = TRUE)
rows[, `:=`(label = "LiDAR x S2 (dynamic annual)", scenario = "DYN_S2_ANNUAL")]

P_TAB <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table_archetype_scores_CHS41.csv"
tab   <- fread(P_TAB)
fwrite(tab, sub("\\.csv$", "_backup_2026-08-31.csv", P_TAB))
tab   <- tab[scenario != "DYN_S2_ANNUAL"]
tab   <- rbind(tab, rows[, names(tab), with = FALSE])
setorder(tab, label, group)
fwrite(tab, P_TAB)
cat("\ntable mise a jour :", P_TAB, "\n")
print(tab[scenario == "DYN_S2_ANNUAL", .(group, n, R2 = round(R2, 3),
                                         bias = round(bias, 2), RMSE = round(RMSE, 2))])
