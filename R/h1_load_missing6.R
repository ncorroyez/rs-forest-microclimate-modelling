# ==============================================================================
# Load the 6 cLHS coalitions missing from DT_daily_cLHS.rds (only 10/16 cached).
# Output : outputs/lovb/data/DT_daily_cLHS_missing6.rds
# Then aggregate + merge into DT_contrib_cLHS_full.rds (all 16 coalitions wide).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
})
source(here::here("R/lovb_01_load.R"))
source(here::here("R/lovb_02_aggregate.R"))

MISSING_BITS <- c("0011", "0101", "0110", "1001", "1010", "1100")

cli_h1("Load 6 missing cLHS coalitions (2400 NetCDF)")
DT_miss <- lovb_load_daily_clhs(
  coalitions = MISSING_BITS,
  cache_path = here::here("outputs/lovb/data/DT_daily_cLHS_missing6.rds"),
  use_cache  = TRUE
)
cli_alert_success("DT_daily_cLHS_missing6 ready : {nrow(DT_miss)} rows")

# Aggregate daily -> mean per (plot, bit_code)
cli_h2("Aggregate missing 6 -> Tmax_mean per (plot, coalition)")
DT_agg_miss <- DT_miss[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
                             Tmax_P90  = quantile(Delta_Tmax, 0.9,
                                                    na.rm = TRUE, type = 7),
                             n_days = .N),
                        by = .(x, y, bit_code)]
saveRDS(DT_agg_miss,
         here::here("outputs/lovb/data/DT_agg_cLHS_missing6.rds"))
cli_alert_success("Aggregated : {nrow(DT_agg_miss)} (plot, coalition) rows")

# Build wide table with all 16 coalitions per plot
cli_h2("Build wide DT_contrib_cLHS_full with all 16 coalitions per plot")
DT_existing <- readRDS(here::here("outputs/lovb/data/DT_agg_cLHS.rds"))
DT_all <- rbind(DT_existing[, .(x, y, bit_code, Tmax_mean, Tmax_P90,
                                   n_days = NULL)],
                  DT_agg_miss, fill = TRUE)
# Make sure no duplicate (plot, bit_code)
setkey(DT_all, x, y, bit_code)
DT_all <- unique(DT_all, by = c("x", "y", "bit_code"))
cli_alert("Combined rows : {nrow(DT_all)}  ({length(unique(DT_all$bit_code))} unique coalitions)")

# Pivot wide
DT_w <- dcast(DT_all, x + y ~ bit_code,
                value.var = c("Tmax_mean", "Tmax_P90"))
cli_alert("Wide DT : {nrow(DT_w)} plots x {ncol(DT_w)} cols")

# Merge with traits (LAI, Hmax, fCover, FPC1, Cluster) from existing contrib
DT_contrib <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
DT_traits <- DT_contrib[, .(x, y, Cluster, LAI, Hmax, fCover, FPC1)]
DT_full <- merge(DT_traits, DT_w, by = c("x", "y"))
saveRDS(DT_full,
         here::here("outputs/lovb/data/DT_contrib_cLHS_16coalitions.rds"))
cli_alert_success("Saved DT_contrib_cLHS_16coalitions.rds : {nrow(DT_full)} plots x 16 coalitions")
