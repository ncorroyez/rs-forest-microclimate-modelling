# ==============================================================================
# Phase 4 — Recompute attribution (Shapley/LOVB/LVA) for floor05 pipeline.
#
# Inputs (from Phase 3 in outputs/lovb_floor05/data/) :
#   DT_contrib_cLHS_floor05.rds              (400 plots, 10 coalitions wide)
#   DT_contrib_cLHS_floor05_16coalitions.rds (400 plots, 16 coalitions wide)
#   DT_contrib_archetypes_floor05.rds        (4 archetypes, 10 coalitions wide)
#   DT_daily_archetypes_floor05.rds          (used for Shapley archetype recompute)
#   DT_contrib_HOBO_floor05.rds              (53 sensors, 10 coalitions wide)
#
# Outputs in outputs/lovb_floor05/tables/ :
#   tab_LOVB_global_floor05.csv
#   tab_attribution_archetypes_pooled_floor05.csv
#   tab_attribution_cLHS_mean_floor05.csv
#   tab_consensus_4methods_floor05.csv
#   tab_HOBO_spearman_floor05.csv
#   tab_HOBO_spearman_dTmax_floor05.csv
# And cache :
#   DT_shapley_per_plot_floor05.rds
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

OUT_TAB  <- here::here("outputs/lovb_floor05/tables")
OUT_DATA <- here::here("outputs/lovb_floor05/data")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

cli_h1("Phase 4 — recompute attribution (floor05)")
t_start <- Sys.time()

VARS <- c("LAI", "Hmax", "fCover", "LAD")

# ---- 1. cLHS LOVB global (RMSE) and per-plot ---------------------------------
cli_h2("cLHS LOVB / LVA")
DT_c <- readRDS(file.path(OUT_DATA, "DT_contrib_cLHS_floor05.rds"))
glob_rows <- list()
for (v in VARS) for (m in c("mean","P90")) {
  d <- DT_c[[paste0("Delta_", v, "_", m)]]
  glob_rows[[length(glob_rows)+1L]] <- data.table(
    Variable = v, Metric = m,
    RMSE_global = sqrt(mean(d^2, na.rm = TRUE)),
    Bias_global = mean(d, na.rm = TRUE),
    pct_above_0.5C = round(100 * mean(abs(d) > 0.5, na.rm = TRUE), 2),
    N = sum(!is.na(d))
  )
}
tab_lovb_global <- rbindlist(glob_rows)
tab_lovb_global[, Ranking := frank(-RMSE_global, ties.method = "min"), by = Metric]
fwrite(tab_lovb_global, file.path(OUT_TAB, "tab_LOVB_global_floor05.csv"))
cli_alert("LOVB global (floor05) :"); print(tab_lovb_global)

# Attribution cLHS (mean) : Shapley + LOVB + LVA from DT_contrib + Shapley csv
DT_full <- readRDS(file.path(OUT_DATA, "DT_contrib_cLHS_floor05_16coalitions.rds"))
bits <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",
           "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
tcols <- paste0("Tmax_mean_", bits)

# Per-plot Shapley
shap_one_plot <- function(tmax_vec) {
  v_null <- as.numeric(tmax_vec["0000"])
  v <- function(b) as.numeric(tmax_vec[b]) - v_null
  make_bin <- function(bits_on, n = 4L) {
    b <- rep("0", n); b[bits_on] <- "1"; paste(b, collapse = "")
  }
  variables <- c("LAI", "Hmax", "fCover", "LAD")
  out <- numeric(4L)
  for (var_bit in seq_len(4L)) {
    others <- setdiff(seq_len(4L), var_bit)
    total  <- 0
    for (s in 0:3) {
      subsets <- if (s == 0L) list(integer(0)) else combn(others, s, simplify = FALSE)
      w <- factorial(s) * factorial(4 - s - 1) / factorial(4)
      for (S in subsets)
        total <- total + w * (v(make_bin(c(S, var_bit))) - v(make_bin(S)))
    }
    out[var_bit] <- total
  }
  setNames(out, variables)
}

cli_alert("Computing per-plot Shapley (16 coalitions, n=400 plots) ...")
rows <- vector("list", nrow(DT_full))
for (i in seq_len(nrow(DT_full))) {
  tmax_vec <- as.numeric(unlist(DT_full[i, ..tcols])); names(tmax_vec) <- bits
  if (any(is.na(tmax_vec))) next
  phi <- shap_one_plot(tmax_vec)
  rows[[i]] <- data.table(
    x = DT_full$x[i], y = DT_full$y[i], Cluster = DT_full$Cluster[i],
    LAI = DT_full$LAI[i], Hmax = DT_full$Hmax[i],
    fCover = DT_full$fCover[i], FPC1 = DT_full$FPC1[i],
    variable = names(phi), phi = unname(phi)
  )
}
DT_shap <- rbindlist(rows)
saveRDS(DT_shap, file.path(OUT_DATA, "DT_shapley_per_plot_floor05.rds"))
cli_alert_success("Saved DT_shapley_per_plot_floor05.rds ({length(unique(DT_shap$variable))*nrow(DT_full)} rows)")

# Aggregate per variable (cLHS global Shapley) : mean phi
shap_global_cLHS <- DT_shap[, .(Shapley_cLHS_value = mean(phi, na.rm = TRUE)),
                              by = variable]
shap_global_cLHS[, phi_abs := abs(Shapley_cLHS_value)]
shap_global_cLHS[, Shapley_cLHS_rank := frank(-phi_abs, ties.method = "min")]
shap_global_cLHS[, phi_abs := NULL]

# LOVB cLHS rank by RMSE_global mean
lovb_cLHS <- tab_lovb_global[Metric == "mean",
                                .(variable = Variable,
                                    LOVB_cLHS_value = RMSE_global)]
lovb_cLHS[, LOVB_cLHS_rank := frank(-LOVB_cLHS_value, ties.method = "min")]

# LVA cLHS : mean of |LVA_v_mean| across 400 plots
DT_c <- readRDS(file.path(OUT_DATA, "DT_contrib_cLHS_floor05.rds"))
lva_rows <- data.table(
  variable = VARS,
  LVA_cLHS_value = c(mean(abs(DT_c$LVA_LAI_mean),    na.rm = TRUE),
                       mean(abs(DT_c$LVA_Hmax_mean),   na.rm = TRUE),
                       mean(abs(DT_c$LVA_fCover_mean), na.rm = TRUE),
                       mean(abs(DT_c$LVA_LAD_mean),    na.rm = TRUE))
)
lva_rows[, LVA_cLHS_rank := frank(-LVA_cLHS_value, ties.method = "min")]

tab_attr_clhs <- merge(merge(shap_global_cLHS, lovb_cLHS, by = "variable"),
                         lva_rows, by = "variable")
tab_attr_clhs[, variable := factor(variable, levels = VARS)]
setorder(tab_attr_clhs, variable)
fwrite(tab_attr_clhs, file.path(OUT_TAB, "tab_attribution_cLHS_mean_floor05.csv"))
cli_alert("Attribution cLHS (floor05) :"); print(tab_attr_clhs)

# ---- 2. Archetypes attribution -----------------------------------------------
cli_h2("Archetypes attribution")
DT_a <- readRDS(file.path(OUT_DATA, "DT_contrib_archetypes_floor05.rds"))
# Need wide 16 coalitions per archetype for Shapley.
# DT_daily_archetypes_floor05 has 16 for C1 (from floor05) but only 10 for C2/C3/C4.
# Load the 6 missing coalitions for C2/C3/C4 from H1_archetypes/ (original sims).
source(here::here("R/lovb_01_load.R"))
DT_daily_a <- readRDS(file.path(OUT_DATA, "DT_daily_archetypes_floor05.rds"))
MISSING_BITS_ARCH <- c("0011","0101","0110","1001","1010","1100")
df_macro_arch <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
arch_extra <- list()
for (arch in c("Arch_C2", "Arch_C3", "Arch_C4")) {
  dir_a <- here::here("out_files/H1_archetypes", arch)
  for (bit in MISSING_BITS_ARCH) {
    sc_nm <- paste0("ARCH_", lovb_bit_to_suffix(bit))
    nc <- file.path(dir_a, sprintf("musica_out_%s_%s.nc", arch, sc_nm))
    if (!file.exists(nc) || file.size(nc) < 1e6) {
      cli_alert_warning("Missing arch NC : {basename(nc)}"); next
    }
    res <- tryCatch(extract_deltatmax_one(nc, df_macro_arch, CFG$date_seq),
                     error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0L) next
    dt <- as.data.table(res)
    dt[, archetype := arch]; dt[, bit_code := bit]
    arch_extra[[length(arch_extra)+1L]] <- dt
  }
}
DT_arch_extra <- rbindlist(arch_extra, fill = TRUE)
DT_arch_extra <- DT_arch_extra[, .(archetype, date, bit_code, Delta_Tmax)]
DT_daily_a <- rbind(DT_daily_a, DT_arch_extra, fill = TRUE)
cli_alert("Arch daily after loading 6 missing per C2/C3/C4 : {nrow(DT_daily_a)} rows ({length(unique(DT_daily_a$bit_code))} coalitions)")
DT_agg_a16 <- DT_daily_a[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE)),
                          by = .(archetype, bit_code)]
DT_w_a <- dcast(DT_agg_a16, archetype ~ bit_code, value.var = "Tmax_mean")
setnames(DT_w_a, bits, paste0("Tmax_mean_", bits))

# Per-archetype Shapley
arch_shap_rows <- list()
for (i in seq_len(nrow(DT_w_a))) {
  arch <- DT_w_a$archetype[i]
  tmax_vec <- as.numeric(unlist(DT_w_a[i, ..tcols])); names(tmax_vec) <- bits
  if (any(is.na(tmax_vec))) { cli_alert_warning("NAs for {arch}"); next }
  phi <- shap_one_plot(tmax_vec)
  arch_shap_rows[[arch]] <- data.table(archetype = arch,
                                          variable = names(phi),
                                          phi = unname(phi))
}
DT_shap_arch <- rbindlist(arch_shap_rows)
# Pooled (mean across 4 archetypes)
shap_pool <- DT_shap_arch[, .(Shapley_arch_value = mean(phi, na.rm = TRUE)),
                            by = variable]
shap_pool[, phi_abs := abs(Shapley_arch_value)]
shap_pool[, Shapley_arch_rank := frank(-phi_abs, ties.method = "min")]
shap_pool[, phi_abs := NULL]

# LOVB arch : mean |Delta_v_mean| across 4 archetypes
lovb_arch <- data.table(
  variable = VARS,
  LOVB_arch_value = c(mean(abs(DT_a$Delta_LAI_mean),    na.rm = TRUE),
                        mean(abs(DT_a$Delta_Hmax_mean),   na.rm = TRUE),
                        mean(abs(DT_a$Delta_fCover_mean), na.rm = TRUE),
                        mean(abs(DT_a$Delta_LAD_mean),    na.rm = TRUE))
)
lovb_arch[, LOVB_arch_rank := frank(-LOVB_arch_value, ties.method = "min")]

# LVA arch
lva_arch <- data.table(
  variable = VARS,
  LVA_arch_value = c(mean(abs(DT_a$LVA_LAI_mean),    na.rm = TRUE),
                       mean(abs(DT_a$LVA_Hmax_mean),   na.rm = TRUE),
                       mean(abs(DT_a$LVA_fCover_mean), na.rm = TRUE),
                       mean(abs(DT_a$LVA_LAD_mean),    na.rm = TRUE))
)
lva_arch[, LVA_arch_rank := frank(-LVA_arch_value, ties.method = "min")]

tab_attr_arch <- merge(merge(shap_pool, lovb_arch, by = "variable"),
                         lva_arch, by = "variable")
tab_attr_arch[, variable := factor(variable, levels = VARS)]
setorder(tab_attr_arch, variable)
fwrite(tab_attr_arch,
        file.path(OUT_TAB, "tab_attribution_archetypes_pooled_floor05.csv"))
cli_alert("Attribution archetypes (floor05) :"); print(tab_attr_arch)

# ---- 3. Consensus tables (Level 1 + Level 2 + 6-method) ----------------------
cli_h2("Consensus")
cons <- merge(tab_attr_arch[, .(variable, Shapley_arch_rank, LOVB_arch_rank, LVA_arch_rank)],
                tab_attr_clhs[, .(variable, Shapley_cLHS_rank, LOVB_cLHS_rank, LVA_cLHS_rank)],
                by = "variable")
cons[, mean_rank := rowMeans(.SD), .SDcols = c("Shapley_arch_rank","LOVB_arch_rank","LVA_arch_rank",
                                                  "Shapley_cLHS_rank","LOVB_cLHS_rank","LVA_cLHS_rank")]
cons[, consensus_rank := frank(mean_rank, ties.method = "min")]
fwrite(cons, file.path(OUT_TAB, "tab_consensus_6methods_floor05.csv"))
cli_alert("Consensus (6 methods, floor05) :"); print(cons)

# 4-method version (cLHS only)
cons4 <- tab_attr_clhs[, .(Variable = variable, Shapley_rank = Shapley_cLHS_rank,
                              LOVB_rank = LOVB_cLHS_rank, LVA_rank = LVA_cLHS_rank)]
cons4[, Consensus_mean := rowMeans(.SD),
        .SDcols = c("Shapley_rank","LOVB_rank","LVA_rank")]
cons4[, Consensus_rank := frank(Consensus_mean, ties.method = "min")]
fwrite(cons4, file.path(OUT_TAB, "tab_consensus_4methods_floor05.csv"))

# ---- 4. HOBO Spearman dTmax floor05 ------------------------------------------
cli_h2("HOBO spearman dTmax (floor05)")
DT_h <- readRDS(file.path(OUT_DATA, "DT_contrib_HOBO_floor05.rds"))
cg   <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
hobo <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
DT_h <- merge(DT_h, hobo, by = "id_plot")
DT_h[, dT_signed := -dTmax_obs]

spearman_boot_ci <- function(x, y, B = 1000L, seed = 42L) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  rho_hat <- suppressWarnings(cor(x, y, method = "spearman"))
  p_val   <- suppressWarnings(cor.test(x, y, method = "spearman",
                                          exact = FALSE)$p.value)
  set.seed(seed); rb <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    rb[b] <- suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }
  ci <- quantile(rb, c(0.025, 0.975), na.rm = TRUE, type = 7)
  list(rho = rho_hat, p = p_val, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
}

rows_sp <- list()
for (v in VARS) {
  DT_h[, LVA_v := get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
  for (method in c("LOVB","LVA","Shapley")) {
    if (method == "LOVB")    y <- DT_h[[paste0("Delta_", v, "_mean")]]
    if (method == "LVA")     y <- DT_h$LVA_v
    if (method == "Shapley") y <- (DT_h[[paste0("Delta_", v, "_mean")]] + DT_h$LVA_v) / 2
    r <- spearman_boot_ci(DT_h$dT_signed, y, B = 1000L,
                            seed = 42L + which(VARS == v))
    rows_sp[[length(rows_sp)+1L]] <- data.table(
      method = method, Variable = v,
      rho_S = r$rho, p_value = r$p,
      ci_lo = r$ci_lo, ci_hi = r$ci_hi
    )
  }
}
tab_sp <- rbindlist(rows_sp)
fwrite(tab_sp, file.path(OUT_TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
cli_alert("HOBO Spearman dTmax (floor05) :"); print(tab_sp)

cli_alert_success("Phase 4 done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')),1)} min")
