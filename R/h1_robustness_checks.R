# ==============================================================================
# Robustness checks for the LAI > Hmax ≈ fCover > LAD hierarchy.
#
# Check 1 : temporal aggregation     -- P90 of daily Tmax instead of mean
# Check 2 : magnitude metric         -- RMSE and R² instead of MAE
# Check 3 : baseline / aggregator    -- Q10 and Q90 across plots instead of mean
#
# All three checks operate on the floor05 caches (cLHS scale, n = 400 plots,
# 16 coalitions, n=4 archetypes).  Exact Shapley is recomputed from the full
# 16-coalition table for each aggregation level (mean/P90).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
})

DATA <- here::here("outputs/lovb_floor05/data")
OUT  <- here::here("outputs/robustness")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

VARS    <- c("LAI","Hmax","fCover","LAD")
BIT_POS <- setNames(seq_along(VARS), VARS)   # LAI=1, Hmax=2, fCover=3, LAD=4

# 16 bit codes, ordered 0000 .. 1111
all_bits <- sprintf("%04d", as.integer(R.utils::intToBin(0:15)))
shapley_weights <- function(n) {
  # weight w(|S|) = |S|! (n-|S|-1)! / n!
  sapply(0:(n-1), function(k) factorial(k) * factorial(n - k - 1) / factorial(n))
}

# ----------------------------------------------------------------------------
# Generic Shapley : returns phi_v per plot for a given aggregation Tmax_<agg>_XXXX
# Convention C : phi_v > 0 means v contributes to buffering (reduces Tmax)
# ----------------------------------------------------------------------------
compute_shapley <- function(DT16, agg = "mean") {
  n <- length(VARS); w <- shapley_weights(n)
  cols <- paste0("Tmax_", agg, "_", all_bits)
  M <- as.matrix(DT16[, ..cols])
  colnames(M) <- all_bits
  out <- data.table(x = DT16$x, y = DT16$y)
  for (v in VARS) {
    pos <- BIT_POS[v]
    phi_plot <- numeric(nrow(M))
    for (i in 0:15) {
      bc <- sprintf("%04d", as.integer(R.utils::intToBin(i)))
      bits <- as.integer(strsplit(bc, "")[[1]])
      if (bits[pos] == 1) next
      bits_with <- bits; bits_with[pos] <- 1L
      bc_with <- paste(bits_with, collapse = "")
      k <- sum(bits)
      # buffering convention : Tmax(S) - Tmax(S ∪ {v})
      phi_plot <- phi_plot + w[k + 1] * (M[, bc] - M[, bc_with])
    }
    out[[paste0("phi_", v)]] <- phi_plot
  }
  out
}

# ----------------------------------------------------------------------------
# Per-plot Δ_v for LOO (= LOVB) and LVA, in convention C (positive = buffering)
# LOO : Δ_v = Tmax_LOVB_v - Tmax_REF  (removing v warmed up the canopy)
# LVA : Δ_v = Tmax_LVA_v  - Tmax_NULL (adding v alone cooled the canopy)
# ----------------------------------------------------------------------------
compute_attrib_C <- function(DT16, agg = "mean") {
  cols <- paste0("Tmax_", agg, "_", all_bits)
  M <- as.matrix(DT16[, ..cols]); colnames(M) <- all_bits
  out <- data.table(x = DT16$x, y = DT16$y)
  for (v in VARS) {
    pos <- BIT_POS[v]
    bc_lovb <- rep("1", 4); bc_lovb[pos] <- "0"; bc_lovb <- paste(bc_lovb, collapse="")
    bc_lva  <- rep("0", 4); bc_lva[pos]  <- "1"; bc_lva  <- paste(bc_lva,  collapse="")
    out[[paste0("LOO_", v)]] <- M[, bc_lovb] - M[, "1111"]
    out[[paste0("LVA_", v)]] <- M[, bc_lva]  - M[, "0000"]
  }
  out
}

# ----------------------------------------------------------------------------
# Ranking helper : given per-plot attributions, return per-variable scalar score
# under one of several aggregators
# ----------------------------------------------------------------------------
scalar_score <- function(x, method = c("MAE","RMSE","Q10","Q90","Median")) {
  method <- match.arg(method)
  x <- x[!is.na(x)]
  switch(method,
    MAE    = mean(abs(x)),
    RMSE   = sqrt(mean(x^2)),
    Q10    = quantile(abs(x), 0.10, names = FALSE),
    Q90    = quantile(abs(x), 0.90, names = FALSE),
    Median = median(abs(x)))
}

rank_by <- function(scores) {
  # input : named numeric vector of per-variable scores ; output : rank 1..4
  # higher score = more contribution = rank 1
  ranks <- rank(-scores, ties.method = "min")
  ranks[VARS]
}

# ----------------------------------------------------------------------------
# Master compute : returns long table (method × variable × aggregator) of
# scalar score and rank
# ----------------------------------------------------------------------------
build_rank_table <- function(DT16, agg = "mean",
                                aggregators = c("MAE","RMSE","Q10","Q90")) {
  att <- compute_attrib_C(DT16, agg = agg)
  shap <- compute_shapley(DT16, agg = agg)
  rows <- list()
  for (m in c("LOO","LVA","Shapley")) {
    for (aggr in aggregators) {
      scores <- sapply(VARS, function(v) {
        col <- if (m == "Shapley") paste0("phi_", v) else paste0(m, "_", v)
        scalar_score(if (m == "Shapley") shap[[col]] else att[[col]], aggr)
      })
      names(scores) <- VARS
      rows[[length(rows) + 1L]] <- data.table(
        Temporal_agg = agg, Aggregator = aggr, Method = m,
        Variable = VARS, Score = scores[VARS], Rank = rank_by(scores))
    }
  }
  rbindlist(rows)
}

# ============================================================================
# R² check : for each (method, v), how well does the masked simulation
#             reconstruct REF.  Lower R² → stronger contribution.
# ============================================================================
r2_check <- function(DT16, agg = "mean") {
  cols <- paste0("Tmax_", agg, "_", all_bits)
  M <- as.matrix(DT16[, ..cols]); colnames(M) <- all_bits
  ref <- M[, "1111"]; null_ <- M[, "0000"]
  rows <- list()
  ss_tot <- sum((ref - mean(ref))^2)
  for (v in VARS) {
    pos <- BIT_POS[v]
    bc_lovb <- rep("1", 4); bc_lovb[pos] <- "0"; bc_lovb <- paste(bc_lovb, collapse="")
    bc_lva  <- rep("0", 4); bc_lva[pos]  <- "1"; bc_lva  <- paste(bc_lva,  collapse="")
    r2_lovb <- 1 - sum((M[, bc_lovb] - ref)^2) / ss_tot
    r2_lva  <- 1 - sum((M[, bc_lva]  - ref)^2) / ss_tot
    rows[[length(rows) + 1L]] <- data.table(
      Variable = v,
      R2_LOO_vs_REF = r2_lovb,
      R2_LVA_vs_REF = r2_lva,
      RMSE_LOO_vs_REF = sqrt(mean((M[, bc_lovb] - ref)^2)),
      RMSE_LVA_vs_REF = sqrt(mean((M[, bc_lva]  - ref)^2)))
  }
  rbindlist(rows)
}

# ============================================================================
# Driver
# ============================================================================
main <- function() {
  if (!requireNamespace("R.utils", quietly = TRUE))
    stop("install.packages('R.utils') first")
  cli_h1("Robustness checks — cLHS scale (n=400)")
  DT16 <- readRDS(file.path(DATA, "DT_contrib_cLHS_floor05_16coalitions.rds"))

  ## ---- Check 1 : temporal aggregation P90 vs mean --------------------------
  tab_mean <- build_rank_table(DT16, "mean")
  tab_P90  <- build_rank_table(DT16, "P90")
  tab_both <- rbind(tab_mean, tab_P90)
  fwrite(tab_both, file.path(OUT, "tab_rank_temporal_aggregator.csv"))

  cli_h2("Check 1 : ranks under mean vs P90")
  print(dcast(tab_both[Aggregator == "MAE"],
                Method + Variable ~ Temporal_agg, value.var = "Rank"))

  ## ---- Check 2 : magnitude metric ------------------------------------------
  cli_h2("Check 2 : MAE vs RMSE ranks (mean Tmax)")
  print(dcast(tab_mean[Aggregator %in% c("MAE","RMSE")],
                Method + Variable ~ Aggregator, value.var = "Rank"))

  cli_h2("Check 2-bis : R² of masked simulation vs REF (lower = bigger effect)")
  r2_tab <- r2_check(DT16, "mean")
  r2_tab[, RMSE_LOO_rank := rank(-RMSE_LOO_vs_REF)]
  r2_tab[, RMSE_LVA_rank := rank(-RMSE_LVA_vs_REF)]
  r2_tab[, R2_LOO_rank   := rank(R2_LOO_vs_REF)]   # lower R² → rank 1
  r2_tab[, R2_LVA_rank   := rank(R2_LVA_vs_REF)]
  fwrite(r2_tab, file.path(OUT, "tab_r2_rmse_check.csv"))
  print(r2_tab)

  ## ---- Check 3 : Q10 / Q90 cross-plot aggregator ---------------------------
  cli_h2("Check 3 : MAE vs Q10 vs Q90 ranks")
  print(dcast(tab_mean[Aggregator %in% c("MAE","Q10","Q90")],
                Method + Variable ~ Aggregator, value.var = "Rank"))

  ## ---- Save a wide summary -------------------------------------------------
  wide <- dcast(tab_mean, Method + Variable ~ Aggregator, value.var = "Rank")
  wide_P90 <- dcast(tab_P90[Aggregator == "MAE"],
                       Method + Variable ~ Temporal_agg, value.var = "Rank")
  setnames(wide_P90, "P90", "MAE_P90")
  summary_wide <- merge(wide, wide_P90, by = c("Method","Variable"))
  fwrite(summary_wide, file.path(OUT, "tab_robustness_summary.csv"))

  cli_h2("FINAL SUMMARY (rank per Method × Variable across robustness axes)")
  print(summary_wide)

  cli_alert_success("Robustness tables saved to {.path {OUT}}")
}

if (sys.nframe() == 0) main()
