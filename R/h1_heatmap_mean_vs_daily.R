# ==============================================================================
# Generate two versions of the heatmap : MEAN (current) vs DAILY.
#   MEAN : attribution computed from Tmax_mean (already averaged over 122 days).
#   DAILY : attribution computed from daily Tmax values, then aggregated.
# Both at the 2 levels (archetypes + cLHS) + consensus column.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")

VARS <- c("LAI","fCover","Hmax","LAD")

# ============================================================================
# Bit-code conventions
# ============================================================================
bit_NULL <- "0000"; bit_REF  <- "1111"
bit_LOVB <- c(LAI="0111", Hmax="1011", fCover="1101", LAD="1110")
bit_LVA  <- c(LAI="1000", Hmax="0100", fCover="0010", LAD="0001")

# ============================================================================
# Shapley exact (16 coalitions) helper
# ============================================================================
all_bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
                "1000","1001","1010","1011","1100","1101","1110","1111")
shap_one <- function(tvec) {
  v_null <- as.numeric(tvec["0000"])
  v <- function(b) as.numeric(tvec[b]) - v_null
  mb <- function(bo, n=4L) { b <- rep("0", n); b[bo] <- "1"; paste(b, collapse="") }
  out <- numeric(4L)
  for (vb in 1:4) {
    others <- setdiff(1:4, vb); total <- 0
    for (s in 0:3) {
      ss <- if (s == 0L) list(integer(0)) else combn(others, s, simplify = FALSE)
      w <- factorial(s) * factorial(4-s-1) / factorial(4)
      for (S in ss) total <- total + w * (v(mb(c(S, vb))) - v(mb(S)))
    }
    out[vb] <- total
  }
  setNames(out, c("LAI","Hmax","fCover","LAD"))
}

# ============================================================================
# Build attribution tables — DAILY version
# ============================================================================
build_daily_attribution <- function() {
  cli_h2("DAILY attribution computation")

  # cLHS daily : DT_daily_cLHS_floor05 has 16 coalitions × 400 plots × 122 days
  DT_d <- readRDS(file.path(DATA, "DT_daily_cLHS_floor05.rds"))

  # ---- LOVB cLHS daily ----
  ref <- DT_d[bit_code == bit_REF, .(x, y, date, Tmax_REF = Delta_Tmax)]
  lovb_daily <- list()
  for (v in VARS) {
    lv <- DT_d[bit_code == bit_LOVB[v], .(x, y, date, Tmax_LOVB = Delta_Tmax)]
    m <- merge(ref, lv, by = c("x","y","date"))
    m[, delta := Tmax_REF - Tmax_LOVB]
    lovb_daily[[v]] <- data.table(
      variable = v,
      LOVB_cLHS_value = sqrt(mean(m$delta^2, na.rm = TRUE)),
      n = nrow(m)
    )
  }
  lovb_clhs_daily <- rbindlist(lovb_daily)
  lovb_clhs_daily[, LOVB_cLHS_rank := frank(-LOVB_cLHS_value, ties.method = "min")]

  # ---- LVA cLHS daily ----
  null_d <- DT_d[bit_code == bit_NULL, .(x, y, date, Tmax_NULL = Delta_Tmax)]
  lva_daily <- list()
  for (v in VARS) {
    lvav <- DT_d[bit_code == bit_LVA[v], .(x, y, date, Tmax_LVA = Delta_Tmax)]
    m <- merge(null_d, lvav, by = c("x","y","date"))
    m[, lva := Tmax_LVA - Tmax_NULL]
    lva_daily[[v]] <- data.table(
      variable = v,
      LVA_cLHS_value = mean(abs(m$lva), na.rm = TRUE),
      n = nrow(m)
    )
  }
  lva_clhs_daily <- rbindlist(lva_daily)
  lva_clhs_daily[, LVA_cLHS_rank := frank(-LVA_cLHS_value, ties.method = "min")]

  # ---- Shapley cLHS daily : compute phi per (plot, day) -----------------------
  # Pivot daily DT to wide (one row per plot×day, 16 coalition columns)
  cli_alert("Computing per-day Shapley (400 plots × 122 days = 48,800 phi computations)...")
  DT_w_d <- dcast(DT_d, x + y + date ~ bit_code, value.var = "Delta_Tmax")
  tcols <- intersect(all_bits, names(DT_w_d))
  if (length(tcols) < 16) cli_alert_warning("Only {length(tcols)}/16 coalitions in daily cache")
  setnames(DT_w_d, tcols, paste0("Tmax_", tcols))
  Tmax_cols <- paste0("Tmax_", all_bits)

  # phi per (plot, day)
  phi_list <- list()
  for (i in seq_len(nrow(DT_w_d))) {
    tvec <- as.numeric(unlist(DT_w_d[i, ..Tmax_cols])); names(tvec) <- all_bits
    if (any(is.na(tvec))) next
    phi <- shap_one(tvec)
    phi_list[[i]] <- data.table(x = DT_w_d$x[i], y = DT_w_d$y[i],
                                  date = DT_w_d$date[i],
                                  variable = names(phi), phi = unname(phi))
  }
  DT_phi_d <- rbindlist(phi_list)
  cli_alert("Computed phi for {length(unique(paste(DT_phi_d$x, DT_phi_d$y, DT_phi_d$date)))} (plot,day) pairs")

  shap_daily <- DT_phi_d[, .(Shapley_cLHS_value = mean(abs(phi), na.rm = TRUE)),
                            by = variable]
  shap_daily[, Shapley_cLHS_rank := frank(-Shapley_cLHS_value, ties.method = "min")]

  # ---- Combine cLHS daily ----
  tab_clhs <- merge(merge(shap_daily, lovb_clhs_daily[, .(variable, LOVB_cLHS_value,
                                                              LOVB_cLHS_rank)],
                            by = "variable"),
                      lva_clhs_daily[, .(variable, LVA_cLHS_value, LVA_cLHS_rank)],
                      by = "variable")
  tab_clhs[, variable := factor(variable, levels = VARS)]
  setorder(tab_clhs, variable)
  fwrite(tab_clhs, file.path(TAB, "tab_attribution_cLHS_daily.csv"))
  cli_alert("cLHS daily attribution :"); print(tab_clhs)

  # ---- Archetype daily ----
  DT_da <- readRDS(file.path(DATA, "DT_daily_archetypes_floor05.rds"))
  # If missing coalitions, load them (as we did before)
  if (length(unique(DT_da$bit_code)) < 16) {
    cli_alert_warning("Archetype daily has only {length(unique(DT_da$bit_code))} coalitions — using 10 available")
  }

  ref_a <- DT_da[bit_code == bit_REF, .(archetype, date, Tmax_REF = Delta_Tmax)]
  lovb_arch_daily <- list()
  for (v in VARS) {
    lv <- DT_da[bit_code == bit_LOVB[v], .(archetype, date, Tmax_LOVB = Delta_Tmax)]
    m <- merge(ref_a, lv, by = c("archetype","date"))
    m[, delta := Tmax_REF - Tmax_LOVB]
    lovb_arch_daily[[v]] <- data.table(
      variable = v,
      LOVB_arch_value = mean(abs(m$delta), na.rm = TRUE)
    )
  }
  lovb_a_daily <- rbindlist(lovb_arch_daily)
  lovb_a_daily[, LOVB_arch_rank := frank(-LOVB_arch_value, ties.method = "min")]

  null_a <- DT_da[bit_code == bit_NULL, .(archetype, date, Tmax_NULL = Delta_Tmax)]
  lva_arch_daily <- list()
  for (v in VARS) {
    lvav <- DT_da[bit_code == bit_LVA[v], .(archetype, date, Tmax_LVA = Delta_Tmax)]
    m <- merge(null_a, lvav, by = c("archetype","date"))
    m[, lva := Tmax_LVA - Tmax_NULL]
    lva_arch_daily[[v]] <- data.table(
      variable = v,
      LVA_arch_value = mean(abs(m$lva), na.rm = TRUE)
    )
  }
  lva_a_daily <- rbindlist(lva_arch_daily)
  lva_a_daily[, LVA_arch_rank := frank(-LVA_arch_value, ties.method = "min")]

  # Shapley arch daily (need full 16 coalitions — they exist for C1, partial for C2-4)
  # Use what's available — pivot wide
  DT_wa <- dcast(DT_da, archetype + date ~ bit_code, value.var = "Delta_Tmax")
  tcols_a <- intersect(all_bits, names(DT_wa))
  setnames(DT_wa, tcols_a, paste0("Tmax_", tcols_a))
  if (length(tcols_a) < 16) {
    cli_alert_warning("Skipping arch Shapley daily — only {length(tcols_a)}/16 coalitions. Using NA.")
    shap_a_daily <- data.table(variable = VARS,
                                  Shapley_arch_value = NA_real_,
                                  Shapley_arch_rank = NA_integer_)
  } else {
    Tmax_cols_a <- paste0("Tmax_", all_bits)
    phi_list_a <- list()
    for (i in seq_len(nrow(DT_wa))) {
      tvec <- as.numeric(unlist(DT_wa[i, ..Tmax_cols_a])); names(tvec) <- all_bits
      if (any(is.na(tvec))) next
      phi <- shap_one(tvec)
      phi_list_a[[i]] <- data.table(archetype = DT_wa$archetype[i],
                                      date = DT_wa$date[i],
                                      variable = names(phi), phi = unname(phi))
    }
    DT_phi_a <- rbindlist(phi_list_a)
    shap_a_daily <- DT_phi_a[, .(Shapley_arch_value = mean(abs(phi), na.rm = TRUE)),
                                by = variable]
    shap_a_daily[, Shapley_arch_rank := frank(-Shapley_arch_value, ties.method = "min")]
  }

  tab_arch <- merge(merge(shap_a_daily, lovb_a_daily[, .(variable, LOVB_arch_value,
                                                              LOVB_arch_rank)],
                            by = "variable"),
                      lva_a_daily[, .(variable, LVA_arch_value, LVA_arch_rank)],
                      by = "variable")
  tab_arch[, variable := factor(variable, levels = VARS)]
  setorder(tab_arch, variable)
  fwrite(tab_arch, file.path(TAB, "tab_attribution_archetypes_pooled_daily.csv"))
  cli_alert("Archetype daily attribution :"); print(tab_arch)

  list(clhs = tab_clhs, arch = tab_arch)
}

# ============================================================================
# Build heatmap from attribution tables (works for both MEAN and DAILY)
# ============================================================================
build_heatmap <- function(arch, clhs, title_text, subtitle_text, out_path) {
  arch[, variable := factor(variable, levels = VARS)]
  clhs[, variable := factor(variable, levels = VARS)]
  methods <- c("LOVB_arch","LVA_arch","Shapley_arch",
                "LOVB_cLHS","LVA_cLHS","Shapley_cLHS")
  labels  <- c("LOVB","LVA","Shapley","LOVB","LVA","Shapley")
  long_list <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    m <- methods[i]
    src <- if (grepl("_arch$", m)) arch else clhs
    rank_col <- paste0(m, "_rank")
    long_list[[i]] <- src[, .(variable, method = labels[i], col_idx = i,
                                Rank = get(rank_col))]
  }
  long <- rbindlist(long_list)
  cons <- long[, .(Rank = mean(Rank, na.rm = TRUE), col_idx = 7L), by = variable]
  cons[, method := "Consensus"]
  long_full <- rbind(long[, .(variable, method, col_idx, Rank)],
                      cons[, .(variable, method, col_idx, Rank)])
  long_full[, fill_rank := as.character(round(Rank))]
  pal <- c("1"="#440154","2"="#3B528B","3"="#5DC863","4"="#FDE725")
  x_labels <- c("LOVB","LVA","Shapley","LOVB","LVA","Shapley","Consensus")
  top_labels_df <- data.table(
    x = c(2, 5, 7),
    label = c("Level 1 — archetypes (4 K-means centroids)",
                "Level 2 — cLHS (400 plots)",
                "Consensus"))

  p <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.0) +
    geom_text(aes(label = ifelse(method == "Consensus",
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(long_full$fill_rank) <= 2,
                                 "white", "grey15"),
               fontface = "bold", size = 5.5) +
    geom_vline(xintercept = c(3.5, 6.5), colour = "white", linewidth = 3) +
    geom_text(data = top_labels_df,
               aes(x = x, y = 4.85, label = label),
               inherit.aes = FALSE, fontface = "bold", size = 3.3,
               colour = "grey20") +
    scale_fill_manual(values = pal, name = "Rank", drop = FALSE) +
    scale_x_continuous(breaks = 1:7, labels = x_labels, expand = c(0,0)) +
    scale_y_discrete(limits = rev, expand = c(0,0)) +
    coord_cartesian(ylim = c(0.5, 5.1), clip = "off") +
    labs(
      title    = title_text,
      subtitle = subtitle_text,
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 14),
           plot.subtitle = element_text(colour = "grey25", size = 11),
           axis.text.x   = element_text(face = "bold", size = 10),
           axis.text.y   = element_text(face = "bold", size = 13),
           plot.margin   = margin(15, 10, 10, 10),
           panel.grid    = element_blank())
  ggsave(out_path, p, width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {out_path}}")
}

# ============================================================================
# Main : generate both versions
# ============================================================================
if (sys.nframe() == 0) {
  # MEAN version (rename existing fig_heatmap_2levels.png to _mean suffix)
  cli_h1("MEAN version — copy current heatmap as fig_heatmap_2levels_mean.png")
  file.copy(file.path(FIG, "fig_heatmap_2levels.png"),
              file.path(FIG, "fig_heatmap_2levels_mean.png"),
              overwrite = TRUE)
  cli_alert_success("Copied to fig_heatmap_2levels_mean.png")

  # DAILY version — compute fresh
  cli_h1("DAILY version — compute attribution per day, aggregate ranks")
  out_daily <- build_daily_attribution()
  build_heatmap(out_daily$arch, out_daily$clhs,
                  title_text    = "Attribution ranking — 2 aggregation levels × 3 methods (DAILY)",
                  subtitle_text = "Attribution computed daily on Δ_Tmax values, then aggregated as mean(|Δ_v|) over (plot × day) pairs",
                  out_path      = file.path(FIG, "fig_heatmap_2levels_daily.png"))

  cli_alert_success("Both heatmaps in {.path {FIG}}")
}
