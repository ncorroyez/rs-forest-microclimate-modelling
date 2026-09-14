# ==============================================================================
# Fill missing HOBO coalition NCs from legacy H1f_* directories.
# Then recompute exact 16-coalition Shapley per sensor (53 sensors).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(ncdf4); library(tidyverse); library(musica.tools)
})
source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")
VARS <- c("LAI","Hmax","fCover","LAD")

cli_h1("Fill missing HOBO coalition NCs from legacy dirs")
DT <- readRDS(file.path(DATA, "DT_daily_HOBO_floor05_16coalitions.rds"))
all_bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
                "1000","1001","1010","1011","1100","1101","1110","1111")

# Find missing (sensor, coalition) pairs
have <- unique(DT[, .(id_plot, bit_code)])
all_sensors <- unique(DT$id_plot)
all_pairs <- CJ(id_plot = all_sensors, bit_code = all_bits)
missing <- fsetdiff(all_pairs, have)
cli_alert("Missing (sensor × coalition) pairs : {nrow(missing)}")
print(missing[, .N, by = bit_code])

# Legacy H1f_* directory map
legacy_map <- c(
  "0000" = "H1f_0_Null_baseline",
  "1000" = "H1f_1_LAI_only",
  "1100" = "H1f_2_LAI_Hmax",
  "1110" = "H1f_3_LAI_Hmax_fCover",
  "1111" = "H1f_4_Full_real"
)
# Also REF_all_real = 1111 alternative
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
hobo_legacy <- here::here("out_files/musica_hobo_validation")
recovered <- list()
for (i in seq_len(nrow(missing))) {
  bit <- missing$bit_code[i]; sensor <- missing$id_plot[i]
  legacy_dir <- legacy_map[bit]
  if (is.na(legacy_dir)) {
    cli_alert_warning("No legacy map for {bit} (sensor {sensor})")
    next
  }
  nc <- file.path(hobo_legacy, legacy_dir,
                    sprintf("musica_out_HOBO_%s.nc", sensor))
  if (!file.exists(nc) || file.size(nc) < 5e4) {
    cli_alert_warning("Missing legacy NC : {basename(nc)}")
    next
  }
  res <- tryCatch(extract_deltatmax_one(nc, df_macro, CFG$date_seq),
                   error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0L) next
  dt <- as.data.table(res)[, .(date, Delta_Tmax)]
  dt[, id_plot := sensor]; dt[, bit_code := bit]
  recovered[[length(recovered)+1L]] <- dt[, .(id_plot, date, bit_code, Delta_Tmax)]
}
if (length(recovered) == 0) stop("No recoveries succeeded")
DT_recovered <- rbindlist(recovered, fill = TRUE)
cat("DT_recovered cols :", names(DT_recovered), "  nrow :", nrow(DT_recovered), "\n")
cli_alert("Recovered : {nrow(DT_recovered)} daily rows ({length(unique(paste(DT_recovered$id_plot, DT_recovered$bit_code)))} sensor×coalition pairs)")

# Merge
DT_full <- rbind(DT, DT_recovered, fill = TRUE)
setkey(DT_full, id_plot, bit_code, date)
saveRDS(DT_full,
         file.path(DATA, "DT_daily_HOBO_floor05_16coalitions.rds"))
cli_alert_success("Updated DT_daily_HOBO_floor05_16coalitions.rds : {nrow(DT_full)} rows")

# Check coverage now
final_pairs <- unique(DT_full[, .(id_plot, bit_code)])
all_pairs2 <- CJ(id_plot = unique(DT_full$id_plot), bit_code = all_bits)
still_missing <- fsetdiff(all_pairs2, final_pairs)
cli_alert("Still missing after legacy fill : {nrow(still_missing)}")
if (nrow(still_missing) > 0) print(still_missing[order(id_plot)])

# ---- Recompute exact Shapley per sensor -------------------------------------
DT_agg <- DT_full[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE)),
                    by = .(id_plot, bit_code)]
DT_w <- dcast(DT_agg, id_plot ~ bit_code, value.var = "Tmax_mean")
bits <- all_bits
setnames(DT_w, bits, paste0("Tmax_mean_", bits))
tcols <- paste0("Tmax_mean_", bits)

shap_one <- function(tvec) {
  v_null <- as.numeric(tvec["0000"])
  v <- function(b) as.numeric(tvec[b]) - v_null
  mb <- function(bo, n = 4L) { b <- rep("0", n); b[bo] <- "1"; paste(b, collapse = "") }
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
shap_rows <- list()
for (i in seq_len(nrow(DT_w))) {
  tvec <- as.numeric(unlist(DT_w[i, ..tcols])); names(tvec) <- bits
  if (any(is.na(tvec))) { cli_alert_warning("Still NAs for {DT_w$id_plot[i]}"); next }
  phi <- shap_one(tvec)
  shap_rows[[i]] <- data.table(id_plot = DT_w$id_plot[i],
                                  variable = names(phi),
                                  phi = unname(phi))
}
DT_shap <- rbindlist(shap_rows)
saveRDS(DT_shap, file.path(DATA, "DT_shapley_per_HOBO_floor05.rds"))
cli_alert_success("Exact per-HOBO Shapley : {nrow(DT_shap)/4} sensors")

# Spearman vs dTmax_obs
cg <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
hobo_t <- as.data.table(cg$HOBO)[, .(id_plot = id, dT = -dTmax_mean)]
shap_w <- dcast(DT_shap, id_plot ~ variable, value.var = "phi")
setnames(shap_w, VARS, paste0("Shapley_", VARS))
DT_eval <- merge(hobo_t, shap_w, by = "id_plot")

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

shap_spear <- list()
for (v in VARS) {
  r <- spearman_boot_ci(DT_eval$dT, DT_eval[[paste0("Shapley_", v)]],
                          B = 1000L, seed = 42L + which(VARS == v))
  shap_spear[[v]] <- data.table(method = "Shapley", Variable = v,
                                  rho_S = r$rho, p_value = r$p,
                                  ci_lo = r$ci_lo, ci_hi = r$ci_hi)
}
tab_shap_new <- rbindlist(shap_spear)
cli_alert("Exact HOBO Shapley Spearman (n = {nrow(DT_eval)}) :"); print(tab_shap_new)

# Update table
tab_old <- fread(file.path(TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
tab_old <- tab_old[method != "Shapley"]
tab_updated <- rbind(tab_old, tab_shap_new, fill = TRUE)
fwrite(tab_updated, file.path(TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
cli_alert_success("Updated tab_HOBO_spearman_dTmax_floor05.csv")

# ---- Regenerate fig_HOBO_facet.png and fig_HOBO_Shapley_spearman_dTmax.png ---
cli_h2("Regenerate HOBO figures with exact Shapley")
DT_h <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
DT_h <- merge(DT_h, hobo_t, by = "id_plot")
for (v in VARS) {
  DT_h[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
  DT_h[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
}
DT_h <- merge(DT_h, shap_w, by = "id_plot")
long_hobo <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m)
  rbindlist(lapply(VARS, function(v) data.table(
    method = m, Variable = v,
    dT = DT_h$dT,
    y  = DT_h[[paste0(m, "_", v)]])))))
long_hobo[, method   := factor(method, levels = c("LOVB","LVA","Shapley"))]
long_hobo[, Variable := factor(Variable, levels = VARS)]
tab_sp <- fread(file.path(TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
tab_sp[, Variable := factor(Variable, levels = VARS)]
tab_sp[, method   := factor(method, levels = c("LOVB","LVA","Shapley"))]
tab_sp[, rho_lab := sprintf("rho == %+.2f", rho_S)]
tab_sp[, p_lab   := ifelse(p_value < 1e-3, "italic(p) < 0.001",
                              sprintf("italic(p) == %.3f", p_value))]
tab_sp[, ci_lab  := sprintf("CI[95] *' '* '[' * %+.2f * '; ' * %+.2f * ']'",
                                ci_lo, ci_hi)]

p_hobo <- ggplot(long_hobo, aes(x = dT, y = y)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                colour = "#FFB400", fill = "#FFB400",
                linewidth = 0.7, alpha = 0.18, span = 0.9) +
  geom_text(data = tab_sp, aes(x = -Inf, y = Inf, label = rho_lab),
             parse = TRUE, hjust = -0.08, vjust = 1.3, size = 4.5,
             fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = tab_sp, aes(x = -Inf, y = Inf, label = p_lab),
             parse = TRUE, hjust = -0.08, vjust = 3.0, size = 3.8,
             inherit.aes = FALSE, colour = "grey25") +
  geom_text(data = tab_sp, aes(x = -Inf, y = Inf, label = ci_lab),
             parse = TRUE, hjust = -0.06, vjust = 5.0, size = 3.3,
             inherit.aes = FALSE, colour = "grey35") +
  facet_grid(method ~ Variable, scales = "free_y", switch = "y") +
  labs(
    title    = bquote("HOBO sensor-level validation (n = 53) — "
                        ~ rho[Spearman] ~ "vs observed" ~ Delta * T[max]),
    subtitle = "Rows : LOVB / LVA / Shapley (exact 16-coalition)     |     LOESS span = 0.9",
    x = bquote(Delta * T[max]^"obs" ~ "(°C)"),
    y = bquote("Simulated  " * Delta[v] ~ "(°C)")
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text       = element_text(face = "bold", size = 13),
         strip.placement  = "outside",
         plot.title       = element_text(face = "bold", size = 16),
         plot.subtitle    = element_text(colour = "grey25", size = 12),
         panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_HOBO_facet.png"),
        p_hobo, width = 16, height = 13, dpi = 300)
cli_alert_success("Refreshed fig_HOBO_facet.png")

# Standalone Shapley HOBO figure
long_s <- rbindlist(lapply(VARS, function(v) data.table(
  Variable = v, dT = DT_h$dT, y = DT_h[[paste0("Shapley_", v)]])))
long_s[, Variable := factor(Variable, levels = VARS)]
ann_s <- tab_sp[method == "Shapley", .(Variable, rho_lab, p_lab, ci_lab)]
ann_s[, Variable := factor(Variable, levels = VARS)]
p_shap <- ggplot(long_s, aes(x = dT, y = y)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(colour = "#2C5F2D", alpha = 0.7, size = 1.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                colour = "#FFB400", fill = "#FFB400",
                linewidth = 0.7, alpha = 0.18, span = 0.9) +
  geom_text(data = ann_s, aes(x = -Inf, y = Inf, label = rho_lab),
             parse = TRUE, hjust = -0.08, vjust = 1.3, size = 4.5,
             fontface = "bold", inherit.aes = FALSE, colour = "grey15") +
  geom_text(data = ann_s, aes(x = -Inf, y = Inf, label = p_lab),
             parse = TRUE, hjust = -0.08, vjust = 3.0, size = 3.8,
             inherit.aes = FALSE, colour = "grey25") +
  geom_text(data = ann_s, aes(x = -Inf, y = Inf, label = ci_lab),
             parse = TRUE, hjust = -0.06, vjust = 5.0, size = 3.3,
             inherit.aes = FALSE, colour = "grey35") +
  facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
  labs(
    title    = bquote("HOBO validation (n = 53) — Shapley exact (16 coalitions)"
                        ~ "   " ~ rho[Spearman] ~ "vs" ~ Delta * T[max]),
    subtitle = "LOESS span = 0.9     |     Spearman ρ + bootstrap CI95 (1000 reps)",
    x = bquote(Delta * T[max]^"obs" ~ "(°C)"),
    y = bquote(Delta[v]^"Shapley" ~ "(°C)")
  ) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
         strip.text       = element_text(face = "bold", size = 12),
         plot.title       = element_text(face = "bold", size = 15),
         plot.subtitle    = element_text(colour = "grey25", size = 11),
         panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_HOBO_Shapley_spearman_dTmax.png"),
        p_shap, width = 14, height = 4.5, dpi = 300)
cli_alert_success("Refreshed fig_HOBO_Shapley_spearman_dTmax.png")
