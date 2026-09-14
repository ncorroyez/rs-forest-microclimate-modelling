# Compute exact HOBO Shapley on the 41 sensors with full 16-coalition data.
# (12 sensors excluded due to MuSICA failures on NULL/LVA_LAI scenarios.)
suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")
VARS <- c("LAI","Hmax","fCover","LAD")

cli_h1("Exact 16-coalition Shapley per HOBO sensor (n = 41)")
DT <- readRDS(file.path(DATA, "DT_daily_HOBO_floor05_16coalitions.rds"))
DT_agg <- DT[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE)),
              by = .(id_plot, bit_code)]
DT_w <- dcast(DT_agg, id_plot ~ bit_code, value.var = "Tmax_mean")
bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
            "1000","1001","1010","1011","1100","1101","1110","1111")
have_cols <- intersect(bits, names(DT_w))
if (length(have_cols) < 16) {
  cli_alert_warning("Only {length(have_cols)}/16 coalitions present : {paste(setdiff(bits, have_cols), collapse=',')} missing")
}
setnames(DT_w, have_cols, paste0("Tmax_mean_", have_cols))
tcols <- paste0("Tmax_mean_", bits)
tcols_have <- tcols[tcols %in% names(DT_w)]

# Keep sensors with all 16 coalitions present and non-NA
ok_sensors <- DT_w$id_plot
for (s in DT_w$id_plot) {
  row_vals <- as.numeric(unlist(DT_w[id_plot == s, ..tcols_have]))
  if (length(row_vals) < 16 || any(is.na(row_vals))) {
    ok_sensors <- setdiff(ok_sensors, s)
  }
}
cli_alert("Sensors with full 16 coalitions : {length(ok_sensors)} / {nrow(DT_w)}")
DT_w <- DT_w[id_plot %in% ok_sensors]

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
  phi <- shap_one(tvec)
  shap_rows[[i]] <- data.table(id_plot = DT_w$id_plot[i],
                                  variable = names(phi),
                                  phi = unname(phi))
}
DT_shap <- rbindlist(shap_rows)
saveRDS(DT_shap, file.path(DATA, "DT_shapley_per_HOBO_floor05.rds"))
cli_alert_success("Exact Shapley per HOBO sensor : {nrow(DT_shap)/4} sensors")

# Sanity check
chk <- merge(DT_shap[, .(sum_phi = sum(phi)), by = id_plot],
               DT_w[, .(id_plot, v_full = Tmax_mean_1111 - Tmax_mean_0000)],
               by = "id_plot")
chk[, dev := round(sum_phi - v_full, 6)]
cli_alert("Efficiency check : max|sum(phi) - v(full)| = {max(abs(chk$dev))}")

# ---- Spearman + bootstrap CI95 -----------------------------------------------
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
                                  ci_lo = r$ci_lo, ci_hi = r$ci_hi,
                                  n_sensors = nrow(DT_eval))
}
tab_shap_new <- rbindlist(shap_spear)
cli_alert("Exact HOBO Shapley Spearman :"); print(tab_shap_new)

# Update CSV
tab_old <- fread(file.path(TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
tab_old <- tab_old[method != "Shapley"]
tab_updated <- rbind(tab_old, tab_shap_new[, .(method, Variable, rho_S,
                                                  p_value, ci_lo, ci_hi)],
                       fill = TRUE)
fwrite(tab_updated, file.path(TAB, "tab_HOBO_spearman_dTmax_floor05.csv"))
cli_alert_success("Updated tab_HOBO_spearman_dTmax_floor05.csv")

# ---- Regenerate figures -----------------------------------------------------
cli_h2("Regenerate HOBO figures")
DT_h <- readRDS(file.path(DATA, "DT_contrib_HOBO_floor05.rds"))
DT_h <- merge(DT_h, hobo_t, by = "id_plot")
for (v in VARS) {
  DT_h[, paste0("LVA_", v) := get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
  DT_h[, paste0("LOVB_", v) := get(paste0("Delta_", v, "_mean"))]
}
DT_h <- merge(DT_h, shap_w, by = "id_plot", all.x = TRUE)

long_hobo <- rbindlist(lapply(c("LOVB","LVA","Shapley"), function(m)
  rbindlist(lapply(VARS, function(v) data.table(
    method = m, Variable = v,
    dT = DT_h$dT,
    y  = DT_h[[paste0(m, "_", v)]])))))
long_hobo[, method   := factor(method, levels = c("LOVB","LVA","Shapley"))]
long_hobo[, Variable := factor(Variable, levels = VARS)]
long_hobo <- long_hobo[!is.na(y)]

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
    title    = bquote("HOBO sensor-level validation — " ~
                        rho[Spearman] ~ "vs observed" ~ Delta * T[max]),
    subtitle = "Rows : LOVB (n=53) / LVA (n=53) / Shapley exact 16-coalition (n=41)     |     LOESS span = 0.9     |     Bootstrap CI95 (1000 reps)",
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
cli_alert_success("Refreshed fig_HOBO_facet.png (Shapley exact n=41)")

# Standalone Shapley
long_s <- rbindlist(lapply(VARS, function(v) data.table(
  Variable = v, dT = DT_h$dT, y = DT_h[[paste0("Shapley_", v)]])))
long_s[, Variable := factor(Variable, levels = VARS)]
long_s <- long_s[!is.na(y)]
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
    title    = bquote("HOBO validation — Shapley exact (16 coalitions, n=41)"
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
