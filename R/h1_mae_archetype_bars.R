# Archetype barplots in |Δ_v| (MAE-compatible absolute values).
suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})
DATA <- here::here("outputs/lovb_floor05/data")
TAB  <- here::here("outputs/lovb_floor05/tables")
FIG  <- here::here("outputs/lovb_floor05/figures")

VARS <- c("LAI","fCover","Hmax","LAD")
.PAL_VAR <- c(LAI = "#440154", fCover = "#31688e",
                Hmax = "#35b779", LAD = "#fde725")
.CLUSTER_ARCH <- c("Arch_C1" = "C1 Open", "Arch_C2" = "C2 Dense",
                    "Arch_C3" = "C3 Low cover", "Arch_C4" = "C4 Inter")

method_yexpr <- function(method) switch(method,
  "LOVB"    = bquote("|" * Delta[v]^"LOVB" * "|" ~ "(°C)"),
  "LVA"     = bquote("|" * Delta[v]^"LVA"  * "|" ~ "(°C)"),
  "Shapley" = bquote("|" * Delta[v]^"Shapley" * "|" ~ "(°C)"))

# Load data
DT_a <- readRDS(file.path(DATA, "DT_contrib_archetypes_floor05.rds"))

# Per-archetype Shapley : need 16 coalitions (load if needed)
DT_daily_a <- readRDS(file.path(DATA, "DT_daily_archetypes_floor05.rds"))
n_per_arch <- DT_daily_a[, length(unique(bit_code)), by = archetype]
if (min(n_per_arch$V1) < 16) {
  suppressMessages({library(ncdf4); library(tidyverse); library(musica.tools)})
  source(here::here("R/config.R")); source(here::here("R/musica.R"))
  source(here::here("R/io.R")); source(here::here("R/lovb_01_load.R"))
  MISSING <- c("0011","0101","0110","1001","1010","1100")
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  extra <- list()
  for (arch in c("Arch_C2","Arch_C3","Arch_C4")) {
    for (bit in MISSING) {
      nc <- here::here("out_files/H1_archetypes", arch,
                         sprintf("musica_out_%s_ARCH_%s.nc", arch, lovb_bit_to_suffix(bit)))
      if (!file.exists(nc) || file.size(nc) < 1e5) next
      res <- tryCatch(extract_deltatmax_one(nc, df_macro, CFG$date_seq), error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0L) next
      dt <- as.data.table(res); dt[, archetype := arch]; dt[, bit_code := bit]
      extra[[length(extra)+1L]] <- dt[, .(archetype, date, bit_code, Delta_Tmax)]
    }
  }
  DT_daily_a <- rbind(DT_daily_a, rbindlist(extra, fill = TRUE), fill = TRUE)
}
DT_agg_a <- DT_daily_a[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE)),
                          by = .(archetype, bit_code)]
all_bits <- c("0000","0001","0010","0011","0100","0101","0110","0111",
                "1000","1001","1010","1011","1100","1101","1110","1111")
DT_w <- dcast(DT_agg_a, archetype ~ bit_code, value.var = "Tmax_mean")
setnames(DT_w, all_bits, paste0("T_", all_bits))
tcols <- paste0("T_", all_bits)
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
shap_rows <- list()
for (i in seq_len(nrow(DT_w))) {
  tvec <- as.numeric(unlist(DT_w[i, ..tcols])); names(tvec) <- all_bits
  if (any(is.na(tvec))) next
  phi <- shap_one(tvec)
  shap_rows[[i]] <- data.table(archetype = DT_w$archetype[i],
                                  variable = names(phi), value = unname(phi))
}
DT_shap_arch <- rbindlist(shap_rows)

make_bars <- function(long_df, method_name, ylab_expr) {
  # Absolute values
  long_df[, abs_val := abs(value)]
  long_df[, variable := factor(variable, levels = VARS)]
  long_df[, archetype_lab := .CLUSTER_ARCH[archetype]]
  long_df[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_ARCH)]
  ymax <- max(long_df$abs_val, na.rm = TRUE) * 1.15
  ggplot(long_df, aes(x = variable, y = abs_val, fill = variable)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_col(colour = "grey20", width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", abs_val), vjust = -0.3),
               fontface = "bold", size = 3.6) +
    facet_wrap(~ archetype_lab, nrow = 1L) +
    scale_fill_manual(values = .PAL_VAR, guide = "none") +
    scale_y_continuous(limits = c(0, ymax),
                        expand = expansion(mult = c(0, 0.08))) +
    labs(title = paste(method_name, "contribution by archetype  (MAE = |Δ_v|)"),
          x = NULL, y = ylab_expr) +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill="grey92", colour=NA),
           strip.text = element_text(face="bold", size=12),
           plot.title = element_text(face="bold", size=14))
}

# LOVB
lovb_long <- rbindlist(lapply(VARS, function(v) data.table(
  archetype = DT_a$archetype, variable = v,
  value = DT_a[[paste0("Delta_", v, "_mean")]])))
ggsave(file.path(FIG, "fig_archetypes_LOVB_mean.png"),
        make_bars(lovb_long, "LOVB", method_yexpr("LOVB")),
        width = 12, height = 5.5, dpi = 300)
cli_alert_success("Saved fig_archetypes_LOVB_mean.png (MAE)")

# LVA
lva_long <- rbindlist(lapply(VARS, function(v) data.table(
  archetype = DT_a$archetype, variable = v,
  value = DT_a[[paste0("LVA_", v, "_mean")]])))
ggsave(file.path(FIG, "fig_archetypes_LVA_mean.png"),
        make_bars(lva_long, "LVA", method_yexpr("LVA")),
        width = 12, height = 5.5, dpi = 300)
cli_alert_success("Saved fig_archetypes_LVA_mean.png (MAE)")

# Shapley
ggsave(file.path(FIG, "fig_archetypes_Shapley_mean.png"),
        make_bars(DT_shap_arch, "Shapley", method_yexpr("Shapley")),
        width = 12, height = 5.5, dpi = 300)
cli_alert_success("Saved fig_archetypes_Shapley_mean.png (MAE)")
