# ==============================================================================
# Heatmap LOO on hourly log(slope), V11 archetypes (model-3.2.3).
# Mirror of make_heatmap_v10_logslope.R, sources V11 archetype outputs.
# Output : fig_heatmap_LOO_v11_logslope.{png,pdf}
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

sc_name_of <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

arch_nc_path <- function(cl, bit) {
  arch_lbl <- sprintf("Arch_C%s", cl)
  file.path("out_files/H1_archetypes_v11_model323", arch_lbl,
            sprintf("musica_out_%s_%s.nc", arch_lbl, sc_name_of[bit]))
}

extract_logslope <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  res <- tryCatch(extract_hourly_slope_one(path, era5_hourly, CFG$date_seq,
                                            z_target = CFG$tair_target_height,
                                            use_nair1 = TRUE, time_shift_hr = 0),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  log(abs(res$slope))
}

cli_h1("V11 heatmap LOO log(slope)")
rows <- list()
for (cl in 1:4) for (bit in unname(.ARCH_COAL_MAP)) {
  rows[[length(rows) + 1L]] <- data.table(
    Cluster = cl, bit = bit,
    val = extract_logslope(arch_nc_path(cl, bit)))
}
W <- dcast(rbindlist(rows), Cluster ~ bit, value.var = "val")

VARS     <- c("LAI", "Hmax", "fCover", "LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

var_lab_pm <- c(
  LAI    = "atop(LAI, \"(m\"^\"2\"~\"m\"^\"-2\"*\")\")",
  Hmax   = "atop(italic(H)[max], \"(m)\")",
  fCover = "atop(fCover, \"(\\u2013)\")",
  LAD    = "atop(LAD, \"(m\"^\"2\"~\"m\"^\"-3\"*\")\")"
)

delta_long <- list()
for (v in VARS) for (i in seq_len(nrow(W))) {
  delta_long[[length(delta_long) + 1L]] <- data.table(
    Cluster = W$Cluster[i], variable = v,
    delta   = W[["1111"]][i] - W[[loo_bits[v]]][i])
}
long <- rbindlist(delta_long)
long[, Profile := as.character(relabel_cluster(Cluster))]
avg_rows <- long[, .(Profile = "Avg", delta = mean(delta, na.rm = TRUE)),
                     by = variable]
long_full <- rbind(long[, .(variable, Profile, delta)], avg_rows)
long_full[, Profile := factor(Profile, levels = c(paste0("P", 1:4), "Avg"))]

row_order <- c("LAD", "Hmax", "fCover", "LAI")
long_full[, var_lab := var_lab_pm[as.character(variable)]]
long_full[, var_lab := factor(var_lab, levels = var_lab_pm[row_order])]

max_abs <- max(abs(long_full$delta), na.rm = TRUE) * 1.05

p <- ggplot(long_full, aes(x = var_lab, y = Profile, fill = delta)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.3f", delta)),
            fontface = "bold", size = 8,
            colour = ifelse(abs(long_full$delta) > max_abs * 0.55, "white", "grey15")) +
  geom_hline(yintercept = 1.5, colour = "white", linewidth = 4) +
  geom_hline(yintercept = 1.5, colour = "grey40", linewidth = 1.0) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                        midpoint = 0, limits = c(-max_abs, max_abs),
                        name = expression(Delta[v]^"LOO" ~ "(log " * slope * ")")) +
  scale_x_discrete(labels = function(x) parse(text = as.character(x))) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_text(face = "bold", size = 24),
        axis.text.y = element_text(face = "bold", size = 24),
        axis.ticks  = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        legend.key.height = unit(1.6, "cm"),
        legend.position = "right")

ggsave(file.path(OUT, "fig_heatmap_LOO_v11_logslope.png"),
        p, width = 16, height = 9, dpi = 600, bg = "white")
ggsave(file.path(OUT, "fig_heatmap_LOO_v11_logslope.pdf"),
        p, width = 16, height = 9, device = cairo_pdf)
cli_alert_success("Saved fig_heatmap_LOO_v11_logslope (png + pdf)")

fwrite(long_full, file.path(OUT, "tab_heatmap_LOO_v11_logslope.csv"))
