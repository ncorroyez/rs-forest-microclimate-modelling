# ==============================================================================
# Figure 3 — Type B Shapley : per-plot phi_v for cLHS (400 plots).
#
# Requires DT_contrib_cLHS_16coalitions.rds (all 16 coalitions per plot).
# Computes phi_v per plot via exact 2^4 enumeration (24 permutations).
# Value function : v(S) = Tmax_mean(S) − Tmax_mean(NULL).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
})

.VAR_ORDER    <- c("LAI", "fCover", "Hmax", "LAD")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")

# Exact Shapley phi_v from a named vector of 16 coalition values (buffering, °C).
# Convention : v(S) = Tmax_mean(S) - Tmax_mean(NULL), so phi_v has units of °C
# and sum over v gives Tmax_mean(REF) - Tmax_mean(NULL).
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
        total <- total + w * (v(make_bin(c(S, var_bit))) -
                                  v(make_bin(S)))
    }
    out[var_bit] <- total
  }
  setNames(out, variables)
}

compute_per_plot_shapley <- function() {
  cli_h1("Compute per-plot Shapley phi_v (exact, 16 coalitions per plot)")
  DT <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS_16coalitions.rds"))
  bits  <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111",
              "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111")
  tcols <- paste0("Tmax_mean_", bits)
  missing <- setdiff(tcols, names(DT))
  if (length(missing) > 0L)
    stop("Missing coalition columns in DT : ",
          paste(missing, collapse = ", "))
  cli_alert("Computing phi_v for {nrow(DT)} plots ...")

  rows <- vector("list", nrow(DT))
  for (i in seq_len(nrow(DT))) {
    tmax_vec <- as.numeric(unlist(DT[i, ..tcols]))
    names(tmax_vec) <- bits
    if (any(is.na(tmax_vec))) next
    phi <- shap_one_plot(tmax_vec)
    rows[[i]] <- data.table(
      x = DT$x[i], y = DT$y[i], Cluster = DT$Cluster[i],
      LAI = DT$LAI[i], Hmax = DT$Hmax[i],
      fCover = DT$fCover[i], FPC1 = DT$FPC1[i],
      variable = names(phi), phi = unname(phi)
    )
  }
  long <- rbindlist(rows)
  cli_alert_success("phi_v computed for {length(unique(long[, .(x,y)])) } plots")

  saveRDS(long, here::here("outputs/lovb/data/DT_shapley_per_plot.rds"))
  cli_alert_success("Saved DT_shapley_per_plot.rds")
  long
}

build_typeB_Shapley <- function() {
  cli_h1("Figure 3 — Type B Shapley (per plot, mean over 122 d)")
  long <- if (file.exists(here::here("outputs/lovb/data/DT_shapley_per_plot.rds")))
              readRDS(here::here("outputs/lovb/data/DT_shapley_per_plot.rds"))
            else compute_per_plot_shapley()
  trait_map <- c(LAI = "LAI", fCover = "fCover", Hmax = "Hmax", LAD = "FPC1")
  long[, trait := mapply(function(v, i)
                            long[[trait_map[v]]][i], variable, seq_len(.N))]
  # vectorized version
  long[, trait := NA_real_]
  for (v in .VAR_ORDER) {
    long[variable == v, trait := get(trait_map[v])]
  }
  long[, variable    := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]

  cli_alert("Shapley Type B : {nrow(long) / length(.VAR_ORDER)} plots per panel")

  p <- ggplot(long, aes(x = trait, y = phi)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
                linewidth = 0.6) +
    geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.6) +
    geom_smooth(method = "loess", span = 0.75, se = TRUE,
                  colour = "grey20", fill = "grey70", alpha = 0.3) +
    facet_wrap(~ variable, nrow = 1, scales = "free_x",
                labeller = labeller(variable = function(x)
                                       paste("Contribution of", x))) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         name = NULL) +
    labs(
      title    = "Contribution phi_v vs trait (MEAN) — Shapley",
      subtitle = "Exact phi_v per plot from 16 coalitions   |   x = real value of v (FPC1 for LAD)   |   LOESS span=0.75 with 95%CI",
      x = "Trait value (real per plot)",
      y = "phi_v   Tmax (mean)  [degree C]"
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
           plot.title      = element_text(face = "bold"),
           legend.position = "bottom")

  fig_path <- here::here("outputs/lovb/figures/fig_typeB_Shapley_mean.png")
  ggsave(fig_path, p, width = 14, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  step <- if (length(args) >= 1L) args[[1L]] else "compute"
  if (step == "compute") compute_per_plot_shapley()
  else if (step == "figure")  build_typeB_Shapley()
  else if (step == "all") {
    compute_per_plot_shapley()
    build_typeB_Shapley()
  } else stop("Unknown : ", step)
}
