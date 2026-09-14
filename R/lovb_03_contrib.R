# ==============================================================================
# LOVB analysis — Module 3 : contributions (REF vs LOVB, REF vs LVA, NULL ref)
#
# Pivots DT_agg from long (bit_code) to wide and computes :
#   - Delta_v   = Tmax_<agg>_REF - Tmax_<agg>_LOVB_v   (max-order contribution)
#   - LVA_v     = Tmax_<agg>_LVA_v - Tmax_<agg>_NULL    (zero-order contribution)
# for v in {LAI, Hmax, fCover, LAD} and <agg> in {mean, P90}.
#
# Bit-code -> variable name mapping :
#   0000 -> NULL          (LVA reference)
#   0001 -> LVA_LAD       1000 -> LVA_LAI
#   0010 -> LVA_fCover    0111 -> LOVB_LAI
#   0100 -> LVA_Hmax      1011 -> LOVB_Hmax
#   1101 -> LOVB_fCover   1110 -> LOVB_LAD
#   1111 -> REF
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
})

# Coalition bit-code -> readable name suffix (used in column names)
.LOVB_BIT_LABEL <- c(
  "0000" = "NULL",
  "0001" = "LVA_LAD",     "0010" = "LVA_fCover",
  "0100" = "LVA_Hmax",    "1000" = "LVA_LAI",
  "0111" = "LOVB_LAI",    "1011" = "LOVB_Hmax",
  "1101" = "LOVB_fCover", "1110" = "LOVB_LAD",
  "1111" = "REF"
)

# Generic pivot + derived metrics
.lovb_contrib_core <- function(DT_agg, id_cols) {
  stopifnot(all(c("bit_code", "Tmax_mean", "Tmax_P90") %in% names(DT_agg)))
  stopifnot(all(id_cols %in% names(DT_agg)))

  # Pivot wide on bit_code for both metrics
  fmla <- as.formula(paste(paste(id_cols, collapse = " + "), "~ bit_code"))
  DT_w <- dcast(DT_agg, fmla, value.var = c("Tmax_mean", "Tmax_P90"))

  # Rename columns Tmax_mean_<bit> -> Tmax_mean_<label>
  bits <- names(.LOVB_BIT_LABEL)
  old_cols <- c(paste0("Tmax_mean_", bits), paste0("Tmax_P90_", bits))
  new_cols <- c(paste0("Tmax_mean_", .LOVB_BIT_LABEL),
                paste0("Tmax_P90_",  .LOVB_BIT_LABEL))
  have <- intersect(old_cols, names(DT_w))
  ok   <- old_cols %in% have
  setnames(DT_w, old_cols[ok], new_cols[ok])

  # Derived metrics : Delta_v (REF - LOVB_v) and LVA_v (LVA_v - NULL)
  vars <- c("LAI", "Hmax", "fCover", "LAD")
  for (v in vars) {
    DT_w[, paste0("Delta_", v, "_mean") :=
            get("Tmax_mean_REF") - get(paste0("Tmax_mean_LOVB_", v))]
    DT_w[, paste0("Delta_", v, "_P90") :=
            get("Tmax_P90_REF") - get(paste0("Tmax_P90_LOVB_", v))]
    DT_w[, paste0("LVA_", v, "_mean") :=
            get(paste0("Tmax_mean_LVA_", v)) - get("Tmax_mean_NULL")]
    DT_w[, paste0("LVA_", v, "_P90") :=
            get(paste0("Tmax_P90_LVA_", v))  - get("Tmax_P90_NULL")]
  }
  DT_w
}

# ---- 1. cLHS contributions ---------------------------------------------------
lovb_contrib_clhs <- function(DT_agg,
                                cache_path = here("outputs/lovb/data/DT_contrib_cLHS.rds"),
                                use_cache  = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    cli_alert_info("Loading cached cLHS DT_contrib from {.path {cache_path}}")
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Build cLHS DT_contrib (wide + derived metrics)")

  id_cols <- c("x", "y", "Cluster", "LAI", "Hmax", "fCover", "FPC1")
  DT_c <- .lovb_contrib_core(DT_agg, id_cols)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT_c, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT_c)}} rows, {.val {ncol(DT_c)}} cols)")
  DT_c
}

# ---- 2. Archetypes contributions ---------------------------------------------
lovb_contrib_archetypes <- function(DT_agg,
                                      cache_path = here("outputs/lovb/data/DT_contrib_archetypes.rds"),
                                      use_cache  = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    cli_alert_info("Loading cached archetypes DT_contrib from {.path {cache_path}}")
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Build archetypes DT_contrib (wide + derived metrics)")

  id_cols <- c("archetype", "Cluster")
  DT_c <- .lovb_contrib_core(DT_agg, id_cols)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT_c, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT_c)}} rows, {.val {ncol(DT_c)}} cols)")
  DT_c
}
