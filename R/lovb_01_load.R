# ==============================================================================
# LOVB analysis — Module 1 : load daily Delta_Tmax for selected coalitions
#
# Reads existing factorial simulations (no MuSICA re-run).
# Uses extract_deltatmax_one() from R/musica.R (interpolates Tair at z=1m,
# default CFG$tair_target_height). That utility uses dplyr internally; we
# convert its output to data.table immediately. No dplyr in LOVB functions.
#
# Output : data.table cached as outputs/lovb/data/DT_daily_<source>.rds
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ncdf4)
  library(lubridate)
  library(musica.tools)   # provides get_variable() used by R/musica.R::get_tair_at_z
  library(rmusica)        # provides force_utc_nc() used by R/io.R
  library(dplyr)          # required by extract_deltatmax_one() internals
})

# Reuse existing project utilities (R/musica.R uses dplyr internally — invisible
# to our pipeline because we setDT() the result).
source(here("R/config.R"))
source(here("R/io.R"))
source(here("R/lad.R"))
source(here("R/musica.R"))    # provides extract_deltatmax_one() with z=1m interpolation

# ---- Bit-code conventions ----------------------------------------------------

# 10 coalitions of interest : NULL + 4 LVA singletons + 4 LOVB + REF
LOVB_COALITIONS <- c(
  "0000",                                 # NULL  (LVA reference)
  "0001", "0010", "0100", "1000",         # 4 LVA singletons (bit=1 on one var only)
  "0111", "1011", "1101", "1110",         # 4 LOVB (bit=0 on one var only, others=1)
  "1111"                                  # REF (all real)
)

# Suffix convention from existing pipeline (R/scenarios.R) :
#   pos1=LAI : m=mean,    r=real
#   pos2=Hmax: m=mean,    r=real
#   pos3=fCover: a=absent (=1), r=real
#   pos4=LAD : u=uniform, r=real
lovb_bit_to_suffix <- function(bit_code) {
  stopifnot(is.character(bit_code), nchar(bit_code) == 4L)
  b <- strsplit(bit_code, "", fixed = TRUE)[[1]]
  paste(
    ifelse(b[1] == "1", "r", "m"),
    ifelse(b[2] == "1", "r", "m"),
    ifelse(b[3] == "1", "r", "a"),
    ifelse(b[4] == "1", "r", "u"),
    sep = "_"
  )
}
lovb_bit_to_h1f  <- function(b) paste0("H1F_",  lovb_bit_to_suffix(b))
lovb_bit_to_arch <- function(b) paste0("ARCH_", lovb_bit_to_suffix(b))

# ---- Helper : load one factorial directory of NCs --> data.table -------------
.lovb_extract_dir <- function(scenario_dir, bit_code, df_macro, date_seq,
                                xy_from_filename = TRUE,
                                archetype_label  = NULL) {
  nc_files <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE)
  nc_files <- nc_files[file.size(nc_files) > 1e6]   # skip truncated
  if (length(nc_files) == 0L) return(NULL)

  rows <- vector("list", length(nc_files))
  for (j in seq_along(nc_files)) {
    f   <- nc_files[j]
    res <- tryCatch(extract_deltatmax_one(f, df_macro, date_seq),
                     error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0L) next
    dt <- as.data.table(res)

    if (xy_from_filename) {
      bn <- basename(f)
      dt[, x := as.integer(sub(".*_X(\\d+)_Y\\d+\\.nc$", "\\1", bn))]
      dt[, y := as.integer(sub(".*_X\\d+_Y(\\d+)\\.nc$", "\\1", bn))]
    } else {
      dt[, archetype := archetype_label]
    }
    dt[, bit_code := bit_code]
    rows[[j]] <- dt
  }
  rbindlist(rows, fill = TRUE)
}

# ---- 1. cLHS : 400 plots x N coalitions x 122 days ---------------------------
lovb_load_daily_clhs <- function(coalitions = LOVB_COALITIONS,
                                  cache_path = here("outputs/lovb/data/DT_daily_cLHS.rds"),
                                  use_cache  = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    cli_alert_info("Loading cached cLHS daily DT from {.path {cache_path}}")
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Extract daily Delta_Tmax — cLHS factorial")
  cli_alert("{length(coalitions)} coalitions x 400 plots NCs")

  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  out <- vector("list", length(coalitions))
  for (i in seq_along(coalitions)) {
    bit <- coalitions[i]
    sc  <- lovb_bit_to_h1f(bit)
    dir <- here("out_files/H1_factorial", sc)
    if (!dir.exists(dir)) {
      cli_alert_warning("Missing dir for {bit} ({sc})")
      next
    }
    cli_alert("[{i}/{length(coalitions)}] {bit} <-> {sc}")
    out[[i]] <- .lovb_extract_dir(dir, bit, df_macro, CFG$date_seq,
                                    xy_from_filename = TRUE)
  }
  DT <- rbindlist(out, fill = TRUE)
  DT <- DT[, .(x, y, date, bit_code, Delta_Tmax)]
  setkey(DT, x, y, bit_code, date)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT)}} rows)")
  DT
}

# ---- 2. Archetypes : 4 archetypes x N coalitions x 122 days ------------------
lovb_load_daily_archetypes <- function(coalitions = LOVB_COALITIONS,
                                        cache_path = here("outputs/lovb/data/DT_daily_archetypes.rds"),
                                        use_cache  = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    cli_alert_info("Loading cached archetypes daily DT from {.path {cache_path}}")
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Extract daily Delta_Tmax — archetypes")
  cli_alert("4 archetypes x {length(coalitions)} coalitions NCs")

  df_macro   <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  archetypes <- c("Arch_C1", "Arch_C2", "Arch_C3", "Arch_C4")
  out <- list()
  for (arch in archetypes) {
    dir <- here("out_files/H1_archetypes", arch)
    if (!dir.exists(dir)) { cli_alert_warning("Missing : {.path {dir}}"); next }
    for (bit in coalitions) {
      sc <- lovb_bit_to_arch(bit)
      nc <- file.path(dir, sprintf("musica_out_%s_%s.nc", arch, sc))
      if (!file.exists(nc) || file.size(nc) < 1e6) {
        cli_alert_warning("Missing/truncated : {basename(nc)}"); next
      }
      res <- tryCatch(extract_deltatmax_one(nc, df_macro, CFG$date_seq),
                       error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0L) next
      dt <- as.data.table(res)
      dt[, archetype := arch]
      dt[, bit_code  := bit]
      out[[length(out) + 1L]] <- dt
    }
  }
  DT <- rbindlist(out, fill = TRUE)
  DT <- DT[, .(archetype, date, bit_code, Delta_Tmax)]
  setkey(DT, archetype, bit_code, date)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT)}} rows)")
  DT
}
