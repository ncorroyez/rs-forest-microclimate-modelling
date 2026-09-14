# ==============================================================================
# Isolation test: STATIC_ALS under v3.2.3 with ABL_flag='none' (NO yoyo, standard
# forcing, parametric phenology = the EXACT v3.2.0 setup but the v3.2.3 binary).
# If the inter-plot ΔTmax spread compresses like the iter run (SD~0.48 vs v3.2.0
# 0.90), the over-smoothing is the BINARY, not the iter/PBLH setup.
#   Rscript run_v323none_staticals.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(parallel); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC <- CFG_C3$forcing_file                         # standard forcing (no PBLH), as v3.2.0
ABL  <- list("abl_flag" = '"none"')                 # no yoyo
NCDIR <- file.path(CFG_C3$out_dir, "nc_v323none", "STATIC_ALS"); dir.create(NCDIR, recursive=TRUE, showWarnings=FALSE)

prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
# STATIC_ALS exactly as v3.2.0: phenology_fn=NULL -> run_musica_one parametric fallback
.lf <- function(col) function(pr) as.numeric(pr[[col]])
SC <- list(name="STATIC_ALS", lai_fn=.lf("LAI_ALS"), hmax_fn=.lf("Hmax"),
           fcover_fn=.lf("fCover"), lad_fn=make_lad_real, phenology_fn=NULL)

run_one <- function(i){ prow <- as.data.frame(df[i,])
  nc <- file.path(NCDIR, sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(nc) && file.size(nc) > 1000) return("skip")
  if (file.exists(nc)) file.remove(nc)
  ok <- tryCatch({ run_musica_one(prow, SC, nc, FORC, BIN, extra_setup = ABL)
                   file.exists(nc) && file.size(nc) > 1000 }, error=function(e) FALSE)
  if (ok) "OK" else "FAIL" }
t0 <- Sys.time()
res <- unlist(mclapply(seq_len(nrow(df)), run_one, mc.cores = 10L, mc.preschedule = FALSE))
cat(sprintf("DONE in %.1f min | %s\n", as.numeric(difftime(Sys.time(),t0,units="mins")),
            paste(names(table(res)), table(res), sep="=", collapse=" ")))
