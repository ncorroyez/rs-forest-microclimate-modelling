# ==============================================================================
# V11 (model-3.2.3) — Run 4 archetypes × 16 coalitions = 64 sims.
# Output : out_files/H1_archetypes_v11_model323/Arch_C{1,2,3,4}/
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(future); library(furrr)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

NEW_BINARY <- here::here("in_files/model-3.2.3/musica")
stopifnot(file.exists(NEW_BINARY))
cli_h1("V11 archetypes — model-3.2.3")
cli_alert("Binary : {NEW_BINARY}")

N_WORKERS <- 4
plan(multisession, workers = N_WORKERS)

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
df_archetypes <- make_synthetic_archetypes(as.data.frame(df_floor))
cli_alert("Archetypes :")
print(df_archetypes[, c("Cluster", "LAI", "Hmax", "fCover")])

fac_scs <- build_factorial_scenarios_archetypes(df_floor)
out_root <- here::here("out_files/H1_archetypes_v11_model323")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

jobs <- list()
for (i in seq_len(nrow(df_archetypes))) {
  arch <- df_archetypes[i, , drop = FALSE]
  arch_lbl <- sprintf("Arch_C%s", arch$Cluster)
  out_dir <- file.path(out_root, arch_lbl)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (sc_nm in names(fac_scs)) {
    sc <- fac_scs[[sc_nm]]
    out_nc <- file.path(out_dir, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs) + 1L]] <- list(arch = arch, sc = sc, out_nc = out_nc,
                                      lbl = paste(arch_lbl, sc_nm))
  }
}
cli_alert("Archetype jobs to run : {length(jobs)}")
if (length(jobs) == 0) {
  cli_alert_success("Nothing to do — all 64 sims already present.")
  quit(save = "no")
}

t0 <- Sys.time()
future_walk(jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  source(here::here("R/h1_shapley_archetypes.R"))
  NEW_BINARY <- here::here("in_files/model-3.2.3/musica")
  tryCatch(
    run_musica_one(job$arch, job$sc, job$out_nc,
                   CFG$forcing_file, NEW_BINARY),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message))
  )
})
dt <- difftime(Sys.time(), t0, units = "mins")
cli_alert_success("Archetypes done in {round(as.numeric(dt), 1)} min")

# Final count
total_ok <- sum(file.exists(unlist(lapply(jobs, function(j) j$out_nc))) &
                file.size(unlist(lapply(jobs, function(j) j$out_nc))) > 1e6)
cli_alert("OK : {total_ok}/{length(jobs)}")
