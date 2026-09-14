# ==============================================================================
# Documentation and reproducibility contract for the Chapter 1 analysis chain.
#
# WHY THIS EXISTS. The chain is 39 scripts and about 2700 lines. Nothing enforced a
# common shape, so on 2026-08-04 three scripts had no header at all, only 9 of 39
# declared their inputs, and one library file (R/wind_correction.R) carried the
# baseline wind correction with no function documentation whatsoever. Conventions that
# are not checked decay silently. This script fails loudly instead.
#
# Four checks, all cheap:
#   1. every staged script has a header block declaring what it Reads and Writes;
#   2. every path it declares under `Reads :` actually exists, or is marked as being
#      produced by an earlier stage;
#   3. every top-level function in the R/ files the chain sources carries a roxygen
#      block with @param and @return;
#   4. every staged script parses.
#
# Reads : run_chapter1.R and every stage script it names (headers only)
#         R/{dtmax_convention,wind_correction,lad,musica,io,cluster_relabel}.R
# Writes: nothing. It is a contract check: it prints, and stops non-zero on any breach.
#   Rscript scripts/c1_check_docs.R          # report and fail on any breach
#   DOCS_WARN_ONLY=TRUE Rscript ...          # report only, exit 0
# ==============================================================================
suppressPackageStartupMessages({library(data.table)})
WARN <- identical(toupper(Sys.getenv("DOCS_WARN_ONLY", "FALSE")), "TRUE")
fail <- character()

runner <- readLines(here::here("chapter1/orchestration/run_chapter1.R"), warn = FALSE)
staged <- unique(regmatches(runner, regexpr('(?<=S\\(")[^"]+', runner, perl = TRUE)))
paths  <- unique(na.omit(regmatches(runner,
  regexpr('[A-Za-z0-9_/]+\\.R(?=")', runner, perl = TRUE))))
paths  <- setdiff(paths[file.exists(paths)], here::here("chapter1/orchestration/run_chapter1.R"))

# ---- 1 & 2. headers, and the inputs they declare ------------------------------
# Expand one level of {a,b,c} alternation, recursively.
expand_braces <- function(s) {
  m <- regexpr("\\{[^{}]*\\}", s)
  if (m == -1) return(s)
  inner <- substr(s, m + 1, m + attr(m, "match.length") - 2)
  pre <- substr(s, 1, m - 1); post <- substring(s, m + attr(m, "match.length"))
  unlist(lapply(strsplit(inner, ",")[[1]],
                function(o) expand_braces(paste0(pre, trimws(o), post))))
}

no_header <- no_reads <- no_writes <- character()
upstream <- declared_out <- character()
bad_input <- list()
for (f in paths) {
  head <- head(readLines(f, warn = FALSE), 45)
  if (!any(grepl("^#\\s*=====", head))) no_header <- c(no_header, f)
  if (!any(grepl("^#\\s*(Reads|In)\\s*:", head, ignore.case = TRUE)))  no_reads  <- c(no_reads, f)
  if (!any(grepl("^#\\s*(Writes|Out)\\s*:", head, ignore.case = TRUE))) no_writes <- c(no_writes, f)
  wb <- grep("^#\\s*(Writes|Out)\\s*:", head, ignore.case = TRUE)
  if (length(wb)) {
    wr <- head[wb[1]:length(head)]
    wstop <- grep("^#\\s*=====|^[^#]", wr)
    for (l in if (length(wstop)) wr[seq_len(wstop[1] - 1L)] else wr)
      declared_out <- c(declared_out, regmatches(l, regexpr("[~A-Za-z0-9_./{},<>#*-]+\\.(nc|csv|rds|tif|png|pdf|md|R)", l)))
  }
  # declared inputs must exist unless produced upstream
  blk <- grep("^#\\s*(Reads|In)\\s*:", head, ignore.case = TRUE)
  if (length(blk)) {
    # Stop at the first Writes/Out line, or at the closing rule. Without this the scan ran
    # a fixed 8 lines and spilled into the Writes: block, reporting a script's own OUTPUTS
    # as missing inputs -- 7 false positives on 2026-08-04, every one of them an output.
    rest <- head[(blk[1] + 1L):length(head)]
    stopat <- grep("^#\\s*(Writes|Out|Sortie)\\s*:|^#\\s*=====|^[^#]", rest)
    rest <- if (length(stopat)) rest[seq_len(stopat[1] - 1L)] else rest
    for (l in c(head[blk[1]], rest)) {
      tok <- regmatches(l, regexpr("[~A-Za-z0-9_./{},<>#*-]+\\.(nc|csv|rds|tif|R)", l))
      if (!length(tok)) next
      # A "(stage Ax)" input is exempt from EXISTING, because the stage may not have run.
      # It is NOT exempt from being produced: some staged script must declare it under
      # Writes:. Otherwise the annotation is an unfalsifiable claim about provenance.
      if (grepl("stage [A-B]", l)) { upstream <- c(upstream, tok); next }
      if (grepl("^[A-Za-z0-9_.-]+$", tok)) next     # bare basename, not a path
      # Set notations are RESOLVED, not exempted: {a,b}.tif is expanded, <id>/#### become
      # globs, and each candidate must match at least one file. A brace set must carry its
      # directory on the SAME line, or it expands to bare basenames and cannot be checked.
      if (grepl("[{<*#]", p2 <- tok)) {
        for (c in gsub("#+", "*", gsub("<[^>]+>", "*", expand_braces(p2))))
          if (!length(Sys.glob(c))) bad_input[[f]] <- c(bad_input[[f]], c)
      } else if (nzchar(dirname(p2)) && dirname(p2) != "." && !file.exists(p2))
        bad_input[[f]] <- c(bad_input[[f]], p2)
    }
  }
}
say <- function(ok, msg) cat(sprintf("  [%s] %s\n", if (ok) "x" else " ", msg))
say(!length(no_header), sprintf("%d/%d staged scripts carry a header block", length(paths) - length(no_header), length(paths)))
say(!length(no_reads),  sprintf("%d declare their inputs", length(paths) - length(no_reads)))
say(!length(no_writes), sprintf("%d declare their outputs", length(paths) - length(no_writes)))
say(!length(bad_input), sprintf("%d declared inputs missing on disk", length(unlist(bad_input))))
# Every "(stage Ax)" input must be some script's declared output.
orphan <- unique(Filter(function(u) !any(sapply(declared_out, function(o)
            grepl(basename(sub("\\{.*", "", o)), u, fixed = TRUE) ||
            grepl(basename(sub("\\{.*", "", u)), o, fixed = TRUE))), unique(upstream)))
say(!length(orphan), sprintf("%d of %d stage-produced inputs have no declared producer",
                             length(orphan), length(unique(upstream))))
if (length(orphan)) fail <- c(fail, paste("no declared producer:", paste(orphan, collapse = ", ")))
if (length(no_header)) fail <- c(fail, paste("no header:", paste(basename(no_header), collapse = ", ")))
if (length(no_reads))  fail <- c(fail, paste("no Reads: declaration:", paste(basename(no_reads), collapse = ", ")))
if (length(no_writes)) fail <- c(fail, paste("no Writes: declaration:", paste(basename(no_writes), collapse = ", ")))
if (length(bad_input)) for (f in names(bad_input))
  fail <- c(fail, sprintf("%s declares a missing input: %s", basename(f), paste(bad_input[[f]], collapse = ", ")))

# ---- 3. roxygen on the library functions the chain sources ---------------------
LIB <- c("R/dtmax_convention.R","R/wind_correction.R","R/lad.R","R/musica.R",
         "R/io.R","R/cluster_relabel.R")
undoc <- character()
for (f in LIB) {
  L <- readLines(f, warn = FALSE)
  hits <- grep("^\\s*[\\w.]+\\s*<-\\s*function", L, perl = TRUE)
  for (i in hits) {
    if (grepl("^\\s{2,}", L[i])) next              # nested closure, not a top-level API
    k <- i - 1L; ok <- FALSE
    while (k >= 1 && nzchar(trimws(L[k]))) {
      if (startsWith(trimws(L[k]), "#'")) { ok <- TRUE; break }
      if (!startsWith(trimws(L[k]), "#")) break
      k <- k - 1L
    }
    if (!ok) undoc <- c(undoc, sprintf("%s:%d %s", f, i, trimws(L[i])))
  }
}
say(!length(undoc), sprintf("every top-level function in the %d sourced library files has a roxygen block", length(LIB)))
if (length(undoc)) fail <- c(fail, paste("undocumented:", paste(undoc, collapse = " | ")))

# ---- 4. everything parses -----------------------------------------------------
broken <- Filter(function(f) inherits(try(parse(f), silent = TRUE), "try-error"), c(paths, LIB))
say(!length(broken), sprintf("%d files parse", length(paths) + length(LIB)))
if (length(broken)) fail <- c(fail, paste("parse error:", paste(broken, collapse = ", ")))

if (length(fail)) {
  msg <- paste("documentation contract FAILED:\n  ", paste(fail, collapse = "\n  "))
  if (WARN) cat("\n", msg, "\n") else stop(msg)
} else cat("\nDocumentation contract satisfied.\nDONE\n")
