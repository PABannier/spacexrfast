suppressPackageStartupMessages({
  library(spacexrfast)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript benchmarks/thread_scaling.R prepared_rctd_after_sigma.rds", call. = FALSE)
}

cores <- as.integer(strsplit(Sys.getenv("THREAD_SCALING_CORES", "1,4,8"), ",", fixed = TRUE)[[1]])
chunk_size <- Sys.getenv("THREAD_SCALING_CHUNK_SIZE", "auto")
if (!identical(chunk_size, "auto")) {
  chunk_size <- as.integer(chunk_size)
}

sample_threads <- function(pid, output, stopfile, interval = 0.2) {
  con <- file(output, open = "wt")
  on.exit(close(con), add = TRUE)
  repeat {
    if (file.exists(stopfile)) break
    nlwp <- suppressWarnings(system2("ps", c("-o", "nlwp=", "-p", pid), stdout = TRUE, stderr = FALSE))
    if (length(nlwp) == 1 && nzchar(trimws(nlwp))) {
      writeLines(paste(Sys.time(), trimws(nlwp)), con)
      flush(con)
    }
    Sys.sleep(interval)
  }
}

run_with_thread_sampler <- function(expr) {
  if (.Platform$OS.type != "unix" || !requireNamespace("parallel", quietly = TRUE)) {
    elapsed <- system.time(value <- force(expr))
    return(list(value = value, elapsed = unname(elapsed[["elapsed"]]), max_nlwp = NA_integer_))
  }

  samples <- tempfile("nlwp_samples_")
  stopfile <- tempfile("nlwp_stop_")
  sampler <- parallel::mcparallel(sample_threads(Sys.getpid(), samples, stopfile), silent = TRUE)
  collected <- FALSE
  on.exit({
    file.create(stopfile)
    if (!collected) suppressWarnings(parallel::mccollect(sampler, wait = TRUE))
    unlink(c(samples, stopfile), force = TRUE)
  }, add = TRUE)

  elapsed <- system.time(value <- force(expr))
  file.create(stopfile)
  suppressWarnings(parallel::mccollect(sampler, wait = TRUE))
  collected <- TRUE
  values <- if (file.exists(samples)) {
    lines <- readLines(samples, warn = FALSE)
    as.integer(sub("^.* ([0-9]+)$", "\\1", lines))
  } else {
    integer()
  }
  max_nlwp <- if (length(values) == 0 || all(is.na(values))) NA_integer_ else max(values, na.rm = TRUE)
  list(value = value, elapsed = unname(elapsed[["elapsed"]]), max_nlwp = max_nlwp)
}

rctd <- readRDS(args[[1]])
cat("cells,", ncol(rctd@spatialRNA@counts), "\n", sep = "")
cat("genes,", length(rctd@internal_vars$gene_list_reg), "\n", sep = "")
cat("types,", length(rctd@cell_type_info$renorm[[2]]), "\n", sep = "")

for (core_count in cores) {
  result <- run_with_thread_sampler(
    run.RCTD.fast.doublet(rctd, max_cores = core_count, chunk_size = chunk_size)
  )
  cat("max_cores,", core_count, "\n", sep = "")
  cat("elapsed_seconds,", result$elapsed, "\n", sep = "")
  cat("max_nlwp,", result$max_nlwp, "\n", sep = "")
}
