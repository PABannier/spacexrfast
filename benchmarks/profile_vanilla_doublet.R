suppressPackageStartupMessages({
  library(spacexr)
  library(spacexrfast)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript benchmarks/profile_vanilla_doublet.R prepared_rctd.rds", call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

rctd <- readRDS(args[[1]])
vanilla <- profile_vanilla_doublet(rctd)
fast_time <- system.time({
  fast <- run.RCTD.fast.doublet(rctd, max_cores = rctd@config$max_cores %||% 8)
})

cat("vanilla_elapsed_seconds,", vanilla$elapsed_seconds, "\n", sep = "")
cat("fast_elapsed_seconds,", unname(fast_time[["elapsed"]]), "\n", sep = "")
cat(
  "spot_class_equal,",
  identical(
    as.character(vanilla$rctd@results$results_df$spot_class),
    as.character(fast@results$results_df$spot_class)
  ),
  "\n",
  sep = ""
)
cat(
  "first_type_equal,",
  identical(
    as.character(vanilla$rctd@results$results_df$first_type),
    as.character(fast@results$results_df$first_type)
  ),
  "\n",
  sep = ""
)
cat(
  "second_type_equal,",
  identical(
    as.character(vanilla$rctd@results$results_df$second_type),
    as.character(fast@results$results_df$second_type)
  ),
  "\n",
  sep = ""
)
