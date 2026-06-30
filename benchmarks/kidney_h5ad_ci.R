suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(spacexrfast)
})

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

stop_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this CI example.", call. = FALSE)
  }
}

read_first_h5_path <- function(path, candidates) {
  for (candidate in candidates) {
    value <- tryCatch(rhdf5::h5read(path, candidate), error = function(e) NULL)
    if (!is.null(value)) {
      return(value)
    }
  }
  NULL
}

read_h5ad_axis_names <- function(path, axis, n) {
  values <- read_first_h5_path(path, c(paste0("/", axis, "/_index"), paste0("/", axis, "/index")))
  if (is.null(values)) {
    values <- paste0(axis, "_", seq_len(n))
  }
  make.unique(as.character(values))
}

read_h5ad_matrix <- function(path, max_cells = 1000) {
  stop_missing("rhdf5")
  listing <- rhdf5::h5ls(path, recursive = 3)
  has_group <- function(group) {
    any(listing$group == dirname(group) & listing$name == basename(group) & listing$otype == "H5I_GROUP")
  }
  has_dataset <- function(dataset) {
    any(listing$group == dirname(dataset) & listing$name == basename(dataset) & listing$otype == "H5I_DATASET")
  }

  obs_names <- read_h5ad_axis_names(path, "obs", 0)
  var_names <- read_h5ad_axis_names(path, "var", 0)
  n_obs <- length(obs_names)
  n_var <- length(var_names)
  n_cells <- min(max_cells, n_obs)

  matrix_group <- if (has_group("/layers/counts")) "/layers/counts" else "/X"

  if (has_group(matrix_group) && has_dataset(paste0(matrix_group, "/data"))) {
    indptr <- as.integer(rhdf5::h5read(path, paste0(matrix_group, "/indptr")))
    sparse_data_path <- paste0(matrix_group, "/data")
    sparse_indices_path <- paste0(matrix_group, "/indices")

    if (length(indptr) == n_obs + 1L) {
      nnz <- indptr[n_cells + 1L]
      if (nnz == 0L) {
        counts <- Matrix(0, nrow = n_var, ncol = n_cells, sparse = TRUE)
      } else {
        data <- as.numeric(rhdf5::h5read(path, sparse_data_path, index = list(seq_len(nnz))))
        indices <- as.integer(rhdf5::h5read(path, sparse_indices_path, index = list(seq_len(nnz)))) + 1L
        ptr <- indptr[seq_len(n_cells + 1L)]
        counts <- sparseMatrix(
          i = indices,
          j = rep(seq_len(n_cells), diff(ptr)),
          x = data,
          dims = c(n_var, n_cells)
        )
      }
    } else if (length(indptr) == n_var + 1L) {
      data <- as.numeric(rhdf5::h5read(path, sparse_data_path))
      indices <- as.integer(rhdf5::h5read(path, sparse_indices_path)) + 1L
      counts <- sparseMatrix(
        i = rep(seq_len(n_var), diff(indptr)),
        j = indices,
        x = data,
        dims = c(n_var, n_obs)
      )[, seq_len(n_cells), drop = FALSE]
    } else {
      stop("Could not infer H5AD sparse orientation from indptr length.", call. = FALSE)
    }
  } else if (has_dataset(matrix_group)) {
    dims <- listing$dim[listing$group == dirname(matrix_group) & listing$name == basename(matrix_group)][1]
    dims <- as.integer(strsplit(dims, " x ", fixed = TRUE)[[1]])
    if (length(dims) != 2L) {
      stop("Could not infer H5AD dense matrix dimensions.", call. = FALSE)
    }
    if (dims[1] == n_obs) {
      dense <- rhdf5::h5read(path, matrix_group, index = list(seq_len(n_cells), seq_len(n_var)))
      counts <- Matrix(t(dense), sparse = TRUE)
    } else {
      dense <- rhdf5::h5read(path, matrix_group, index = list(seq_len(n_var), seq_len(n_cells)))
      counts <- Matrix(dense, sparse = TRUE)
    }
  } else {
    stop("Could not find /layers/counts or /X in H5AD file.", call. = FALSE)
  }

  rownames(counts) <- make.unique(var_names)
  colnames(counts) <- make.unique(obs_names[seq_len(n_cells)])
  counts
}

build_pseudo_reference <- function(counts, max_types, cells_per_type = 25L) {
  type_names <- paste0("pseudo_type_", seq_len(max_types))
  profile <- matrix(
    1e-4,
    nrow = nrow(counts),
    ncol = length(type_names),
    dimnames = list(rownames(counts), type_names)
  )
  gene_groups <- split(seq_len(nrow(counts)), rep(seq_along(type_names), length.out = nrow(counts)))
  for (i in seq_along(type_names)) {
    profile[gene_groups[[i]], i] <- Matrix::rowMeans(counts[gene_groups[[i]], , drop = FALSE]) + 1e-3
  }
  profile <- sweep(profile, 2, colSums(profile), "/")

  ref_counts <- do.call(cbind, lapply(seq_along(type_names), function(i) {
    replicate(cells_per_type, rpois(nrow(profile), lambda = profile[, i] * 5000))
  }))
  rownames(ref_counts) <- rownames(profile)
  colnames(ref_counts) <- paste0("ref_", seq_len(ncol(ref_counts)))
  cell_types <- factor(rep(type_names, each = cells_per_type), levels = type_names)
  names(cell_types) <- colnames(ref_counts)
  Reference(Matrix(ref_counts, sparse = TRUE), cell_types, min_UMI = 10)
}

matrix_path <- path.expand(Sys.getenv("KIDNEY_H5AD", "~/Documents/10x-examples/atlases/rcc-kidney-census-v1_reference.h5ad"))
max_cells <- as.integer(Sys.getenv("KIDNEY_EXAMPLE_CELLS", "1000"))
max_genes <- as.integer(Sys.getenv("KIDNEY_EXAMPLE_GENES", "80"))
max_types <- as.integer(Sys.getenv("KIDNEY_EXAMPLE_TYPES", "4"))
max_cores <- as.integer(Sys.getenv("KIDNEY_EXAMPLE_CORES", "2"))
weight_tolerance <- as.numeric(Sys.getenv("KIDNEY_EXAMPLE_WEIGHT_TOL", "2e-4"))
set.seed(as.integer(Sys.getenv("KIDNEY_EXAMPLE_SEED", "1")))

if (!file.exists(matrix_path)) {
  stop("KIDNEY_H5AD does not exist: ", matrix_path, call. = FALSE)
}

log_step("reading H5AD matrix: ", matrix_path)
counts <- read_h5ad_matrix(matrix_path, max_cells = max_cells)
counts <- counts[Matrix::rowSums(counts) > 0, Matrix::colSums(counts) > 0, drop = FALSE]
genes <- head(rownames(counts)[order(Matrix::rowSums(counts), decreasing = TRUE)], max_genes)
cells <- head(colnames(counts)[order(Matrix::colSums(counts), decreasing = TRUE)], max_cells)
counts <- counts[genes, cells, drop = FALSE]

log_step("building RCTD object for ", ncol(counts), " cells, ", nrow(counts), " genes, ", max_types, " pseudo types")
coords <- data.frame(x = seq_along(cells), y = seq_along(cells), row.names = cells)
spatial <- SpatialRNA(coords, counts, require_int = FALSE)
reference <- build_pseudo_reference(counts, max_types = max_types)

timer <- proc.time()
rctd <- create.RCTD(
  spatial, reference,
  max_cores = max_cores, UMI_min = 1, counts_MIN = 1, UMI_min_sigma = 1,
  CELL_MIN_INSTANCE = 5, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20
)
rctd <- fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)
preprocess_elapsed <- unname((proc.time() - timer)[["elapsed"]])
log_step("preprocessing elapsed seconds: ", round(preprocess_elapsed, 3))

log_step("running upstream spacexr doublet mode")
vanilla_time <- system.time(vanilla <- spacexr:::fitPixels(rctd, doublet_mode = "doublet"))

log_step("running fast doublet mode with progress enabled")
fast_time <- system.time(fast <- run.RCTD.fast.doublet(
  rctd,
  max_cores = max_cores,
  chunk_size = "auto",
  progress = TRUE
))

spot_class_equal <- identical(
  as.character(vanilla@results$results_df$spot_class),
  as.character(fast@results$results_df$spot_class)
)
first_type_equal <- identical(
  as.character(vanilla@results$results_df$first_type),
  as.character(fast@results$results_df$first_type)
)
second_type_equal <- identical(
  as.character(vanilla@results$results_df$second_type),
  as.character(fast@results$results_df$second_type)
)
weights_diff <- max(abs(as.matrix(vanilla@results$weights) - as.matrix(fast@results$weights)))
weights_doublet_diff <- max(abs(as.matrix(vanilla@results$weights_doublet) - as.matrix(fast@results$weights_doublet)))

speedup <- unname(vanilla_time[["elapsed"]] / fast_time[["elapsed"]])
cat("cells,", ncol(counts), "\n", sep = "")
cat("genes,", nrow(counts), "\n", sep = "")
cat("types,", max_types, "\n", sep = "")
cat("vanilla_elapsed_seconds,", unname(vanilla_time[["elapsed"]]), "\n", sep = "")
cat("fast_elapsed_seconds,", unname(fast_time[["elapsed"]]), "\n", sep = "")
cat("speedup,", speedup, "\n", sep = "")
cat("spot_class_equal,", spot_class_equal, "\n", sep = "")
cat("first_type_equal,", first_type_equal, "\n", sep = "")
cat("second_type_equal,", second_type_equal, "\n", sep = "")
cat("weights_max_abs_diff,", weights_diff, "\n", sep = "")
cat("weights_doublet_max_abs_diff,", weights_doublet_diff, "\n", sep = "")

if (!spot_class_equal || !first_type_equal || !second_type_equal) {
  stop("Fast labels do not match upstream RCTD on the kidney H5AD example.", call. = FALSE)
}
if (!is.finite(weights_diff) || weights_diff > weight_tolerance) {
  stop("Fast all-type weights exceed tolerance: ", weights_diff, call. = FALSE)
}
if (!is.finite(weights_doublet_diff) || weights_doublet_diff > weight_tolerance) {
  stop("Fast doublet weights exceed tolerance: ", weights_doublet_diff, call. = FALSE)
}
