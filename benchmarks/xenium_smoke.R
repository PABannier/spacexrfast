suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
  library(spacexrfast)
})

read_10x_h5 <- function(path) {
  if (requireNamespace("Seurat", quietly = TRUE) && requireNamespace("hdf5r", quietly = TRUE)) {
    counts <- Seurat::Read10X_h5(path)
    if (is.list(counts)) {
      counts <- counts[[1]]
    }
    return(counts)
  }

  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Install either hdf5r+Seurat or rhdf5 to read 10x HDF5 matrices.", call. = FALSE)
  }
  data <- rhdf5::h5read(path, "/matrix/data")
  indices <- rhdf5::h5read(path, "/matrix/indices")
  indptr <- rhdf5::h5read(path, "/matrix/indptr")
  shape <- as.integer(rhdf5::h5read(path, "/matrix/shape"))
  barcodes <- rhdf5::h5read(path, "/matrix/barcodes")
  feature_names <- rhdf5::h5read(path, "/matrix/features/name")
  counts <- new(
    "dgCMatrix",
    x = as.numeric(data),
    i = as.integer(indices),
    p = as.integer(indptr),
    Dim = shape
  )
  rownames(counts) <- make.unique(as.character(feature_names))
  colnames(counts) <- as.character(barcodes)
  counts
}

matrix_path <- Sys.getenv(
  "XENIUM_MATRIX",
  "~/Documents/10x-examples/kidney/outs/cell_feature_matrix.h5"
)
matrix_path <- path.expand(matrix_path)
max_cells <- as.integer(Sys.getenv("XENIUM_SMOKE_CELLS", "50"))
max_genes <- as.integer(Sys.getenv("XENIUM_SMOKE_GENES", "80"))
max_cores <- as.integer(Sys.getenv("XENIUM_SMOKE_CORES", "2"))
max_types <- as.integer(Sys.getenv("XENIUM_SMOKE_TYPES", "3"))
set.seed(as.integer(Sys.getenv("XENIUM_SMOKE_SEED", "1")))

message("Reading ", matrix_path)
counts <- read_10x_h5(matrix_path)
counts <- counts[Matrix::rowSums(counts) > 0, Matrix::colSums(counts) > 0, drop = FALSE]
genes <- head(rownames(counts)[order(Matrix::rowSums(counts), decreasing = TRUE)], max_genes)
cells <- head(colnames(counts)[order(Matrix::colSums(counts), decreasing = TRUE)], max_cells)
counts <- counts[genes, cells, drop = FALSE]

coords <- data.frame(x = seq_along(cells), y = seq_along(cells), row.names = cells)
spatial <- SpatialRNA(coords, counts, require_int = FALSE)

type_names <- paste0("pseudo_type_", seq_len(max_types))
profile <- matrix(1e-4, nrow = nrow(counts), ncol = length(type_names),
  dimnames = list(rownames(counts), type_names)
)
gene_groups <- split(seq_len(nrow(counts)), rep(seq_along(type_names), length.out = nrow(counts)))
for (i in seq_along(type_names)) {
  profile[gene_groups[[i]], i] <- Matrix::rowMeans(counts[gene_groups[[i]], , drop = FALSE]) + 1e-3
}
profile <- sweep(profile, 2, colSums(profile), "/")

ref_counts <- do.call(cbind, lapply(seq_along(type_names), function(i) {
  replicate(30, rpois(nrow(profile), lambda = profile[, i] * 5000))
}))
rownames(ref_counts) <- rownames(profile)
colnames(ref_counts) <- paste0("ref_", seq_len(ncol(ref_counts)))
cell_types <- factor(rep(type_names, each = 30), levels = type_names)
names(cell_types) <- colnames(ref_counts)
reference <- Reference(Matrix(ref_counts, sparse = TRUE), cell_types, min_UMI = 10)

rctd <- create.RCTD(
  spatial, reference,
  max_cores = max_cores, UMI_min = 1, counts_MIN = 1, UMI_min_sigma = 1,
  CELL_MIN_INSTANCE = 5, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20
)
rctd <- fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)

vanilla_time <- system.time(vanilla <- spacexr:::fitPixels(rctd, doublet_mode = "doublet"))
fast_time <- system.time(fast <- run.RCTD.fast.doublet(rctd, max_cores = max_cores, chunk_size = "auto"))

cat("cells,", ncol(counts), "\n", sep = "")
cat("genes,", nrow(counts), "\n", sep = "")
cat("types,", length(type_names), "\n", sep = "")
cat("vanilla_elapsed_seconds,", unname(vanilla_time[["elapsed"]]), "\n", sep = "")
cat("fast_elapsed_seconds,", unname(fast_time[["elapsed"]]), "\n", sep = "")
cat(
  "spot_class_equal,",
  identical(as.character(vanilla@results$results_df$spot_class), as.character(fast@results$results_df$spot_class)),
  "\n",
  sep = ""
)
cat(
  "first_type_equal,",
  identical(as.character(vanilla@results$results_df$first_type), as.character(fast@results$results_df$first_type)),
  "\n",
  sep = ""
)
cat(
  "second_type_equal,",
  identical(as.character(vanilla@results$results_df$second_type), as.character(fast@results$results_df$second_type)),
  "\n",
  sep = ""
)
cat(
  "weights_max_abs_diff,",
  max(abs(as.matrix(vanilla@results$weights) - as.matrix(fast@results$weights))),
  "\n",
  sep = ""
)
cat(
  "weights_doublet_max_abs_diff,",
  max(abs(as.matrix(vanilla@results$weights_doublet) - as.matrix(fast@results$weights_doublet))),
  "\n",
  sep = ""
)
