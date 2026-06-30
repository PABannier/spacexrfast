library(spacexrfast)

test_that("C++ doublet chunk matches upstream categorical outputs on synthetic data", {
  skip_if_not_installed("spacexr")

  set.seed(1)
  ns <- asNamespace("spacexr")
  q_list <- readRDS(system.file("extdata/Qmat/Q_mat_1.rds", package = "spacexr"))
  q_mat <- q_list[[1]]
  x_vals <- readRDS(system.file("extdata/Qmat/X_vals.rds", package = "spacexr"))
  spacexr:::set_likelihood_vars(q_mat, x_vals)
  sq_mat <- get("SQ_mat", envir = ns)
  k_val <- get("K_val", envir = ns)

  n_genes <- 30
  n_types <- 4
  n_cells <- 5
  genes <- paste0("g", seq_len(n_genes))
  types <- paste0("t", seq_len(n_types))
  profiles <- matrix(runif(n_genes * n_types, 0.0001, 0.02), n_genes, n_types,
    dimnames = list(genes, types)
  )
  profiles <- sweep(profiles, 2, colSums(profiles), "/")
  cell_type_info <- list(profiles, types, n_types)
  beads <- matrix(rpois(n_cells * n_genes, lambda = 2), n_cells, n_genes,
    dimnames = list(paste0("c", seq_len(n_cells)), genes)
  )
  nUMI <- rowSums(beads)
  nUMI[nUMI == 0] <- 1
  class_df <- data.frame(class = types, row.names = types)

  all_weights <- matrix(NA_real_, n_cells, n_types, dimnames = list(rownames(beads), types))
  conv_all <- logical(n_cells)
  vanilla <- vector("list", n_cells)
  for (i in seq_len(n_cells)) {
    vanilla[[i]] <- spacexr:::process_bead_doublet(
      cell_type_info, genes, nUMI[[i]], beads[i, ],
      class_df = class_df, constrain = FALSE, MIN.CHANGE = 0.001,
      CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20
    )
    all_weights[i, ] <- vanilla[[i]]$all_weights[types]
    conv_all[[i]] <- vanilla[[i]]$conv_all
  }

  fast <- spacexrfast:::fast_doublet_chunk_cpp(
    beads = beads,
    nUMI = as.numeric(nUMI),
    ref_profiles = data.matrix(profiles),
    class_ids = seq_len(n_types),
    Q_mat = q_mat,
    SQ_mat = sq_mat,
    X_vals = x_vals,
    K_val = k_val,
    constrain = FALSE,
    min_change = 0.001,
    confidence_threshold = 5,
    doublet_threshold = 20,
    return_diagnostics = FALSE,
    max_cores = 2
  )

  expect_identical(fast$spot_class, vapply(vanilla, function(x) as.character(x$spot_class), character(1)))
  expect_identical(fast$first_type, match(vapply(vanilla, `[[`, character(1), "first_type"), types))
  expect_identical(fast$second_type, match(vapply(vanilla, `[[`, character(1), "second_type"), types))
  expect_lt(max(abs(fast$all_weights - all_weights)), 1e-4)
  expect_lt(max(abs(fast$first_weight - vapply(vanilla, function(x) unname(x$doublet_weights[[1]]), numeric(1)))), 1e-4)
  expect_lt(max(abs(fast$second_weight - vapply(vanilla, function(x) unname(x$doublet_weights[[2]]), numeric(1)))), 1e-4)
  expect_lt(max(abs(fast$min_score - vapply(vanilla, `[[`, numeric(1), "min_score"))), 1e-3)
})

test_that("run.RCTD.fast.doublet preserves the doublet result schema", {
  skip_if_not_installed("spacexr")

  set.seed(2)
  n_genes <- 60
  n_types <- 3
  n_ref_per <- 35
  n_spots <- 6
  genes <- paste0("g", seq_len(n_genes))
  types <- paste0("type", seq_len(n_types))
  base <- matrix(0.1, n_genes, n_types)
  for (type in seq_len(n_types)) {
    base[((type - 1) * 15 + 1):(type * 15), type] <- 5
  }

  ref_counts <- do.call(cbind, lapply(seq_len(n_types), function(type) {
    replicate(n_ref_per, rpois(n_genes, lambda = base[, type] * 25))
  }))
  rownames(ref_counts) <- genes
  colnames(ref_counts) <- paste0("r", seq_len(ncol(ref_counts)))
  cell_types <- factor(rep(types, each = n_ref_per), levels = types)
  names(cell_types) <- colnames(ref_counts)
  reference <- spacexr::Reference(Matrix::Matrix(ref_counts, sparse = TRUE), cell_types, min_UMI = 10)

  spot_types <- sample(seq_len(n_types), n_spots, replace = TRUE)
  spatial_counts <- sapply(spot_types, function(type) rpois(n_genes, lambda = base[, type] * 35))
  rownames(spatial_counts) <- genes
  colnames(spatial_counts) <- paste0("s", seq_len(n_spots))
  coords <- data.frame(x = runif(n_spots), y = runif(n_spots), row.names = colnames(spatial_counts))
  spatial <- spacexr::SpatialRNA(coords, Matrix::Matrix(spatial_counts, sparse = TRUE))

  rctd <- spacexr::create.RCTD(
    spatial, reference,
    max_cores = 1, test_mode = FALSE, UMI_min = 10, counts_MIN = 1,
    UMI_min_sigma = 10, CELL_MIN_INSTANCE = 5,
    CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20
  )
  rctd <- spacexr::fitBulk(rctd)
  rctd <- spacexr:::choose_sigma_c(rctd)
  vanilla <- spacexr:::fitPixels(rctd, doublet_mode = "doublet")
  fast <- run.RCTD.fast.doublet(rctd, max_cores = 2, chunk_size = 3)

  expect_named(fast@results, c("results_df", "weights", "weights_doublet"))
  expect_named(fast@results$results_df, names(vanilla@results$results_df))
  expect_identical(rownames(fast@results$results_df), rownames(vanilla@results$results_df))
  expect_identical(colnames(fast@results$weights_doublet), colnames(vanilla@results$weights_doublet))
  expect_identical(as.character(fast@results$results_df$spot_class), as.character(vanilla@results$results_df$spot_class))
  expect_identical(as.character(fast@results$results_df$first_type), as.character(vanilla@results$results_df$first_type))
  expect_identical(as.character(fast@results$results_df$second_type), as.character(vanilla@results$results_df$second_type))
  expect_lt(max(abs(as.matrix(fast@results$weights) - as.matrix(vanilla@results$weights))), 1e-4)
  expect_lt(max(abs(as.matrix(fast@results$weights_doublet) - as.matrix(vanilla@results$weights_doublet))), 1e-4)
})
