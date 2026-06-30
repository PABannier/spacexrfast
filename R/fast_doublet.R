spacexr_ns <- function() {
  asNamespace("spacexr")
}

spacexr_fun <- function(name) {
  get(name, envir = spacexr_ns(), inherits = FALSE)
}

validate_fast_rctd_input <- function(rctd) {
  if (!methods::is(rctd, "RCTD")) {
    stop("run.RCTD.fast.doublet: `rctd` must be a spacexr RCTD object.", call. = FALSE)
  }
  required_internal <- c("gene_list_reg", "Q_mat", "X_vals")
  missing_internal <- setdiff(required_internal, names(rctd@internal_vars))
  if (length(missing_internal) > 0) {
    stop(
      "run.RCTD.fast.doublet: RCTD object is not prepared; missing internal_vars: ",
      paste(missing_internal, collapse = ", "),
      ". Run fitBulk()/choose_sigma_c() or run.RCTD() preprocessing first.",
      call. = FALSE
    )
  }
  if (is.null(rctd@cell_type_info$renorm)) {
    stop("run.RCTD.fast.doublet: missing `cell_type_info$renorm`; run fitBulk() first.", call. = FALSE)
  }
  if (length(rctd@cell_type_info$renorm[[2]]) < 2) {
    stop("run.RCTD.fast.doublet: doublet mode requires at least two cell types.", call. = FALSE)
  }
  invisible(TRUE)
}

estimate_chunk_size <- function(n_genes, target_gb = 2) {
  max(1L, floor((target_gb * 1024^3) / (n_genes * 8)))
}

normalize_chunk_size <- function(chunk_size, n_genes, n_cells) {
  if (is.null(chunk_size) || identical(chunk_size, "auto")) {
    chunk_size <- estimate_chunk_size(n_genes)
  }
  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1L) {
    stop("run.RCTD.fast.doublet: `chunk_size` must be a positive integer, NULL, or 'auto'.", call. = FALSE)
  }
  min(chunk_size, n_cells)
}

extract_likelihood_tables <- function(rctd) {
  set_likelihood_vars <- spacexr_fun("set_likelihood_vars")
  set_likelihood_vars(rctd@internal_vars$Q_mat, rctd@internal_vars$X_vals)
  ns <- spacexr_ns()
  list(
    Q_mat = get("Q_mat", envir = ns),
    SQ_mat = get("SQ_mat", envir = ns),
    X_vals = get("X_vals", envir = ns),
    K_val = get("K_val", envir = ns)
  )
}

extract_class_ids <- function(class_df, cell_type_names) {
  if (is.null(class_df)) {
    return(rep.int(-1L, length(cell_type_names)))
  }
  if (!("class" %in% colnames(class_df))) {
    stop("run.RCTD.fast.doublet: `class_df` must contain a `class` column.", call. = FALSE)
  }
  missing_types <- setdiff(cell_type_names, rownames(class_df))
  if (length(missing_types) > 0) {
    stop(
      "run.RCTD.fast.doublet: `class_df` is missing rows for cell types: ",
      paste(missing_types, collapse = ", "),
      call. = FALSE
    )
  }
  classes <- as.factor(class_df[cell_type_names, "class"])
  as.integer(classes)
}

run_full_fit_chunk <- function(ref_profiles, beads, nUMI, type_names, constrain, min_change) {
  decompose_full <- spacexr_fun("decompose_full")
  n_cells <- nrow(beads)
  n_types <- length(type_names)
  all_weights <- matrix(NA_real_, nrow = n_cells, ncol = n_types)
  conv_all <- logical(n_cells)
  colnames(all_weights) <- type_names

  for (i in seq_len(n_cells)) {
    cell_type_profiles <- data.matrix(ref_profiles * nUMI[[i]])
    results_all <- decompose_full(
      cell_type_profiles,
      nUMI[[i]],
      beads[i, ],
      constrain = constrain,
      verbose = FALSE,
      MIN_CHANGE = min_change
    )
    all_weights[i, ] <- results_all$weights[type_names]
    conv_all[[i]] <- isTRUE(results_all$converged)
  }

  list(all_weights = all_weights, conv_all = conv_all)
}

refit_selected_doublets_chunk <- function(fast, ref_profiles, beads, nUMI, type_names, constrain, min_change) {
  decompose_sparse <- spacexr_fun("decompose_sparse")
  for (i in seq_len(nrow(beads))) {
    first <- type_names[[fast$first_type[[i]]]]
    second <- type_names[[fast$second_type[[i]]]]
    cell_type_profiles <- data.matrix(ref_profiles * nUMI[[i]])
    result <- decompose_sparse(
      cell_type_profiles,
      nUMI[[i]],
      beads[i, ],
      type1 = first,
      type2 = second,
      constrain = constrain,
      MIN.CHANGE = min_change
    )
    fast$first_weight[[i]] <- unname(result$weights[[1]])
    fast$second_weight[[i]] <- unname(result$weights[[2]])
    fast$conv_doublet[[i]] <- isTRUE(result$converged)
  }
  fast
}

assemble_fast_results <- function(rctd, chunk_results, return_diagnostics) {
  cell_type_names <- rctd@cell_type_info$renorm[[2]]
  barcodes <- colnames(rctd@spatialRNA@counts)
  n_cells <- length(barcodes)
  n_types <- length(cell_type_names)

  all_weights <- matrix(0, nrow = n_cells, ncol = n_types)
  weights_doublet <- matrix(0, nrow = n_cells, ncol = 2)
  spot_levels <- c("reject", "singlet", "doublet_certain", "doublet_uncertain")
  empty_cell_types <- factor(character(n_cells), levels = cell_type_names)
  results_df <- data.frame(
    spot_class = factor(character(n_cells), levels = spot_levels),
    first_type = empty_cell_types,
    second_type = empty_cell_types,
    first_class = logical(n_cells),
    second_class = logical(n_cells),
    min_score = numeric(n_cells),
    singlet_score = numeric(n_cells),
    conv_all = logical(n_cells),
    conv_doublet = logical(n_cells)
  )

  diagnostics <- if (return_diagnostics) {
    list(score_mat = vector("list", n_cells), singlet_scores = vector("list", n_cells))
  } else {
    NULL
  }

  for (chunk in chunk_results) {
    idx <- chunk$indices
    all_weights[idx, ] <- chunk$all_weights
    weights_doublet[idx, ] <- cbind(chunk$first_weight, chunk$second_weight)
    results_df$spot_class[idx] <- factor(chunk$spot_class, levels = spot_levels)
    results_df$first_type[idx] <- factor(cell_type_names[chunk$first_type], levels = cell_type_names)
    results_df$second_type[idx] <- factor(cell_type_names[chunk$second_type], levels = cell_type_names)
    results_df$first_class[idx] <- chunk$first_class
    results_df$second_class[idx] <- chunk$second_class
    results_df$min_score[idx] <- chunk$min_score
    results_df$singlet_score[idx] <- chunk$singlet_score
    results_df$conv_all[idx] <- chunk$conv_all
    results_df$conv_doublet[idx] <- chunk$conv_doublet

    if (return_diagnostics) {
      diagnostics$score_mat[idx] <- chunk$score_mat
      diagnostics$singlet_scores[idx] <- chunk$singlet_scores
    }
  }

  rownames(all_weights) <- barcodes
  colnames(all_weights) <- cell_type_names
  rownames(weights_doublet) <- barcodes
  colnames(weights_doublet) <- c("first_type", "second_type")
  rownames(results_df) <- barcodes

  rctd@results <- list(
    results_df = results_df,
    weights = Matrix::Matrix(all_weights, sparse = FALSE),
    weights_doublet = Matrix::Matrix(weights_doublet, sparse = FALSE)
  )
  if (return_diagnostics) {
    rctd@results$score_mat <- diagnostics$score_mat
    rctd@results$singlet_scores <- diagnostics$singlet_scores
  }
  rctd
}

#' Fast drop-in replacement for RCTD doublet pixel fitting
#'
#' This function assumes upstream RCTD preprocessing has already been run. It
#' preserves the RCTD doublet result schema while using chunked count extraction
#' and a C++ backend for singlet scoring, candidate-pair scoring, selected-pair
#' refitting, and classification.
#'
#' @param rctd A prepared `spacexr::RCTD` object.
#' @param max_cores Number of worker threads for the C++ doublet backend.
#' @param chunk_size Number of cells per dense chunk, `NULL`, or `"auto"`.
#' @param exact Reserved for future approximate modes. Must be `TRUE`.
#' @param return_diagnostics If `TRUE`, include per-cell pair score matrices and
#'   singlet scores.
#' @param validate_inputs Whether to validate the RCTD object before running.
#' @return The input RCTD object with `@results` populated.
run.RCTD.fast.doublet <- function(rctd,
                                  max_cores = 8,
                                  chunk_size = NULL,
                                  exact = TRUE,
                                  return_diagnostics = FALSE,
                                  validate_inputs = TRUE) {
  if (!isTRUE(exact)) {
    stop("run.RCTD.fast.doublet: approximate mode is not implemented.", call. = FALSE)
  }
  if (validate_inputs) {
    validate_fast_rctd_input(rctd)
  }

  rctd@internal_vars$cell_types_assigned <- TRUE
  rctd@config$RCTDmode <- "doublet"

  cell_type_info <- rctd@cell_type_info$renorm
  gene_list <- rctd@internal_vars$gene_list_reg
  counts <- rctd@spatialRNA@counts[gene_list, , drop = FALSE]
  nUMI <- rctd@spatialRNA@nUMI
  cell_type_names <- cell_type_info[[2]]
  ref_profiles <- data.matrix(cell_type_info[[1]][gene_list, cell_type_names, drop = FALSE])
  class_ids <- extract_class_ids(rctd@internal_vars$class_df, cell_type_names)
  likelihood <- extract_likelihood_tables(rctd)

  n_cells <- ncol(counts)
  chunk_size <- normalize_chunk_size(chunk_size, length(gene_list), n_cells)
  max_cores <- as.integer(max_cores)
  if (is.na(max_cores) || max_cores < 1L) {
    stop("run.RCTD.fast.doublet: `max_cores` must be a positive integer.", call. = FALSE)
  }
  min_change <- rctd@config$MIN_CHANGE_REG %||% 0.001
  confidence_threshold <- rctd@config$CONFIDENCE_THRESHOLD %||% 10
  doublet_threshold <- rctd@config$DOUBLET_THRESHOLD %||% 25
  constrain <- FALSE

  starts <- seq.int(1L, n_cells, by = chunk_size)
  chunk_results <- vector("list", length(starts))

  for (chunk_idx in seq_along(starts)) {
    start <- starts[[chunk_idx]]
    end <- min(start + chunk_size - 1L, n_cells)
    idx <- seq.int(start, end)
    beads <- t(as.matrix(counts[, idx, drop = FALSE]))
    storage.mode(beads) <- "double"

    full_fit <- run_full_fit_chunk(
      ref_profiles = ref_profiles,
      beads = beads,
      nUMI = nUMI[idx],
      type_names = cell_type_names,
      constrain = constrain,
      min_change = min_change
    )

    fast <- fast_doublet_chunk_cpp(
      beads = beads,
      nUMI = as.numeric(nUMI[idx]),
      ref_profiles = ref_profiles,
      all_weights = full_fit$all_weights,
      conv_all = full_fit$conv_all,
      class_ids = as.integer(class_ids),
      Q_mat = likelihood$Q_mat,
      SQ_mat = likelihood$SQ_mat,
      X_vals = likelihood$X_vals,
      K_val = as.integer(likelihood$K_val),
      constrain = constrain,
      min_change = min_change,
      confidence_threshold = confidence_threshold,
      doublet_threshold = doublet_threshold,
      return_diagnostics = return_diagnostics,
      max_cores = max_cores
    )
    fast <- refit_selected_doublets_chunk(
      fast = fast,
      ref_profiles = ref_profiles,
      beads = beads,
      nUMI = nUMI[idx],
      type_names = cell_type_names,
      constrain = constrain,
      min_change = min_change
    )
    fast$indices <- idx
    fast$all_weights <- full_fit$all_weights
    fast$conv_all <- full_fit$conv_all
    chunk_results[[chunk_idx]] <- fast
  }

  assemble_fast_results(rctd, chunk_results, return_diagnostics)
}

#' Convenience wrapper for fast RCTD modes
#'
#' @param rctd A prepared `spacexr::RCTD` object.
#' @param doublet_mode Currently only `"doublet"` is supported.
#' @param ... Passed to [run.RCTD.fast.doublet()].
#' @return The updated RCTD object.
run.RCTD.fast <- function(rctd, doublet_mode = "doublet", ...) {
  if (!identical(doublet_mode, "doublet")) {
    stop("run.RCTD.fast: only doublet_mode = 'doublet' is implemented.", call. = FALSE)
  }
  run.RCTD.fast.doublet(rctd, ...)
}

#' Profile vanilla spacexr doublet fitting
#'
#' Runs upstream `fitPixels(doublet_mode = "doublet")` and returns elapsed time
#' with the resulting RCTD object. This is intentionally small and dependency
#' free so it can be used on validation slices.
#'
#' @param rctd A prepared `spacexr::RCTD` object.
#' @return A list with `elapsed_seconds` and `rctd`.
profile_vanilla_doublet <- function(rctd) {
  fitPixels <- spacexr_fun("fitPixels")
  elapsed <- system.time({
    result <- fitPixels(rctd, doublet_mode = "doublet")
  })
  list(elapsed_seconds = unname(elapsed[["elapsed"]]), rctd = result)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
