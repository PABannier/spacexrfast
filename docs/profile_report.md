# Fast RCTD Doublet Profile Notes

This package implements the first production-oriented fast path from
`IMPLEMENTATION_PLAN.md`.

Implemented:

- `run.RCTD.fast.doublet()` drop-in API.
- Adaptive chunked dense extraction over the RCTD regression gene list.
- C++ likelihood interpolation, all-type IRWLS for candidate generation, one-
  and two-type IRWLS, candidate pair loop, singlet scoring, final selected-pair
  refit, and classification.
- Parallel processing over cells via standard C++ threads.
- Optional compact diagnostics through `return_diagnostics`.

Validation:

- `tests/testthat/test-fast-doublet.R` compares categorical outputs, selected
  types, all-type weights, and doublet weights against upstream
  `process_bead_doublet()` and `fitPixels()`.
- `benchmarks/xenium_smoke.R` samples a local Xenium H5 matrix from
  `~/Documents/10x-examples` and compares vanilla and fast outputs on a small
  smoke slice.
- NSCLC smoke slice, 100 cells x 200 genes x 20 types: vanilla `22.612s`,
  fast `0.266s`, with matching labels and selected types.

Known exactness boundary:

- WLS QP updates use C++ active-set solvers instead of calling `quadprog` in
  workers. Categorical labels and selected types match the validation cases;
  weights are tested with max absolute difference below `1e-4`.
