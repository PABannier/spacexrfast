# spacexrfast

[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: R + C++](https://img.shields.io/badge/built%20with-R%20%2B%20C%2B%2B-198ce7.svg)](https://www.rcpp.org/)
[![Version](https://img.shields.io/badge/version-0.1.0-brightgreen.svg)](DESCRIPTION)

**A drop-in, exact-preserving doublet-mode backend for [spacexr](https://github.com/dmcable/spacexr) RCTD — much faster on Xenium slides, with identical cell-type calls.**

RCTD doublet mode is the accuracy gold standard for cell-type deconvolution of spatial
transcriptomics, but it was never built for million-cell Xenium slides — a single slide can
take ~10 hours on 8 cores. `spacexrfast` keeps RCTD's preprocessing, statistics,
and output schema, and moves the per-cell doublet fitting path into parallel C++.

```r
library(spacexr)
library(spacexrfast)

# Standard RCTD preprocessing — unchanged
rctd <- create.RCTD(puck, reference, max_cores = 8)
rctd <- fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)

# Swap only the slow doublet pixel fitting:
rctd <- run.RCTD.fast.doublet(rctd, max_cores = 8, chunk_size = "auto", exact = TRUE)
```

The returned object is a normal RCTD object with the standard doublet schema —
`rctd@results$results_df`, `rctd@results$weights`, and `rctd@results$weights_doublet` —
so it slots straight into any existing downstream pipeline.

## Why spacexrfast?

Vanilla RCTD doublet mode is structurally expensive. For **every** cell it runs a full
all-type decomposition, thresholds candidates, scores each candidate type as a singlet,
fits **every** candidate pair with 25 IRWLS iterations, applies the doublet/singlet/reject
logic, then refits the selected pair with 50 more iterations. With ~20 reference types that
is dozens of two-type optimizations per cell, each routed through R and `quadprog`. Across
millions of cells, that overhead dominates.

`spacexrfast` attacks exactly that bottleneck — and nothing else:

- **Exact-preserving by design.** Categorical outputs (`spot_class`, `first_type`, `second_type`,
  `first_class`, `second_class`) match upstream RCTD bit-for-bit on the validation suite. The
  candidate threshold (`all_weights > 0.01`), fallback logic, tie-breaking, factor levels, and
  confidence/doublet thresholds are all reproduced exactly.
- **Hot path in C++.** Likelihood-table interpolation, all-type IRWLS, two-type IRWLS, singlet
  scoring, the candidate-pair loop, the final selected-pair refit, and per-cell classification all
  run as compiled code — no R or `quadprog` calls in the inner loop.
- **Parallel over cells.** The C++ backend fans out across `max_cores` threads with no R API
  access from workers, so it scales cleanly on multi-core machines.
- **Active-set QP solvers.** The repeated WLS QP subproblems that RCTD hands to `quadprog` are
  replaced by thread-safe C++ active-set solvers for all-type and two-type fits.
- **Adaptive chunking.** Cells are processed in dense gene × chunk blocks sized to a memory
  budget (`chunk_size = "auto"`), keeping peak RAM bounded even on million-cell slides.
- **Drop-in API, drop-in schema.** Keep upstream gene filtering, platform normalization, sigma
  selection, and S4 assembly. Only the pixel-fitting step changes.

## Quick Start

**Prerequisites:** R (≥ 4.0), a C++ compiler with `GNU make`, and
[`spacexr`](https://github.com/dmcable/spacexr) installed.

Install from this checkout:

```sh
R CMD INSTALL --preclean .
```

Then run the fast doublet backend on a prepared RCTD object:

```r
library(spacexr)
library(spacexrfast)

rctd <- fitBulk(rctd)
rctd <- spacexr:::choose_sigma_c(rctd)

rctd <- run.RCTD.fast.doublet(
  rctd,
  max_cores = 8,
  chunk_size = "auto",
  exact = TRUE
)

head(rctd@results$results_df)
#>      spot_class first_type second_type ... min_score singlet_score conv_doublet
#> s1      singlet      typeA       typeB ...    12.43         38.91         TRUE
#> s2  doublet_cert      typeC       typeA ...     8.10         11.02         TRUE
```

## Usage

### `run.RCTD.fast.doublet()`

The main entry point. Assumes upstream RCTD preprocessing (`fitBulk()` + sigma selection) has
already run.

```r
run.RCTD.fast.doublet(
  rctd,                       # a prepared spacexr::RCTD object
  max_cores = 8,              # worker threads for the C++ backend
  chunk_size = NULL,          # cells per dense chunk; NULL or "auto" picks adaptively
  exact = TRUE,               # exact mode (approximate mode is not yet implemented)
  return_diagnostics = FALSE, # attach per-cell pair-score matrices + singlet scores
  validate_inputs = TRUE,     # validate the RCTD object before running
  progress = interactive()    # print chunk progress, elapsed time, ETA, cells/s
)
```

For long Xenium runs, force progress output even in batch jobs:

```r
rctd <- run.RCTD.fast.doublet(
  rctd,
  max_cores = 4,
  chunk_size = 10000,
  progress = TRUE
)
```

Example progress output:

```text
run.RCTD.fast.doublet: 465524 cells, 360 genes, 21 types, 47 chunks, max_cores=4
doublet chunks [================>...........] 27/47 cells 270000/465524 elapsed 03:35 ETA 02:39 rate 1256 cells/s
```

### `run.RCTD.fast()`

Convenience wrapper mirroring upstream's `run.RCTD()` signature:

```r
rctd <- run.RCTD.fast(rctd, doublet_mode = "doublet", max_cores = 8)
```

### `profile_vanilla_doublet()`

Times upstream `fitPixels(doublet_mode = "doublet")` on a slice so you can measure the speedup
on your own data:

```r
baseline <- profile_vanilla_doublet(rctd)
baseline$elapsed_seconds
```

## Benchmarks

On a 100-cell Xenium slice (NSCLC slide, 200 genes, 20 reference types) the fast backend produced
**identical** `spot_class`, `first_type`, and `second_type` while running **85× faster**:

| Backend                    | Cells | Genes | Types | Elapsed (s) | Labels match |
|----------------------------|-------|-------|-------|-------------|--------------|
| Vanilla RCTD (`fitPixels`) | 100   | 200   | 20    | 22.612      | —            |
| **`run.RCTD.fast.doublet`**| 100   | 200   | 20    | **0.266**   | ✓ exact      |

> Speedup scales with the number of candidate pairs per cell (more reference types → bigger win).
> Reproduce on your own Xenium data with the smoke benchmark below.

### Run the Xenium benchmark yourself

Download a Xenium example bundle from 10x Genomics and point the benchmark at its
`cell_feature_matrix.h5`:

- [Renal cell carcinoma Xenium dataset](https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma)
- [Xenium file-format documentation](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html)

```sh
XENIUM_MATRIX="/path/to/cell_feature_matrix.h5" \
XENIUM_SMOKE_CELLS=100 \
XENIUM_SMOKE_GENES=200 \
XENIUM_SMOKE_TYPES=20 \
XENIUM_SMOKE_CORES=4 \
Rscript benchmarks/xenium_smoke.R
```

The benchmark builds a small pseudo-reference from the count matrix, runs both backends, and
prints runtime plus label agreement.

The full local kidney Xenium query (`465,524` cells, 21 pseudo-reference types, 4 cores) completed
the `run.RCTD.fast.doublet()` phase in `09:02` with the built-in chunk progress bar enabled.

To verify thread scaling on a prepared RCTD object that has already run `fitBulk()` and
`choose_sigma_c()`, save it with `saveRDS()` and run:

```sh
THREAD_SCALING_CORES=1,4,8 \
THREAD_SCALING_CHUNK_SIZE=auto \
Rscript benchmarks/thread_scaling.R /path/to/prepared_rctd_after_sigma.rds
```

This prints elapsed time and sampled `ps -o nlwp` thread counts for each `max_cores` setting.

## How it works

`spacexrfast` splits responsibility cleanly between R and C++:

- **R owns** RCTD object creation, gene filtering, platform normalization, sigma selection, dense
  chunk extraction, and S4 result assembly.
- **C++ owns** likelihood interpolation, all-type IRWLS for candidate generation, two-type IRWLS,
  singlet scoring, the candidate-pair loop, the final selected-pair refit, and per-cell
  classification.

Cells are streamed in adaptively-sized chunks; each chunk is materialized dense over the RCTD
regression gene list and handed to a thread pool that processes cells independently. No worker
ever touches the R API, which keeps the backend both fast and safe.

See [`docs/profile_report.md`](docs/profile_report.md) for profiling notes and the known
exactness boundary.

## Testing

Run the unit and integration tests:

```sh
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
```

Run a full package check:

```sh
R CMD build .
R CMD check --no-manual --no-vignettes spacexrfast_0.1.0.tar.gz
```

The test suite compares categorical outputs and selected types against upstream
`process_bead_doublet()` and `fitPixels()` on synthetic fixtures.

## Notes & limitations

- **Exact mode only.** The `exact = TRUE` API rejects approximate mode for now.
- **Weight tolerance.** Categorical labels and selected types match exactly in the validation
  suite. Weight matrices are tested against upstream RCTD with max absolute difference below
  `1e-4`.

## License

Distributed under the [GPL-3](https://www.gnu.org/licenses/gpl-3.0) license, matching upstream
`spacexr`.

## Acknowledgments

Built on top of [spacexr](https://github.com/dmcable/spacexr) (Cable et al., *Nature
Biotechnology* 2022). `spacexrfast` reuses spacexr's preprocessing and statistical model
wholesale — it only accelerates the doublet-mode inner loop.
