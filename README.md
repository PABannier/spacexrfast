# spacexrfast

`spacexrfast` provides a faster drop-in doublet-mode backend for
`spacexr` RCTD. It keeps upstream RCTD preprocessing, all-type fitting, result
assembly, and final selected-pair refitting, while moving the hot doublet
candidate pair scoring path into C++.

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
```

The returned object uses the standard RCTD doublet schema:

- `rctd@results$results_df`
- `rctd@results$weights`
- `rctd@results$weights_doublet`

## Install

From this checkout:

```sh
R CMD INSTALL --preclean .
```

## Run Tests

Run the unit and integration tests:

```sh
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
```

Run a full package check:

```sh
R CMD build .
R CMD check --no-manual --no-vignettes spacexrfast_0.1.0.tar.gz
```

## Xenium Benchmark

Download a Xenium example from the 10x Genomics website, then point the
benchmark at its `cell_feature_matrix.h5`.

Useful 10x pages:

- Renal cell carcinoma Xenium dataset:
  <https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma>
- Xenium file format documentation, including `cell_feature_matrix.h5`:
  <https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html>

After downloading and extracting an output bundle, run:

```sh
XENIUM_MATRIX="/path/to/cell_feature_matrix.h5" \
XENIUM_SMOKE_CELLS=100 \
XENIUM_SMOKE_GENES=200 \
XENIUM_SMOKE_TYPES=20 \
XENIUM_SMOKE_CORES=4 \
Rscript benchmarks/xenium_smoke.R
```

The benchmark builds a small pseudo-reference from the downloaded Xenium count
matrix, runs vanilla RCTD doublet fitting and `run.RCTD.fast.doublet()`, and
prints runtime plus label agreement.

Example result on the local NSCLC Xenium slide in `~/Documents/10x-examples`:

```text
cells,100
genes,200
types,20
vanilla_elapsed_seconds,21.546
fast_elapsed_seconds,1.399
spot_class_equal,TRUE
first_type_equal,TRUE
second_type_equal,TRUE
```

That slice showed a 15.4x speedup with matching `spot_class`, `first_type`, and
`second_type`.

## Notes

- The exact-mode API rejects approximate mode for now.
- Candidate generation still uses upstream RCTD all-type fitting.
- The final reported selected-pair weights are refit through upstream RCTD to
  preserve `weights_doublet` compatibility.
- The C++ backend parallelizes over cells and does not call the R API from
  worker threads.
