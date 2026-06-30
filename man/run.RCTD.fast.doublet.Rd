\name{run.RCTD.fast.doublet}
\alias{run.RCTD.fast.doublet}
\title{Fast drop-in replacement for RCTD doublet pixel fitting}
\usage{
run.RCTD.fast.doublet(
  rctd,
  max_cores = 8,
  chunk_size = NULL,
  exact = TRUE,
  return_diagnostics = FALSE,
  validate_inputs = TRUE,
  progress = interactive()
)
}
\arguments{
\item{rctd}{A prepared \code{spacexr::RCTD} object.}
\item{max_cores}{Number of C++ worker threads.}
\item{chunk_size}{Number of cells per dense chunk, \code{NULL}, or \code{"auto"}.}
\item{exact}{Reserved for future approximate modes. Must be \code{TRUE}.}
\item{return_diagnostics}{Whether to include score matrices and singlet scores.}
\item{validate_inputs}{Whether to validate the RCTD object before running.}
\item{progress}{Whether to print chunk-level progress with elapsed time, ETA,
and processing rate.}
}
\value{The input RCTD object with \code{@results} populated.}
\description{
Preserves upstream RCTD preprocessing and output schema while moving all
per-cell doublet-mode fitting and classification into C++.
}
