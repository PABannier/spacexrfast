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
  validate_inputs = TRUE
)
}
\arguments{
\item{rctd}{A prepared \code{spacexr::RCTD} object.}
\item{max_cores}{Number of C++ worker threads.}
\item{chunk_size}{Number of cells per dense chunk, \code{NULL}, or \code{"auto"}.}
\item{exact}{Reserved for future approximate modes. Must be \code{TRUE}.}
\item{return_diagnostics}{Whether to include score matrices and singlet scores.}
\item{validate_inputs}{Whether to validate the RCTD object before running.}
}
\value{The input RCTD object with \code{@results} populated.}
\description{
Preserves upstream RCTD preprocessing and output schema while moving repeated
doublet-mode scoring and selected-pair refitting into C++.
}
