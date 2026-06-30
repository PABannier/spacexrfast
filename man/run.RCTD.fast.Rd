\name{run.RCTD.fast}
\alias{run.RCTD.fast}
\title{Convenience wrapper for fast RCTD modes}
\usage{
run.RCTD.fast(rctd, doublet_mode = "doublet", ...)
}
\arguments{
\item{rctd}{A prepared \code{spacexr::RCTD} object.}
\item{doublet_mode}{Currently only \code{"doublet"} is supported.}
\item{...}{Passed to \code{\link{run.RCTD.fast.doublet}}.}
}
\value{The updated RCTD object.}
\description{
Dispatches to the fast doublet backend.
}
