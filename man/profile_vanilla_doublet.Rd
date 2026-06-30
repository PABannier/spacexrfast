\name{profile_vanilla_doublet}
\alias{profile_vanilla_doublet}
\title{Profile vanilla spacexr doublet fitting}
\usage{
profile_vanilla_doublet(rctd)
}
\arguments{
\item{rctd}{A prepared \code{spacexr::RCTD} object.}
}
\value{A list with \code{elapsed_seconds} and the resulting \code{rctd}.}
\description{
Runs upstream \code{fitPixels(doublet_mode = "doublet")} for validation slices.
}
