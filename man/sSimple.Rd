\name{sSimple}
\alias{sSimple}
\title{Simplify the seasons.}
\usage{
  sSimple(clust, win = 3, tol = 1)
}
\arguments{
  \item{clust}{The result of a season clustering: A vector
  of integers (from \code{1:k}) indicating the cluster to
  which each day is allocated.}

  \item{win}{The length of the moving window.}

  \item{tol}{The tolerance to be used.}
}
\value{
  A vector of the same length as \code{seasons}.
}
\description{
  Simplify the seasons after the initial clustering, by
  removing the smallest seasons.
}
\details{
  The function works on a moving window of length (current
  day + \code{win}) days. For a given day, if all other
  days (with a tolerance of \code{tol} days) have the same
  value as the focus day, the day is kept as is; otherwise,
  the day takes the value of the day before.
}
\examples{
data(caribou)
set.seed(1)
seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
sPrint(seasons)
sPrint(sSimple(seasons))
}
\author{
  Mathieu Basille \email{basille@ase-research.org}
}

