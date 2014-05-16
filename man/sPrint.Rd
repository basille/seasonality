\name{sPrint}
\alias{sPrint}
\title{Print seasons}
\usage{
  sPrint(seasons, ndays = FALSE)
}
\arguments{
  \item{seasons}{The result of a season clustering: A
  vector of integers (from \code{1:k}) indicating the
  cluster to which each day is allocated.}

  \item{ndays}{Logical. Returns the rank of the days at
  which a new season starts.}
}
\value{
  A vector indicating the dates at which a new season
  starts.
}
\description{
  Print a sequence of seasons in a friendly way.
}
\examples{
data(caribou)
set.seed(1)
seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
sPrint(seasons)
}
\author{
  Mathieu Basille \email{basille@ase-research.org}
}

