\name{sFormat}
\alias{sFormat}
\title{Reorder the seasons.}
\usage{
  sFormat(seasons)
}
\arguments{
  \item{seasons}{The result of a season clustering: A
  vector of integers (from \code{1:k}) indicating the
  cluster to which each day is allocated.}
}
\value{
  A vector of the same length as \code{seasons}, with the
  index of the clusters reordered.
}
\description{
  Reorder the seasons as a succession of unique numbers,
  from 1 to the last season (useful in case of duplicated
  clusters, as duplicates get a new index).
}
\examples{
data(caribou)
set.seed(1)
seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
sPrint(seasons)
sPrint(sFormat(seasons))
}
\author{
  Mathieu Basille \email{basille@ase-research.org}
}

