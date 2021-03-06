\name{bsSeasons}
\alias{bsCriterion}
\alias{bsPlot}
\alias{bsSeasons}
\alias{bsWeights}
\title{Season bootstrap}
\usage{
  bsSeasons(data, ind, dataNA, nclust, iter = 100,
    simplify = FALSE, win = 6, tol = 2)

  bsWeights(bsSeasons)

  bsCriterion(seasons, bsWeights, threshold = 0.75)

  bsPlot(bsSeasons, seasons = NULL, bsWeights, title)
}
\arguments{
  \item{data}{A data frame indicating the initial data on
  which to run the bootstrap.}

  \item{ind}{A individual-year table indicating the name of
  the individual (column \code{id}), repeated as many times
  as monitoring periods.}

  \item{dataNA}{The data to use in case of NAs, since the
  k-means can not deal with NAs.}

  \item{nclust}{The number of clusters to apply to the
  k-means.}

  \item{iter}{The number of iterations of the bootstrap.}

  \item{simplify}{Logical. Whether to simplify the
  resulting seasons.}

  \item{win}{If \code{simplify}, the length of the moving
  window.}

  \item{tol}{If \code{simplify}, the tolerance to be used.}

  \item{bsSeasons}{Bootstrap seasons.}

  \item{seasons}{The result of a season clustering: A
  vector of integers (from \code{1:k}) indicating the
  cluster to which each day is allocated.}

  \item{bsWeights}{The bootsrap weights, as given by
  \code{bsWeights}.}

  \item{threshold}{The weight threshold.}

  \item{title}{A title for the strip.}
}
\value{
  A list of length \code{iter}, each element of which
  giving the clustering of one bootstrap iteration.

  \code{bsSeasons} returns a list of vectors of the same
  length as \code{seasons}, giving the seasons for each
  bootstrap loop.

  \code{bsWeights} returns a vector of the same length as
  \code{seasons}, with the weight of each day.

  \code{bsCriterion} returns a vector of the same length as
  \code{seasons}, with the index of the clusters kept.
}
\description{
  Bootstrap procedures to remove the less important
  seasons.
}
\details{
  The function \code{bsSeasons} samples individual
  animal-years with replacement, and run the K-means
  clustering with the estimated number of clusters for the
  whole data set.

  The weight is then given by the function
  \code{bsWeights}, which gives, for each day of the year
  (from \code{1:365}) the number of changes in the last and
  next two days. This weight is then used by the function
  \code{bsCriterion} to retain only seasons which are
  within a given threshold of weight (based on the
  bootstrap data set).

  \code{bsPlot} plots the result of the bootstrap
  procedure.
}
\examples{
### Load the data
data(caribou)

### Recompute the bootstrap seasons:
\dontrun{
caribou$bs <- bsSeasons(data = caribou$move, ind = caribou$ind,
    dataNA = caribou$window, nclust = 8)}

### Compute the weights, and identify the final seasons:
set.seed(1)
seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
weights <- bsWeights(caribou$bs)
seasonsbs <- bsCriterion(seasons, weights, threshold = .8)

### Visualize the final seasons:
bsPlot(seasonsbs, seasons, weights, title = "Caribou")
}
\author{
  Mathieu Basille

  Mathieu Basille \email{basille@ase-research.org}
}

