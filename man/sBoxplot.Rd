\name{sBoxplot}
\alias{sBoxplot}
\title{Season boxplots}
\usage{
  sBoxplot(data, seasons, temporal = TRUE,
    months = c("rectangles", "lines"), cluster = TRUE,
    multi = FALSE, samescale = TRUE)
}
\arguments{
  \item{data}{The original data on which the clustering was
  made (see \code{\link{gap}}).}

  \item{seasons}{The result of a season clustering: A
  vector of integers (from \code{1:k}) indicating the
  cluster to which each day is allocated.}

  \item{temporal}{Logical. If \code{TRUE}, produces
  boxplots along the year (X-axis); if \code{FALSE},
  produces boxplots for each cluster using their index.}

  \item{months}{Draws the months with background rectangles
  (\code{rectangle}) or dotted lines (\code{lines}).}

  \item{cluster}{Logical. Indicates the cluster index above
  the graph.}

  \item{multi}{Logical. Allows for comparison between
  several clusterings, by displaying them side by side. If
  \code{yes}, requires a list of \code{data} and
  \code{seasons}, corresponding to each clustering.}

  \item{samescale}{Logical. In case of comparison, use the
  same scale for common variables.}
}
\description{
  Season boxplots.
}
\examples{
data(caribou)
set.seed(1)
seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
sBoxplot(caribou$window, sSimple(seasons))
}
\author{
  Mathieu Basille \email{basille@ase-research.org}
}

