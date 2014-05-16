\name{sPlot}
\alias{sPlot}
\title{Plot the seasons}
\usage{
  sPlot(seasons, add.lines = FALSE, months = FALSE,
    main = "Seasons", ylab = substitute(seasons), ...)
}
\arguments{
  \item{seasons}{The result of a season clustering: A
  vector of integers (from \code{1:k}) indicating the
  cluster to which each day is allocated.}

  \item{add.lines}{Logical. Adds dotted lines delineating
  the seasons.}

  \item{months}{Logical. Draws monthtly delineations.}

  \item{main}{An overall title for the plot.}

  \item{ylab}{A title for the y axis.}

  \item{...}{Further arguments passed to the \code{lines}
  call.}
}
\description{
  Plot the seasons.
}
\examples{
data(caribou)
set.seed(1)
seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
sPlot(sSimple(seasons))
}
\author{
  Mathieu Basille \email{basille@ase-research.org}
}

