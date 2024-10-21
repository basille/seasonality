## caribou
##
##' Caribou data set used in Basille et al. (2012).
##'
##' @format ## `caribou`
##' A list with 4 elements:
##' \describe{
##'   \item{ind}{A 91×2 data frame, providing the code of the individual (`ind`)
##'     and the year (`year`).}
##'   \item{move}{A list with 9 data frames of environmental variables, each with
##'     365 rows (julian day) and 91 columns (individual-year), range-standardised
##'     by individual-year.}
##'   \item{window}{A 365×9 data frame, providing the value of environmental
##'     variables over a 15-day moving window for all caribou combined.}
##'   \item{bs}{A list of 100 repetitions, providing season numbers as integers
##'     from 0 to i for each day of the year.}
##' }
##' @source Basille M., Fortin D., Dussault C., Ouellet J.-P., Courtois
##'     R. Ecologically based definition of seasons clarifies predator-prey
##'     interactions. Ecography, early view. DOI:
##'     10.1111/j.1600-0587.2011.07367.x
##' @examples
##' ## Load the caribou dataset:
##' data("caribou")
##'
##' ## `caribou$window` can be retrieved from `caribou$move`:
##' ##
##' ## We first compute the weights for each individual-year:
##' caribou$ind$weights <- rep(1/table(caribou$ind$id),
##'     times = table(caribou$ind$id))
##' ## We then compute the weighted mean for each environmental variable:
##' cari_win <- data.frame(lapply(caribou$move, apply, 1, weighted.mean,
##'     w = caribou$ind$weights, na.rm = TRUE))
##' ## And finally range-standardise the resulting data:
##' cari_win <- data.frame(scale(cari_win,
##'     center = apply(cari_win, 2, min, na.rm = TRUE),
##'     scale = apply(cari_win, 2, \(x) diff(range(x, na.rm = TRUE)))),
##'     row.names = as.character(1:365))
##' head(cari_win)
##' summary(cari_win)
##' all.equal(cari_win, caribou$window)
##'
##' ## We then compute the gap statistic for 1–10 clusters:
##' (caribou$gap <- gap(caribou$window)) # Different results from Basille et al.
##' plot(caribou$gap)
##' ## And compute seasons with the K-means for 8 clusters:
##' set.seed(1)
##' caribou$seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
##' sPrint(caribou$seasons)
##'
##' ## Bootstrap approach:
##' caribou$bs <- bsSeasons(data = caribou$move, ind = caribou$ind,
##'     nclust = 8)
##' ## Compute the bootstrap weights, and identify 'true' seasons:
##' caribou$bsweights <- bsWeights(caribou$bs)
##' caribou$seasonsbs <- bsCriterion(caribou$seasons, caribou$bsweights,
##'     threshold = .9)
##' sPrint(caribou$seasonsbs)
##' ## Check visually:
##' bsPlot(caribou$seasonsbs, caribou$seasons, caribou$bsweights,
##'     title = "Bootstrap on Caribou seasons")
##'
##' ## Simplify to get the final seasons:
##' sPrint(sFormat(sSimple(caribou$seasonsbs)))
##' ## We end up with the following seasons:
##' ## - December 28–January 2 (winter 1)
##' ## - January 3–April 15 (winter 2)
##' ## - April 16–May 21 (spring dispersal)
##' ## - May 22–May 27 (pre-calving)
##' ## - May 28–September 17 (calving and summer)
##' ## - September 18–December 27 (autumn)
##' ##
##' ## We can visualize the characteristics of the final seasons:
##' sBoxplot(caribou$window, sFormat(sSimple(caribou$seasonsbs)))
"caribou"
