## bsSeasons
##
##' Bootstrap procedures to remove the less important seasons.
##'
##' The function \code{bsSeasons} samples individual animal-years with
##' replacement, and run the K-means clustering with the estimated
##' number of clusters for the whole data set.
##'
##' The weight is then given by the function \code{bsWeights}, which
##' gives, for each day of the year (from \code{1:365}) the number of
##' changes in the last and next two days. This weight is then used by
##' the function \code{bsCriterion} to retain only seasons which are
##' within a given threshold of weight (based on the bootstrap data
##' set).
##'
##' \code{bsPlot} plots the result of the bootstrap procedure.
##'
##' @title Season bootstrap
##' @param data A data frame indicating the initial data on which to
##' run the bootstrap.
##' @param ind A individual-year table indicating the name of the
##' individual (column \code{id}), repeated as many times as
##' monitoring periods.
##' @param dataNA The data to use in case of NAs, since the k-means
##' can not deal with NAs.
##' @param nclust The number of clusters to apply to the k-means.
##' @param iter The number of iterations of the bootstrap.
##' @param simplify Logical. Whether to simplify the resulting seasons.
##' @param win If \code{simplify}, the length of the moving window.
##' @param tol If \code{simplify}, the tolerance to be used.
##' @return A list of length \code{iter}, each element of which giving
##' the clustering of one bootstrap iteration.
##' @author Mathieu Basille
##' @return \code{bsSeasons} returns a list of vectors of the same
##' length as \code{seasons}, giving the seasons for each bootstrap
##' loop.
##'
##' \code{bsWeights} returns a vector of the same length as
##' \code{seasons}, with the weight of each day.
##'
##' \code{bsCriterion} returns a vector of the same length as
##' \code{seasons}, with the index of the clusters kept.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ### Load the data
##' data(caribou)
##'
##' ### Recompute the bootstrap seasons:
##' \dontrun{
##' caribou$bs <- bsSeasons(data = caribou$move, ind = caribou$ind,
##'     dataNA = caribou$window, nclust = 8)}
##'
##' ### Compute the weights, and identify the final seasons:
##' set.seed(1)
##' seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
##' weights <- bsWeights(caribou$bs)
##' seasonsbs <- bsCriterion(seasons, weights, threshold = .8)
##'
##' ### Visualize the final seasons:
##' bsPlot(seasonsbs, seasons, weights, title = "Caribou")
bsSeasons <- function(data, ind, dataNA, nclust, iter = 100,
    simplify = FALSE, win = 6, tol = 2) {
    bs <- list()
    datanorm <- bsi <- function(i) {
        set.seed(i)
        rd <- sort(sample(1:nrow(ind), nrow(ind), replace = TRUE))
        names <- ind[rd, ]
        names$weights <- rep(1/table(names$id), times = table(names$id))
        norm <- data.frame(row.names = 1:365)
        for (j in names(data)) norm[, j] <- apply(data[[j]][,
            rd], 1, weighted.mean, w = names$weights, na.rm = TRUE)
        norm <- as.data.frame(scale(norm, center = apply(norm,
            2, min, na.rm = TRUE), scale = apply(norm, 2, function(x) diff(range(x,
            na.rm = TRUE)))))
        if (any(is.na(norm)))
            norm[is.na(norm)] <- data$norm[is.na(norm)]
        set.seed(1)
        seasons <- kmeans(norm, nclust, iter.max = 100)$cluster
        if (simplify)
            return(sFormat(sSimple(seasons, win, tol)))
        else return(sFormat(seasons))
    }
    bs <- lapply(1:iter, bsi)
    names(bs) <- 1:iter
    return(bs)
}


## bsWeights
##
##' @param bsSeasons Bootstrap seasons.
##' @export
##' @rdname bsSeasons
bsWeights <- function(bsSeasons)
{
    changes <- function(x) {
        seasons <- sPrint(x, ndays = TRUE)
        seasons <- seasons[-length(seasons)]
        if (x[1] == x[365])
            seasons <- seasons[-1]
        return(seasons)
    }
    bsw <- unlist(lapply(bsSeasons, changes))
    bsw <- as.numeric(table(factor(bsw, levels = 1:365)))
    bswm2 <- c(bsw[364:365], bsw[1:363])
    bswm1 <- c(bsw[365], bsw[1:364])
    bswp1 <- c(bsw[2:365], bsw[1])
    bswp2 <- c(bsw[3:365], bsw[1:2])
    return(bswm2 + bswm1 + bsw + bswp1 + bswp2)
}


## bsCriterion
##
##' @param seasons The result of a season clustering: A vector of
##' integers (from \code{1:k}) indicating the cluster to which each
##' day is allocated.
##' @param bsWeights The bootsrap weights, as given by
##' \code{bsWeights}.
##' @param threshold The weight threshold.
##' @export
##' @rdname bsSeasons
bsCriterion <- function(seasons, bsWeights, threshold = 0.75)
{
    ch <- sPrint(seasons, ndays = TRUE)
    ch <- ch[-c(1, length(ch))]
    chw <- as.numeric(names(ch)[names(ch) %in% which(bsWeights >
        quantile(bsWeights, threshold))])
    seasw <- rep(1, 365)
    for (i in 1:(length(chw) - 1)) seasw[chw[i]:(chw[i + 1] -
        1)] <- i + 1
    names(seasw) <- 1:365
    return(seasw)
}


## bsPlot
##
##' @param title A title for the strip.
##' @export
##' @rdname bsSeasons
bsPlot <- function(bsSeasons, seasons = NULL, bsWeights, title)
{
    plot(1:365, rep(1, 365), ylim = c(0, 1), type = "n", axes = FALSE,
        ann = FALSE)
    abline(v = c(32, 60, 91, 121, 152, 182, 213, 244, 274, 305,
        335), lty = 2, col = grey(0.5))
    rect(1:365, 0.3, 2:366, 0.95, col = grey((max(bsWeights) -
        bsWeights)/max(bsWeights)), border = NA)
    if (!is.null(seasons)) {
        xori <- sPrint(seasons, ndays = TRUE)
        xori <- xori[-c(1, length(xori))]
        segments(x0 = xori, y0 = 0.05, y1 = 0.7, lwd = 2, lty = 3)
    }
    xs <- sPrint(bsSeasons, ndays = TRUE)
    xs <- xs[-c(1, length(xs))]
    segments(x0 = xs, y0 = 0.05, y1 = 0.7, lwd = 2)
    axis(1, at = c(16, 46, 75.5, 106, 136.5, 167, 197.5, 228.5,
        259, 289.5, 320, 350), labels = c("J", "F", "M", "A",
        "M", "J", "J", "A", "S", "O", "N", "D"), tick = FALSE,
        line = -0.5)
    axis(3, at = c(16, 46, 75.5, 106, 136.5, 167, 197.5, 228.5,
        259, 289.5, 320, 350), labels = c("J", "F", "M", "A",
        "M", "J", "J", "A", "S", "O", "N", "D"), tick = FALSE,
        line = -0.5)
    title(ylab = title, cex.lab = 2, line = 1)
}
