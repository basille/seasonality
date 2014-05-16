##' seasonality package
##'
##' This package provides a set of functions to split year-round
##' space-use measurements into biological seasons, completed with
##' additional functions to explore and simplify these
##' seasons. Reference: Basille M., Fortin D., Dussault C., Ouellet
##' J.-P., Courtois R. Ecologically based definition of seasons
##' clarifies predator-prey interactions. Ecography, early view. DOI:
##' 10.1111/j.1600-0587.2011.07367.x
##'
##' @name seasonality
##' @docType package
##' @author Mathieu Basille \email{basille@@ase-research.org}
NULL


## caribou
##
##' Data set used in Basille et al. (2012).
##'
##' @name caribou
##' @docType data
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @references Basille M., Fortin D., Dussault C., Ouellet J.-P.,
##' Courtois R. Ecologically based definition of seasons clarifies
##' predator-prey interactions. Ecography, early view. DOI:
##' 10.1111/j.1600-0587.2011.07367.x
##' @keywords data
NULL


## gap
##
##' Compute the gap statistic (weighted by default).
##'
##' The package \code{clusterSim} proposes a
##' \code{\link[clusterSim]{index.Gap}} function to compute the gap
##' statistic. It can be used with many different clustering methods
##' ("ward", "single", "complete", "average", "mcquitty", "median",
##' "centroid", "pam", "k-means", "diana"), and with uniform or
##' pc-based reference distributions.
##'
##' Bram Van Moorter modified it into \code{index.gap.modif}
##' (\url{http://ase-research.org/moorter/p7_gap.statistic.r}), which
##' uses k-means as a default, returns values when only one large
##' cluster is made, and instead of calculating gap-differences, it
##' now returns the original gap-value.
##'
##' It seems however that the algorithm to compute \eqn{W_k} in
##' \code{\link[clusterSim]{index.Gap}} is not correct; in addition
##' the \code{\link[clusterSim]{index.Gap}} function is quite poorly
##' written and thus difficult to understand; last but not least, it
##' does not allow to compute the weighted gap statistic. The weighted
##' gap statistic have been shown to provide more robust and
##' consistent results, and allows in a multi-layer approach to derive
##' nested clusters.
##' @title Gap statistic
##' @param data A matrix, or a data frame coercible to a matrix. Input
##' data should be of the form \code{obs} \eqn{\times} \code{var}.
##' @param from The minimal number of clusters for which the gap
##' statistic is computed.
##' @param to The maximal number of clusters for which the gap
##' statistic is computed.
##' @param nsim The number of simulations used to compute the gap
##' statistic.
##' @param ref.dist A character string specifying the reference
##' distribution:
##' \describe{
##' \item{\code{unif}}{Generates each reference variable uniformly
##' over the range of the observed values for that variable;}
##' \item{\code{pc}}{Generates the reference variables from a uniform
##' distribution over a box aligned with the principal components of
##' the data.}}
##' @param clust.method A character string specifying the cluster
##' analysis method to be used. This should be one of: \code{"ward"},
##' \code{"single"}, \code{"complete"}, \code{"average"},
##' \code{"mcquitty"}, \code{"median"}, \code{"centroid"},
##' \code{"pam"}, \code{"k-means"}, \code{"diana"}. Only tested for
##' \code{"k-means"}, which is the default.
##' @param dist.method The distance measure to be used. Only tested
##' for \code{"euclidean"}. See \code{\link[stats]{dist}} for other
##' metrics.
##' @param weighted Logical. Whether the gap statistic should be
##' weighted or not (default is \code{TRUE}).
##' @param tol An number specifying the multiplier to reject the null
##' model. The tolerance is analogous to setting the alpha level in
##' the standard hypothesis testing framework, where increased
##' tolerance is similar to selecting a smaller alpha rejection
##' region.  Tibshirani et al. (2001) used a tolerance of 1 (default
##' behaviour), but larger values of tolerance increase the strength
##' of evidence required to include additional clusters;
##' @param seed A single value, interpreted as an integer, used a seed
##' in the clustering method.
##' @return \code{gap} returns a \code{k x p} data frame of class
##' \code{gap} with the following variables:
##' \describe{
##' \item{nCluster}{The number of clusters k;}
##' \item{logWk0}{\eqn{\log(W_k)} (from the data) where \eqn{W_k =
##' \sum^k_{m=1} \frac{1}{2n_m} D_m}, or \eqn{W_k = \sum^k_{m=1}
##' \frac{1}{2n_m(n_m -1)} D_m} if weighted, \eqn{D_m} being the
##' (complete) sum of pairwise distances;}
##' \item{logWk}{\eqn{E^*_n{\log(W_k)}} (from the simulated data
##' sets);}
##' \item{Gap}{The gap statistic as \eqn{\mathrm{Gap}_n(k) =
##' E^*_n\{\log(W_k)\} - \log(W_k) = (1/B)\sum^B_b=1 \log(W^*_{kb}) -
##' \log(W_k)}, \code{B} being the number of simulated data sets;}
##' \item{sdGap}{The standard deviation of the gap statistic, as
##' \eqn{[s_k = (1/B) \{\sum^B_b=1 \log(W^*_{kb}) - (1/B) \sum^B_b=1
##' \log(W^*_{kb})\}^2]^{1/2} \sqrt{(1+1/B)}};}
##' \item{k}{The estimated number of clusters with the classical
##' approach, indicated by an asterisk, with a tolerance \code{T},
##' such as \eqn{\mathrm{Gap}(k) \geq \mathrm{Gap}(k+1) - T *
##' s_{k+1}};}
##' \item{D}{Differences of gap, as \eqn{D\mathrm{Gap}_n(k) =
##' \mathrm{Gap}_n(k) - \mathrm{Gap}_n(k-1)};}
##' \item{DD}{differences of Dgap, as \eqn{DD\mathrm{Gap}_n(k) =
##' D\mathrm{Gap}_n(k) - D\mathrm{Gap}_n(k + 1)};}
##' \item{DDk}{The estimated number of clusters with the DD-weighted
##' approach, indicated by an asterisk; the number of clusters
##' \code{k} is given when DDGap is maximum.}}
##' @references Tibshirani, R.; Walther, G. & Hastie, T. (2001)
##' Estimating the number of clusters in a data set via the gap
##' statistic. Journal of the Royal Statistical Society: Series B
##' (Statistical Methodology), Blackwell Publishers Ltd., 63: 411-423,
##' DOI: 10.1111/1467-9868.00293
##'
##' Yan, M. & Ye, K. (2007) Determining the number of clusters using
##' the weighted gap statistic. Biometrics, 63: 1031-1037, DOI:
##' 10.1111/j.1541-0420.2007.00784.x
##'
##' Basille, M.; Fortin, D.; Dussault, C.; Ouellet, J.-P. & Courtois,
##' R. (2012) Ecologically based definition of seasons clarifies
##' predator-prey interactions. Ecography, early view, DOI:
##' 10.1111/j.1600-0587.2011.07367.x
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ### Simple simulation
##' set.seed(1)
##' X <- matrix(rnorm(30, mean = 5), ncol = 3)
##' set.seed(1)
##' Y <- rbind(matrix(rnorm(300, mean = 5), ncol = 3),
##'     matrix(rnorm(300, mean = 10), ncol = 3))
##'
##' ### K-means tests
##' ## Beware of the case of only 1 group:
##' (GG1 <- gap(X, to = 9, ref.dist = "unif"))
##' plot(GG1)
##' ## Two groups:
##' (GG2 <- gap(Y))
##' plot(GG2)
##'
##' ### Caribou data
##' data(caribou)
##' carigap <- gap(caribou$window)
##' plot(carigap)
gap <- function(data, from = 1, to = 10, nsim = 50, ref.dist = c("pc",
    "unif"), clust.method = "k-means", dist.method = "euclidean",
    weighted = TRUE, tol = 1, seed = 1)
{
    clust <- function(mat, nclust, clust.method) {
        set.seed(seed)
        if (nclust == 1)
            return(rep(1, nrow(mat)))
        else if (clust.method == "pam") {
            require(cluster)
            return(pam(mat, nclust)$cluster)
        }
        else if (clust.method == "k-means")
            return(kmeans(mat, nclust, 100)$cluster)
        else if (clust.method == "diana")
            return(cutree(as.hclust(diana(dist(mat))), k = nclust))
        else if (clust.method %in% c("single", "complete", "average",
            "ward", "mcquitty", "median", "centroid"))
            return(cutree(hclust(dist(mat), method = clust.method),
                nclust))
        else stop("Wrong clustering method")
    }
    computeWk <- function(mat, cl, dist.method = "euclidean",
        weighted = TRUE) {
        Wk <- 0
        for (cli in unique(cl)) {
            if (!weighted) {
                if (dist.method == "euclidean")
                  Wk <- Wk + sum(diag(var(mat[cl == cli, ]))) *
                    (length(cl[cl == cli]) - 1)
                else {
                  D <- sum(dist(mat[cl == cli, ], method = dist.method)^2)
                  Wk <- Wk + D/(nrow(mat[cl == cli, , drop = FALSE]))
                }
            }
            else {
                if (dist.method == "euclidean")
                  Wk <- Wk + sum(diag(var(mat[cl == cli, ])))
                else {
                  D <- sum(dist(mat[cl == cli, ], method = dist.method)^2)
                  Wk <- Wk + D/(nrow(mat[cl == cli, , drop = FALSE]) *
                    (nrow(mat[cl == cli, , drop = FALSE]) - 1))
                }
            }
        }
        return(Wk)
    }
    simX <- function(mat, seed = seed, ref.dist = ref.dist) {
        set.seed(seed)
        simunif <- function(x) return(apply(x, 2, function(w) runif(length(w),
            min = min(w), max = max(w))))
        simpc <- function(x) {
            x <- as.matrix(x)
            xmm <- apply(x, 2, mean)
            for (k in (1:dim(x)[2])) x[, k] <- x[, k] - xmm[k]
            ss <- svd(x)
            xsim <- simunif(x %*% ss$v) %*% t(ss$v)
            for (k in (1:dim(x)[2])) xsim[, k] <- xsim[, k] +
                xmm[k]
            return(xsim)
        }
        if (ref.dist == "unif")
            return(simunif(mat))
        if (ref.dist == "pc")
            return(simpc(mat))
    }
    ref.dist <- match.arg(ref.dist)
    data <- as.matrix(data)
    Gap <- data.frame(matrix(NA, nrow = to - from + 1, ncol = 9))
    names(Gap) <- c("nCluster", "logWk0", "logWk", "Gap", "sdGap",
        "k", "D", "DD", "DDk")
    for (i in from:to) {
        Gap[i, "nCluster"] <- i
        Gap[i, "logWk0"] <- log(computeWk(data, clust(data, i,
            clust.method = clust.method), dist.method = dist.method,
            weighted = weighted))
        Sim <- lapply(1:nsim, function(seed) simX(data, seed = seed,
            ref.dist = ref.dist))
        logWksim <- log(unlist(lapply(Sim, function(xl) computeWk(xl,
            clust(xl, i, clust.method = clust.method), dist.method = dist.method,
            weighted = weighted))))
        Gap[i, "logWk"] <- mean(logWksim)
        Gap[i, "Gap"] <- Gap[i, "logWk"] - Gap[i, "logWk0"]
        Gap[i, "sdGap"] <- sqrt(1 + 1/nsim) * sqrt(var(logWksim) *
            (nsim - 1)/nsim)
    }
    Gap$k <- ""
    Gap$k[min(which(Gap$Gap[-nrow(Gap)] >= c(Gap$Gap - (tol *
        Gap$sdGap))[-1]))] <- "*"
    Gap$D[2:nrow(Gap)] <- diff(Gap$Gap)
    Gap$DD[2:(nrow(Gap) - 1)] <- Gap$D[2:(nrow(Gap) - 1)] - Gap$D[3:nrow(Gap)]
    Gap$DDk <- ""
    Gap$DDk[which(Gap$DD == max(Gap$DD, na.rm = TRUE))] <- "*"
    class(Gap) <- c("gap", "data.frame")
    return(Gap)
}


## plot.gap
##
##' @param x An object of class \code{gap}.
##' @param ... Further arguments passed to or from other methods.
##' @rdname gap
##' @export
plot.gap <- function(x, ...) {
    if (!inherits(x, "gap"))
        stop("Object of class 'gap' expected")
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mfrow = c(2, 2))
    k <- seq(1:nrow(x))
    plot(k, exp(x$logWk0), xlab = "Number of clusters k", ylab = "Wk",
        type = "b", ...)
    plot(k, x$Gap, xlab = "Number of clusters k", ylab = "GAP",
        type = "b", ...)
    segments(k, c(x$Gap - x$sdGap), k, c(x$Gap + x$sdGap))
    points(x$nCluster[x$k == "*"], x$Gap[x$k == "*"], col = "black",
        pch = 19, cex = 1.5)
    plot(k, x$D, xlab = "Number of clusters k", ylab = "D", type = "b",
        ...)
    plot(k, x$DD, xlab = "Number of clusters k", ylab = "DD",
        type = "b", ...)
    points(x$nCluster[x$DDk == "*"], x$DD[x$DDk == "*"], col = "black",
        pch = 19, cex = 1.5)
}


## sPrint
##
##' Print a sequence of seasons in a friendly way.
##'
##' @title Print seasons
##' @param seasons The result of a season clustering: A vector of
##' integers (from \code{1:k}) indicating the cluster to which each
##' day is allocated.
##' @param ndays Logical. Returns the rank of the days at which a new
##' season starts.
##' @return A vector indicating the dates at which a new season
##' starts.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(caribou)
##' set.seed(1)
##' seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
##' sPrint(seasons)
sPrint <- function(seasons, ndays = FALSE) {
    br <- c(1, which(diff(seasons) != 0) + 1, length(seasons))
    brs <- seasons[br]
    names(brs) <- format(strptime(paste(2011, br + as.numeric(names(seasons[1])) -
        1), format = "%Y %j"), "%d/%m")
    if (ndays)
        return(br + as.numeric(names(seasons[1])) - 1)
    else {
        if (brs[length(brs)] == brs[length(brs) - 1])
            brs[length(brs)] <- NA
        return(brs)
    }
}


## sFormat
##
##' Reorder the seasons as a succession of unique numbers, from 1 to
##' the last season (useful in case of duplicated clusters, as
##' duplicates get a new index).
##'
##' @title Reorder the seasons.
##' @param seasons The result of a season clustering: A vector of
##' integers (from \code{1:k}) indicating the cluster to which each
##' day is allocated.
##' @return A vector of the same length as \code{seasons}, with the
##' index of the clusters reordered.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(caribou)
##' set.seed(1)
##' seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
##' sPrint(seasons)
##' sPrint(sFormat(seasons))
sFormat <- function(seasons)
{
    first <- names(seasons)[1]
    last <- seasons[1] == seasons[length(seasons)]
    seasons <- cumsum(c(1, diff(seasons) != 0))
    names(seasons)[1] <- first
    if (last)
        seasons[seasons == max(seasons)] <- 1
    return(seasons)
}


## sSimple
##
##' Simplify the seasons after the initial clustering, by removing the
##' smallest seasons.
##'
##' The function works on a moving window of length (current day +
##' \code{win}) days. For a given day, if all other days (with a
##' tolerance of \code{tol} days) have the same value as the focus
##' day, the day is kept as is; otherwise, the day takes the value of
##' the day before.
##'
##' @title Simplify the seasons.
##' @param clust The result of a season clustering: A vector of
##' integers (from \code{1:k}) indicating the cluster to which each
##' day is allocated.
##' @param win The length of the moving window.
##' @param tol The tolerance to be used.
##' @return A vector of the same length as \code{seasons}.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(caribou)
##' set.seed(1)
##' seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
##' sPrint(seasons)
##' sPrint(sSimple(seasons))
sSimple <- function(clust, win = 3, tol = 1)
{
    clust <- c(clust[length(clust)], clust, clust[1:win])
    for (i in 2:(length(clust) - win)) if (sum(clust[(i + 1):(i +
        win)] == clust[i]) < (win - tol))
        clust[i] <- clust[i - 1]
    ## 2nd identical loop to ensure that there are only 3-days
    ## seasons, in case of 22112122 which gives 22112222
    for (i in 2:(length(clust) - win)) if (sum(clust[(i + 1):(i +
        win)] == clust[i]) < (win - tol))
        clust[i] <- clust[i - 1]
    return(clust[-c(1, (length(clust) - win + 1):length(clust))])
}


## sPlot
##
##' Plot the seasons.
##'
##' @title Plot the seasons
##' @param seasons The result of a season clustering: A vector of
##' integers (from \code{1:k}) indicating the cluster to which each
##' day is allocated.
##' @param add.lines Logical. Adds dotted lines delineating the
##' seasons.
##' @param months Logical. Draws monthtly delineations.
##' @param main An overall title for the plot.
##' @param ylab A title for the y axis.
##' @param ... Further arguments passed to the \code{lines} call.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(caribou)
##' set.seed(1)
##' seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
##' sPlot(sSimple(seasons))
sPlot <- function(seasons, add.lines = FALSE, months = FALSE,
    main = "Seasons", ylab = substitute(seasons), ...)
{
    if (add.lines)
        lines(seasons, ...)
    else {
        par(mfrow = c(1, 1), mar = c(4, 5, 3, 1) + 0.1)
        plot(seasons, type = "n", axes = FALSE, main = main,
            xlab = "", ylab = ylab)
        Range <- range(seasons)
        par(usr = c(1, 365, Range[1] - 0.03 * diff(Range), Range[2] +
            0.03 * diff(Range)))
        if (months)
            rect(c(32, 91, 152, 213, 274, 335), Range[1] - 0.03 *
                diff(Range), c(60, 121, 182, 244, 305, 365),
                Range[2] + 0.03 * diff(Range), col = grey(0.9),
                border = NA)
        else abline(v = as.numeric(names(which(diff(seasons) !=
            0))), lty = 3, lwd = 0.5)
        lines(seasons, ...)
        if (months)
            axis(1, at = c(16, 46, 75.5, 106, 136.5, 167, 197.5,
                228.5, 259, 289.5, 320, 350), labels = c("J",
                "F", "M", "A", "M", "J", "J", "A", "S", "O",
                "N", "D"), tick = FALSE, line = -1, cex.axis = 0.8)
        else axis(1, at = as.numeric(names(which(diff(seasons) !=
            0))), labels = format(strptime(names(which(diff(seasons) !=
            0)), format = "%j"), "%d/%m"), cex.axis = 0.8, las = 2)
        axis(2)
        box()
    }
}


## sBoxplot
##
##' Season boxplots.
##'
##' @title Season boxplots
##' @param data The original data on which the clustering was made
##' (see \code{\link{gap}}).
##' @param seasons The result of a season clustering: A vector of
##' integers (from \code{1:k}) indicating the cluster to which each
##' day is allocated.
##' @param temporal Logical. If \code{TRUE}, produces boxplots along
##' the year (X-axis); if \code{FALSE}, produces boxplots for each
##' cluster using their index.
##' @param months Draws the months with background rectangles
##' (\code{rectangle}) or dotted lines (\code{lines}).
##' @param cluster Logical. Indicates the cluster index above the
##' graph.
##' @param multi Logical. Allows for comparison between several
##' clusterings, by displaying them side by side. If \code{yes},
##' requires a list of \code{data} and \code{seasons}, corresponding
##' to each clustering.
##' @param samescale Logical. In case of comparison, use the same
##' scale for common variables.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(caribou)
##' set.seed(1)
##' seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
##' sBoxplot(caribou$window, sSimple(seasons))
sBoxplot <- function(data, seasons, temporal = TRUE, months = c("rectangles",
    "lines"), cluster = TRUE, multi = FALSE, samescale = TRUE) {
    old.par <- par(no.readonly = TRUE)
    if (multi)
        par(mfcol = c(ncol(data[[1]]), length(data)), mar = c(1,
            2, 3, 0) + 0.1)
    else {
        par(mfrow = n2mfrow(ncol(data)), mar = c(1, 2, 3, 0) +
            0.1)
        data <- list(data)
        seasons <- list(seasons)
    }
    on.exit(par(old.par))
    months <- match.arg(months)
    if (temporal) {
        for (j in 1:length(data)) {
            datarb <- do.call(rbind, data)
            datatmp <- data[[j]]
            seasonstmp <- seasons[[j]]
            changes <- as.numeric(c(names(seasonstmp)[1], names(which(diff(seasonstmp) !=
                0)), names(seasonstmp)[length(seasonstmp)]))
            seas <- seasonstmp[as.character(changes[-length(changes)])]
            at <- changes[-length(changes)] + diff(changes)/2
            for (i in 1:ncol(datatmp)) {
                summ <- do.call(rbind, lapply(1:max(seasonstmp),
                  function(j) summary(datatmp[seasonstmp == j,
                    i])))
                if (samescale) {
                  plot(as.numeric(row.names(datatmp)), datatmp[,
                    i], type = "n", xlim = c(1, 365), ylim = range(datarb[,
                    i]), axes = FALSE, main = names(datatmp)[i])
                  par(usr = c(1, 365, min(datarb[, i]), max(datarb[,
                    i])))
                }
                else {
                  plot(as.numeric(row.names(datatmp)), datatmp[,
                    i], type = "n", xlim = c(1, 365), ylim = range(datatmp[,
                    i]), axes = FALSE, main = names(datatmp)[i])
                  par(usr = c(1, 365, min(datatmp[, i]), max(datatmp[,
                    i])))
                }
                if (months == "lines")
                  abline(v = c(32, 60, 91, 121, 152, 182, 213,
                    244, 274, 305, 335), lty = 3, col = grey(0.6))
                if (months == "rectangles")
                  rect(c(32, 91, 152, 213, 274, 335), min(datarb[,
                    i]), c(60, 121, 182, 244, 305, 365), max(datarb[,
                    i]), col = grey(0.9), border = NA)
                segments(at, summ[seas, "Min."], y1 = summ[seas,
                  "Max."], lty = 2)
                for (j in 1:length(seas)) polygon(x = changes[c(j,
                  j + 1, j + 1, j)] + c(1, -1, -1, 1), y = rep(summ[seas[j],
                  c("1st Qu.", "3rd Qu.")], each = 2), col = "white")
                segments(changes[-length(changes)] + 1, summ[seas,
                  "Median"], changes[-1] - 1, lwd = 3)
                axis(2)
                axis(1, at = c(16, 46, 75.5, 106, 136.5, 167,
                  197.5, 228.5, 259, 289.5, 320, 350), labels = c("J",
                  "F", "M", "A", "M", "J", "J", "A", "S", "O",
                  "N", "D"), tick = FALSE, line = -1, cex.axis = 0.8)
                if (cluster)
                  axis(3, at = at, labels = seas, tick = FALSE,
                    line = -1, cex.axis = 0.8)
                box()
            }
        }
    }
    else {
        for (j in 1:length(data)) {
            datatmp <- data[[j]]
            seasonstmp <- seasons[[j]]
            for (i in 1:ncol(datatmp)) {
                boxplot(datatmp[, i] ~ seasonstmp, axes = FALSE,
                  main = names(datatmp)[i], ylim = c(0, 1))
                axis(1)
                box()
            }
        }
    }
}


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
