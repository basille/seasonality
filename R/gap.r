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
