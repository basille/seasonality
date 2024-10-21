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
##' @export
##' @examples
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
