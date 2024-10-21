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
##' @export
##' @examples
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
