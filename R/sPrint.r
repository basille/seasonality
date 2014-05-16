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
