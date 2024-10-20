* bsSeasons: remove 'datanorm' (Ivar Herfindal 2014/02 + Endre Grüner Ofstad 2016/02)

library(seasonality)

### Load the data
data(caribou)

### Recompute the bootstrap seasons:
debug(bsSeasons)
caribou$bs <- bsSeasons(data = caribou$move, ind = caribou$ind, dataNA = caribou$window, nclust = 8)
caribou$bs <- bsSeasons2(data = caribou$move, ind = caribou$ind, dataNA = caribou$window, nclust = 8)

bla <- caribou


bsSeasons <- function(data, ind, dataNA, nclust, iter = 100,
                      simplify = FALSE, win = 6, tol = 2) {
    bs <- list()
    ## 'datanorm <-' for no reason
    ## datanorm <- bsi <- function(i) {
    bsi <- function(i) {
        set.seed(i)
        rd <- sort(sample(1:nrow(ind), nrow(ind), replace = TRUE))
        names <- ind[rd, ]
        names$weights <- rep(1 / table(names$id), times = table(names$id))
        norm <- data.frame(row.names = 1:365)
        for (j in names(data)) {
            norm[, j] <- apply(data[[j]][
                ,
                rd
            ], 1, weighted.mean, w = names$weights, na.rm = TRUE)
        }
        norm <- as.data.frame(scale(norm, center = apply(norm,
            2, min,
            na.rm = TRUE
        ), scale = apply(norm, 2, function(x) diff(range(x,
                na.rm = TRUE
            )))))
        if (any(is.na(norm))) {
            ## data$norm doesn't exist...
            ## norm[is.na(norm)] <- data$norm[is.na(norm)]
            norm[is.na(norm)] <- dataNA[is.na(norm)]
        }
        set.seed(1)
        seasons <- kmeans(norm, nclust, iter.max = 100)$cluster
        if (simplify) {
            return(sFormat(sSimple(seasons, win, tol)))
        } else {
            return(sFormat(seasons))
        }
    }
    bs <- lapply(1:iter, bsi)
    names(bs) <- 1:iter
    return(bs)
}

### Compute the weights, and identify the final seasons:
set.seed(1)
seasons <- kmeans(caribou$window, 8, iter.max = 100)$cluster
weights <- bsWeights(caribou$bs)
seasonsbs <- bsCriterion(seasons, weights, threshold = .8)

seasons2 <- kmeans(bla$window, 8, iter.max = 100)$cluster
weights2 <- bsWeights(bla$bs)
seasonsbs2 <- bsCriterion(seasons2, weights2, threshold = .8)


### Visualize the final seasons:
bsPlot(seasonsbs, seasons, weights, title = "Caribou")
bsPlot(seasonsbs2, seasons2, weights2, title = "Caribou2")

