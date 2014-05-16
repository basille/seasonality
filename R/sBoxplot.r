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
