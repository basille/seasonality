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
##' @export
##' @examples
##' data("caribou")
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
