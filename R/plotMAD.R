#' @title Plots of Median Absolute Deviation
#'
#' @description For each pipeline, two quantification replicates
#' are compared and log scale absolute deviations of signals are
#' calculated. Then, loess smooths on absolute deviation are
#' plotted stratefied by detrended log signals.
#'
#' @param dat A \code{rnaseqcomp} S4 class object.
#' @param type Plot types (default: 'l').
#' @param lwd Plot line weights (default: 2).
#' @param col Plot colors (default: NULL, colors are assigned
#' by package \code{RColorBrewer}).
#' @param xlim Plot limits of x-axis (default: NULL, limits are
#' estimated automatically).
#' @param ylim Plot limits of y-axis (default: NULL, limits are
#' estimated automatically).
#' @param xlab Plot label of x-axis
#' (default: 'Detrended logSignal').
#' @param ylab Plot label of y-axis (default: 'MAD').
#' @param ... Other parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{MAD plot}{MAD plots for all the quantification pipelines.}
#' \item{MAD}{A numeric vector of median absolute deviations.}
#'
#' @export
#' @examples
#' data(encodeCells)
#' txFIdx <- encodeCells$genemeta$type == "protein_coding"
#' hkIdx <- encodeCells$genemeta$housekeeping
#' unitFIdx <- grepl("Cufflinks",encodeCells$repInfo)
#' dat <- matrixFilter(encodeCells$gm12878,encodeCells$repInfo,
#' txFIdx,hkIdx,unitFIdx)
#' plotMAD(dat)

plotMAD <- function(dat, type='l', lwd = 2, col = NULL,
                   xlim = NULL, ylim = NULL,
                    xlab = "Detrended logSignal", ylab = "MAD",
                   ...){
    if(!is(dat, 'rnaseqcomp'))
        stop('"plotMAD" only plots class "rnaseqcomp".')
    cdList <- lapply(levels(dat@repInfo), function(i)
                     dat@quantData[, dat@repInfo == i])
    hkmed <- lapply(levels(dat@repInfo), function(i)
                    dat@hkmed[dat@repInfo == i])
    # positive units to detrended log signal
    cdList <- lapply(seq_len(length(cdList)), function(i){
        tmp <- log2(cdList[[i]][which(apply(cdList[[i]], 1, min) > 0), ])
        t(t(tmp) - hkmed[[i]]) + dat@scaler
    })
    names(cdList) <- names(hkmed) <- levels(dat@repInfo)
    # M-|A|
    sdlist <- lapply(cdList, function(x)
                     cbind(rowMeans(x), abs(x[,1] - x[,2])))
    if(is.null(xlim))
        xlim <- c(min(sapply(sdlist, function(x) min(x[ ,1]))),
                  max(sapply(sdlist, function(x) max(x[ ,1]))))
    if(is.null(ylim))
        ylim <- c(min(sapply(sdlist, function(x) quantile(x[ ,2], 0.3))),
                  max(sapply(sdlist, function(x) quantile(x[ ,2], 0.75))))
    if(is.null(xlab))  xlab <- 'Detrended logSignal'
    if(is.null(ylab))  ylab <- 'MAD'
    if(is.null(col))   col <- brewer.pal(min(length(sdlist), 9), "Set1")
    col <- rep_len(col, length(sdlist))
    type <- rep_len(type, length(sdlist))
    lwd <- rep_len(lwd, length(sdlist))
    for(i in seq_len(length(sdlist))){
        # loess smooth
        x <- loess.smooth(sdlist[[i]][ ,1], sdlist[[i]][ ,2], span = 2/3,
                          degree = 1, family = "symmetric", evaluation = 1000)
        if(i == 1) {
            plot(x$x, x$y, type = type[i], lwd = lwd[i], col = col[i],
                 ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, ...)
        }else {
            points(x$x, x$y, type = type[i], lwd = lwd[i], col = col[i])
        }
    }
    legend('topright', names(cdList), lwd = lwd, col = col, cex=0.5)
    # MAD
    cat("One number statistics: MAD\n")
    return(sapply(sdlist, function(x) round(median(x[ ,2]), 3)))
}

