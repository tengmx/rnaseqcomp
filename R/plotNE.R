#' @title Plots of Express And Non-express Features
#'
#' @description For each pipeline, two quantification replicates
#' are compared and proportions of both-express, both-non-express
#' and either-or-express features are calculated. Then, reverse
#' proportion accumulation for either-or-express features are
#' plotted stratified by detrended log signals.
#'
#' @param dat A \code{rnaseqcomp} S4 class object.
#' @param step Plot steps on x-axis.
#' @param type Plot types (default: 'l').
#' @param lwd Plot line weights (default: 2).
#' @param col Plot colors (default: NULL, colors are assigned
#' by package \code{RColorBrewer}).
#' @param lty Plot line styles (default: 1).
#' @param xlim Plot limits of x-axis (default: NULL, limits are
#' estimated automatically).
#' @param ylim Plot limits of y-axis (default: NULL, limits are
#' estimated automatically).
#' @param cex.leg Legend size (default: 0.6).
#' @param xlab Plot label of x-axis
#' (default: 'Detrended logSignal').
#' @param ylab Plot label of y-axis
#' (default: 'Reverse Accumulation Proportion of NE'').
#' @param ... Other parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{Either-or-express plot}{Plots for all the quantification
#' pipelines.}
#' \item{matrix}{A proportion matrix of express, non-express and
#' either-or-express for all pipelines.}
#'
#' @export
#' @examples
#' data(encodeCells)
#' evaluationFeature <- encodeCells$genemeta$type == "protein_coding"
#' calibrationFeature <- encodeCells$genemeta$housekeeping
#' unitReference <- grepl("Cufflinks",encodeCells$repInfo)
#' dat <- matrixFilter(encodeCells$gm12878,encodeCells$repInfo,
#'     evaluationFeature,calibrationFeature,unitReference)
#' plotNE(dat)

plotNE <- function(dat, step = 0.1, type = 'l', lwd = 2, col = NULL,
                   lty = 1, xlim = NULL, ylim = NULL, cex.leg = 0.6,
                   xlab = "Detrended logSignal",
                   ylab = "Reverse Accumulation Proportion of NE",
                   ...){
    if(!is(dat,'rnaseqcomp'))
        stop('"plotNE" only plots class "rnaseqcomp".')
    cdList <- lapply(levels(dat@repInfo), function(i)
                     dat@quantData[ ,dat@repInfo == i])
    refMed <- lapply(levels(dat@repInfo), function(i)
                    dat@refMed[dat@repInfo == i])
    # proportion of EE, NN & NE, missing data excluded
    pEE <- round(sapply(cdList, function(x)
                        mean(apply(x, 1, min) > 0, na.rm = TRUE)), 3)
    pNN <- round(sapply(cdList, function(x)
                        mean(apply(x, 1, max) <= 0, na.rm = TRUE)), 3)
    pNE <- 1 - pEE - pNN
    # NE data to detrended log signal
    neList <- lapply(seq_len(length(cdList)), function(i){
        tmp1 <- log2(cdList[[i]][which(apply(cdList[[i]], 1, min) <= 0 &
                                       apply(cdList[[i]], 1, max) > 0), ])
        tmp2 <- t(t(tmp1) - refMed[[i]]) + dat@scaler
        apply(tmp2 * !is.infinite(tmp2), 1, sum, na.rm = TRUE)
    })
    if(is.null(xlim))
        xlim <- c(min(sapply(neList, function(x) min(x))),
                  max(sapply(neList, function(x) max(x))))
    # detrended signal steps & reverse accumulated proportions
    k <- seq(xlim[1], xlim[2], step)
    pnelist <- lapply(seq_len(length(neList)), function(i)
                      sapply(k, function(x)
                             sum(neList[[i]] > x) / nrow(cdList[[i]])))
    names(pnelist) <- names(refMed) <- levels(dat@repInfo)
    if(is.null(ylim))
        ylim <- c(0, max(sapply(pnelist, function(x) max(x))))
    if(is.null(xlab))  xlab <- 'Detrended logSignal'
    if(is.null(ylab))  ylab <- 'Reverse Accumulation Proportion of NE'
    if(is.null(col))   col <- brewer.pal(min(length(pnelist), 9), "Set1")
    col <- rep_len(col, length(pnelist))
    type <- rep_len(type, length(pnelist))
    lwd <- rep_len(lwd, length(pnelist))
    lty <- rep_len(lty, length(pnelist))
    for(i in seq_len(length(pnelist))){
        if(i == 1) {
            plot(k, pnelist[[i]], type = type[i], lwd = lwd[i], col = col[i],
                 lty = lty[i],
                 ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, ...)
        }else {
            points(k, pnelist[[i]], type = type[i], lwd = lwd[i],col = col[i],
                   lty = lty[i])
        }
    }
    legend('topright', names(pnelist), lwd = lwd, col = col,
           lty = lty, cex = cex.leg)
    dat <- cbind(pEE, pNE, pNN)
    rownames(dat) <- names(pnelist)
    return(dat)
}

