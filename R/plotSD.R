#' @title Estimate And Plot Median Standard Deviation
#'
#' @description For each cell line in each pipeline,
#' the standard deviation of detrend logsignals are calculated
#' for individual features. Then, loess smooth on standard
#' deviation are plotted stratified by detrended log signals
#' for select cell line. The median of standard deviation at
#' three different levels of detrend logsignals are reported.
#'
#' @param dat A \code{rnaseqcomp} S4 class object.
#' @param constant A numeric pseudo-constant to be added on all
#' the signals before transferred to log scale. (default: 0.5)
#' @param loessspan A numeric number indicating span used
#' for loess smooth. Details see
#' \code{loess.smooth} function. (Default: 1/3)
#' @param thresholds A vector of two numbers define cutoffs for
#' three levels of detreded log signals. (default: c(1, 6))
#' @param plotcell Either 1 or 2 indicating which cell line
#' will be plotted. This won't affect estimation for both
#' cell lines. (default: 1)
#' @param lwd Plot line weights (default: 2).
#' @param col Plot colors (default: NULL, colors are assigned
#' by package \code{RColorBrewer}).
#' @param lty Plot line styles (default: 1).
#' @param cex.leg Legend size (default: 1).
#' @param xlim Plot limits of x-axis (default: NULL, limits are
#' estimated automatically).
#' @param ylim Plot limits of y-axis (default: NULL, limits are
#' estimated automatically).
#' @param xlab Plot label of x-axis
#' (default: 'Detrended logSignal').
#' @param ylab Plot label of y-axis (default: 'SD').
#' @param ... Other parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#' @import methods
#'
#' @return
#' \item{plot}{SD plots of quantification pipelines for
#' selected cell line by \code{plotcell}.}
#' \item{SD}{A list of two matrices of median standard deviations.}
#'
#' @export
#' @examples
#' data(simdata)
#' condInfo <- factor(simdata$samp$condition)
#' repInfo <- factor(simdata$samp$replicate)
#' evaluationFeature <- rep(TRUE, nrow(simdata$meta))
#' calibrationFeature <- simdata$meta$house & simdata$meta$chr == 'chr1'
#' unitReference <- 1
#' dat <- signalCalibrate(simdata$quant, condInfo, repInfo, evaluationFeature,
#' calibrationFeature, unitReference, calibrationFeature2 = calibrationFeature)
#' plotSD(dat)


plotSD <- function(dat, constant = 0.5, loessspan = 1/3,
                   thresholds = c(1, 6), plotcell = 1,
                   lwd = 2, col = NULL, lty=1, cex.leg=1,
                   xlim = NULL, ylim = NULL,
                   xlab = "Detrended logSignal", ylab = "SD",
                   ...){
    if(!is(dat, 'rnaseqcomp'))
        stop('"plotSD" only plots class "rnaseqcomp".')
    dat@quantData <- lapply(dat@quantData, function(x) x + constant)
    cdList <- list()
    for(i in 1:2){
        cdList[[i]] <- lapply(dat@quantData, function(x)
                        x[,dat@condInfo == levels(dat@condInfo)[i], drop=F])
    }
    sdlist1 <- lapply(seq_along(dat@quantData),function(j){
        tmp <- cdList[[1]][[j]]
        idx <- which(apply(tmp, 1, max) > constant)
        tmp <- log2(tmp[idx,])
        cbind(rowMeans(tmp), apply(tmp, 1, sd))
    })
    sdlist2 <- lapply(seq_along(dat@quantData),function(j){
        tmp <- cdList[[2]][[j]]
        idx <- which(apply(tmp, 1, max) > constant)
        tmp <- log2(tmp[idx,])
        cbind(rowMeans(tmp), apply(tmp, 1, sd))
    })
    names(sdlist1) <- names(sdlist2) <- names(dat@quantData)
    if(is.null(xlim))
        xlim <- c(min(c(sapply(sdlist1, function(x) min(x[ ,1])),
                        sapply(sdlist2, function(x) min(x[ ,1])))),
                  max(c(sapply(sdlist1, function(x) max(x[ ,1])),
                        sapply(sdlist2, function(x) max(x[ ,1])))))
    if(is.null(ylim))
        ylim <- c(min(c(sapply(sdlist1, function(x) quantile(x[ ,2], 0.3)),
                        sapply(sdlist2, function(x) quantile(x[ ,2], 0.3)))),
                  max(c(sapply(sdlist1, function(x) quantile(x[ ,2], 0.75)),
                        sapply(sdlist2, function(x) quantile(x[ ,2], 0.75)))))
    if(is.null(col))   col <- brewer.pal(min(length(sdlist1), 8), "Set2")
    col <- rep_len(col, length(sdlist1))
    lty <- rep_len(lty, length(sdlist1))
    for(i in seq_len(length(sdlist1))){
        if(plotcell==1){
            fitx <- sdlist1[[i]][ ,1]
            fity <- sdlist1[[i]][ ,2]
        }else if(plotcell==2){
            fitx <- sdlist2[[i]][ ,1]
            fity <- sdlist2[[i]][ ,2]
        }else{
            fitx <- sdlist1[[i]][ ,1]
            fity <- sdlist1[[i]][ ,2]
            fitx2 <- sdlist2[[i]][ ,1]
            fity2 <- sdlist2[[i]][ ,2]
        }
        x <- loess.smooth(fitx, fity,span = loessspan, degree = 1,
                              family = "symmetric", evaluation = 1000)
        if(i == 1) {
            plot(x$x, x$y, type = 'l', lwd = lwd, col = col[i],
                 lty = lty[i], xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab, ...)
        }else {
            lines(x$x, x$y, lwd = lwd, col = col[i], lty = lty[i])
        }
        if(!(plotcell %in% 1:2)){
            y <- loess.smooth(fitx2, fity2, span = loessspan, degree = 1,
                              family = "symmetric", evaluation = 1000)
            lines(y$x, y$y, lwd = lwd, col = col[i], lty = lty[i] + 2)
        }
    }
    box()
    if(plotcell %in% 1:2){
        legend('topright', names(dat@quantData), lwd = lwd, col = col,
               lty = lty, bty = "n", cex = cex.leg)
    }else{
        cells <- levels(dat@condInfo)
        legend('topright', c(names(dat@quantData), cells),
               lwd = lwd, col = c(col, rep("black", length(cells))),
               lty = c(lty, lty[1], lty[1] + 2), bty = "n", cex = cex.leg)
    }
    cell1 <- sapply(sdlist1,function(x){
        idx1 <- x[,1] <= thresholds[1] & x[,1] > log2(constant+0.1)
        idx2 <- x[,1] < thresholds[2] & x[,1] > thresholds[1]
        idx3 <- x[,1] >= thresholds[2]
        round(c(median(x[idx1, 2]), median(x[idx2, 2]), median(x[idx3, 2])),3)
       })
    cell2 <- sapply(sdlist2,function(x){
        idx1 <- x[,1] <= thresholds[1] & x[,1] > log2(constant+0.1)
        idx2 <- x[,1] < thresholds[2] & x[,1] > thresholds[1]
        idx3 <- x[,1] >= thresholds[2]
        round(c(median(x[idx1, 2]), median(x[idx2, 2]), median(x[idx3, 2])),3)
       })
    rownames(cell1) <- rownames(cell2) <- c(paste0("A<=", thresholds[1]),
                       paste0(thresholds[1], "<A<", thresholds[2]),
                       paste0("A>=", thresholds[2]))
    SD <- list(cell1 = cell1, cell2 = cell2)
    names(SD) <- levels(dat@condInfo)
    return(SD)
}
