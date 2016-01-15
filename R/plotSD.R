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
#' @param plotcell 1 or 2 indicating which cell line
#' will be plotted. If values other than 1 and 2, both cell
#' lines will be plotted.  This value won't affect estimation for both
#' cell lines. (default: 1)
#' @param ... Parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#' @import methods
#'
#' @return
#' \item{plot}{SD plots of quantification pipelines for
#' selected cell line by \code{plotcell}.}
#' \item{SD}{One  matrix of median standard deviations.}
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
                   ...){
    if(!is(dat, 'rnaseqcomp'))
        stop('"plotSD" only plots class "rnaseqcomp".')
    dat@quantData <- lapply(dat@quantData, function(x) x + constant)
    sdlist <- list()
    for(i in 1:2){
        tmpsig <- lapply(dat@quantData, function(x)
                        x[,dat@condInfo == levels(dat@condInfo)[i], drop=F])
        sdlist[[i]] <- lapply(seq_along(dat@quantData), function(j){
            tmp <- tmpsig[[j]]
            idx <- which(apply(tmp, 1, max) > constant)
            tmp <- log2(tmp[idx,])
            cbind(rowMeans(tmp), apply(tmp, 1, sd))
        })

    }
    para <- list(...)
    if(!('xlab' %in% names(para)))  xlab <- 'Detrended logSignal'
    else xlab <- para$xlab
    if(!('ylab' %in% names(para)))  ylab <- 'SD'
    else ylab <- para$ylab
    if(!('xlim' %in% names(para)))  xlim <- c(-1, 12)
    else xlim <- para$xlim
    if(!('ylim' %in% names(para)))  ylim <- c(0, 1.5)
    else ylim <- para$ylim
    if(!('lty' %in% names(para))) lty <- 1
    else lty <- para$lty
    if(!('lwd' %in% names(para))) lwd <- 2
    else lwd <- para$lwd
    if(!('col' %in% names(para))) {
        if(length(dat@quantData)<3)
            col <- c("blue","orange")[seq_along(dat@quantData)]
        else {
            col <- brewer.pal(min(length(dat@quantData), 8), "Set2")
        }
    }else col <- para$col
    lty <- rep_len(lty, length(dat@quantData))
    col <- rep_len(col, length(dat@quantData))
    for(i in seq_along(dat@quantData)){
        if(plotcell==1){
            fitx <- sdlist[[1]][[i]][ ,1]
            fity <- sdlist[[1]][[i]][ ,2]
        }else if(plotcell==2){
            fitx <- sdlist[[2]][[i]][ ,1]
            fity <- sdlist[[2]][[i]][ ,2]
        }else{
            fitx <- sdlist[[1]][[i]][ ,1]
            fity <- sdlist[[1]][[i]][ ,2]
            fitx2 <- sdlist[[2]][[i]][ ,1]
            fity2 <- sdlist[[2]][[i]][ ,2]
        }
        x <- loess.smooth(fitx, fity,span = loessspan, degree = 1,
                          family = "symmetric", evaluation = 1000)
        if(i == 1) {
            plot(x$x, x$y, type = 'l', lwd = lwd, col = col[i],
                 lty = lty[i], xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab)
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
               lty = lty, bty = "n", cex = 1)
    }else{
        cells <- levels(dat@condInfo)
        legend('topright', c(names(dat@quantData), cells),
               lwd = lwd, col = c(col, rep("black", length(cells))),
               lty = c(lty, lty[1], lty[1] + 2), bty = "n", cex = 1)
    }
    SDs <- lapply(sdlist, function(y)
                  sapply(y, function(x){
                      idx1 <- x[,1] <= thresholds[1] &
                          x[,1] > log2(constant+0.1)
                      idx2 <- x[,1] < thresholds[2] &
                          x[,1] > thresholds[1]
                      idx3 <- x[,1] >= thresholds[2]
                      c(median(x[idx1, 2]), median(x[idx2, 2]),
                        median(x[idx3, 2]))
                  }))
    SD <- sqrt((SDs[[1]]^2 + SDs[[2]]^2)/2)
    colnames(SD) <- names(dat@quantData)
    rownames(SD) <- c(paste0("A<=", thresholds[1]),
                      paste0(thresholds[1], "<A<", thresholds[2]),
                      paste0("A>=", thresholds[2]))
    return(round(SD, 3))
}
