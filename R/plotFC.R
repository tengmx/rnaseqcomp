#' @title Estimate And Plot Fold Change Accuracy
#'
#' @description For each pipeline, differential expression is
#' estimated by fold change on mean signals across replicates of
#' cell lines. For features that are truely differential
#' expressed, their fold changes levels are summarized based on
#' different levels of detrended logsignals.
#'
#' @param dat A \code{rnaseqcomp} S4 class object.
#' @param positive A logical vector with length equivalent to row
#' number of matrices in \code{dat@quantData}. \code{TRUE} means true
#' differential and \code{FALSE} means true non-differential, while
#' missing value \code{NA} means unknown.
#' @param fcsign A numeric vector with length equivalent to row
#' number of matrices in \code{dat@quantData}. Only values {1, -1, 0, NA}
#' are allowed. 1 means upregulated in second cell line, -1 means
#' downregulated in second cell line, and 0 means no change. If elements
#' in \code{fcsign} is NA or correspond to \code{NA} in \code{positive},
#' these elements will be ignored in estimation.
#' @param constant A numeric constant that is added to
#' quantifications before fold changes calculation. (default: 0.5)
#' @param loessspan A numeric number indicating span used
#' for loess smooth. Details see
#' \code{loess.smooth} function. (Default: 1/3)
#' @param thresholds A numeric vector defining cutoffs on fold changes
#' as the points to make threshold averaging on ROC curves.
#' (default: seq(12, 0, len = 300))
#' @param ... Parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{plot}{Fold change plots for all the quantification pipelines.}
#' \item{FC}{A numeric vector indicating median fold changes in three
#' different levels of detrended logsignals.}
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
#' ## only select the true differential that have exact fold changes
#' simdata$meta$fcsign[simdata$meta$fcstatus == "off.on"] <- NA
#' plotFC(dat,simdata$meta$positive,simdata$meta$fcsign)

plotFC <- function(dat, positive, fcsign, constant = 0.5, loessspan=1/3,
                   thresholds = c(1, 6), ...){
    if(!is(dat, 'rnaseqcomp'))
        stop('"plotSD" only plots class "rnaseqcomp".')
    dat@quantData <- lapply(dat@quantData,function(x) x + constant)
    cdList <- list()
    for(i in 1:2){
        cdList[[i]] <- lapply(dat@quantData, function(x)
                   rowMeans(log2(x[, dat@condInfo ==
                                   levels(dat@condInfo)[i], drop=F])))
    }
    fcList <- lapply(seq_along(dat@quantData), function(i){
                cbind((cdList[[2]][[i]] + cdList[[1]][[i]])/2,
                      (cdList[[2]][[i]] - cdList[[1]][[i]]) * fcsign)
            })
    fcList1 <- lapply(fcList, function(x){
        x <- x[which(positive & !is.na(fcsign)), ]
        x
    })
    para <- list(...)
    if(!('xlab' %in% names(para)))  xlab <- 'Detrended logSignal'
    else xlab <- para$xlab
    if(!('ylab' %in% names(para)))  ylab <- 'log2FoldChange'
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
    for(i in seq_along(fcList1)){
        loessfit <- loess.smooth(fcList1[[i]][,1], fcList1[[i]][,2],
                                 span = loessspan, degree = 1,
                                 family = "symmetric", evaluation = 1000)
        if(i == 1) {
            plot(loessfit$x, loessfit$y, type = 'l', lwd = lwd,
                 col = col[i], ylim = ylim, xlim = xlim, xlab = xlab,
                 ylab = ylab, lty = lty[i])
        }else {
            lines(loessfit$x, loessfit$y, lwd = lwd,
                   col = col[i], lty = lty[i])
        }
    }
    legend("bottomright", names(dat@quantData), lwd = lwd, col = col,
           lty = lty, cex = 1, bty = "n")
    FC <- sapply(fcList1,function(x){
        idx1 <- x[,1] <= thresholds[1] & x[,1] > log2(constant+0.1)
        idx2 <- x[,1] < thresholds[2] & x[,1] > thresholds[1]
        idx3 <- x[,1] >= thresholds[2]
        c(median(x[idx1, 2]), median(x[idx2, 2]), median(x[idx3, 2]))
    })
    colnames(FC) <- names(dat@quantData)
    rownames(FC) <- c(paste0("A<=",thresholds[1]),
                      paste0(thresholds[1],"<A<",thresholds[2]),
                      paste0("A>=",thresholds[2]))
    return(round(FC ,3))
}
