#' @title Estimate And Plot Differential Expression
#'
#' @description For each pipeline, differential expression is
#' first estimated by fold change on 1 vs. 1 comparison between
#' cell lines. ROC curves then are made by comparing fold changes
#' with predefined true differentials. Then, ROC curves from multiple
#' 1 vs. 1 comparisons are averaged using threshold averaging
#' strategy. Standardized partial area under the curve (pAUC) is
#' reported for each pipeline.
#'
#' @param dat A \code{rnaseqcomp} S4 class object.
#' @param positive A logical vector with length equivalent to row
#' number of matrices in \code{dat@quantData}. \code{TRUE} means true
#' differential and \code{FALSE} means true non-differential, while
#' missing value \code{NA} means unknown.
#' @param fcsign A numeric vector with length equivalent to row
#' number of matrices in \code{dat@quantData}. Only values {1, -1, 0}
#' are allowed. 1 means upregulated in second cell line, -1 means
#' downregulated in second cell line, and 0 means no change. If elements
#' in \code{fcsign} correspond to \code{NA} in \code{positive}, these
#' elements will be ignored in estimation.
#' @param cut A numeric cutoff used to decide if fold change should be
#' estimated. For a 1 vs 1 comparison, if features have signals less than
#' \code{cut} in both samples, their fold changes will be set to 0.
#' (default: 1)
#' @param constant A numeric constant that is added to
#' quantifications before fold changes calculation. (default: 0.5)
#' @param thresholds A numeric vector defining cutoffs on fold changes
#' as the points to make threshold averaging on ROC curves.
#' (default: seq(12, 0, len = 300))
#' @param arrow A logical indicating if error bars should be added to
#' the averaged ROC curves. (default: FALSE)
#' @param ... Parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{plot}{ROC plots for all the quantification pipelines.}
#' \item{pAUC}{A numeric vector indicating pipeline accuracy.
#' This is standardized partial AUC based on ranges chosen on false
#' positive rate.}
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
#' plotROC(dat,simdata$meta$positive,simdata$meta$fcsign)



plotROC <- function(dat, positive, fcsign, cut = 1, constant = 0.5,
                    thresholds = seq(12, 0, len = 300), arrow = FALSE,
                    ...){
    if(!is(dat, 'rnaseqcomp'))
        stop('"plotSD" only plots class "rnaseqcomp".')
    dat@quantData <- lapply(dat@quantData, function(x) x + constant)
    cdList <- list()
    for(i in 1:2){
        cdList[[i]] <- lapply(dat@quantData, function(x)
                   log2(x[, dat@condInfo == levels(dat@condInfo)[i], drop=F]))
    }
    log2cut <- log2(constant + cut)
    fcList <- lapply(seq_along(dat@quantData), function(i){
        tmp <- matrix(NA, nrow(dat@quantData[[1]]),
                      ncol(cdList[[1]][[i]]) * ncol(cdList[[2]][[i]]))
        for(j in seq_len(ncol(cdList[[1]][[i]])))
            for(k in seq_len(ncol(cdList[[2]][[i]]))){
                idxcol <- (j - 1) * ncol(cdList[[2]][[i]]) + k
                tmp[, idxcol] <- cdList[[2]][[i]][, k] - cdList[[1]][[i]][, j]
                idxrow <- cdList[[1]][[i]][, j] < log2cut &
                    cdList[[2]][[i]][, k] < log2cut
                tmp[idxrow, idxcol] <- 0
            }
        tmp
    })
    fcList1 <- lapply(fcList, function(x){
        x <- x[!is.na(positive), ]
        x
    })
    positivesub <- positive[!is.na(positive)]
    fcsignsub <- fcsign[!is.na(positive)]
    getroc <- function(x, p){
        p2 <- sign(x) == sign(p) * abs(p)
        o <- order(abs(x), decreasing = TRUE)
        fp <- cumsum(!p2[o])
        tp <- cumsum(p2[o])
        fn <- sum(abs(p)) - cumsum(abs(p)[o])
        tn <- sum(p == 0) - cumsum(p[o]==0)
        fpr <- fp / (fp + tn)
        tpr <- tp / (tp + fn)
        cbind(tpr = tpr, fpr = fpr, fc = abs(x)[o])
    }
    proplist <- lapply(seq_along(fcList1), function(i){
        lapply(seq_len(ncol(fcList1[[i]])), function(j)
               getroc(fcList1[[i]][, j], positivesub * fcsignsub))
    })
    fprs <- tprs <- sdfprs <- sdtprs <- list()
    for(i in seq_along(proplist)){
        tprs[[i]] <- fprs[[i]] <- sdfprs[[i]] <- sdtprs[[i]] <-
            array(0, dim=length(thresholds))
        for(t in seq_along(thresholds)){
            fpr <- tpr <- c()
            threshold <- thresholds[t]
            for(j in seq_along(proplist[[i]])){
                idx <- sum(proplist[[i]][[j]][,3] >= threshold)
                if(idx==0){
                    tpr <- c(tpr, 0)
                    fpr <- c(fpr, 0)
                }else{
                    tpr <- c(tpr, proplist[[i]][[j]][idx, 1])
                    fpr<- c(fpr, proplist[[i]][[j]][idx, 2])
                }
            }
            tprs[[i]][t] <- median(tpr)
            fprs[[i]][t] <- median(fpr)
            sdtprs[[i]][t] <- sd(tpr)
            sdfprs[[i]][t] <- sd(fpr)
        }
    }
    if(!exists('xlab'))  xlab <- 'FP'
    if(!exists('ylab'))  ylab <- 'TP'
    if(!exists('xlim')) xlim <- c(0, 0.2)
    if(!exists('ylim')) ylim <- c(0, 1)
    if(!exists('lty')) lty <- 1
    if(!exists('lwd')) lwd <- 2
    if(is.function(col)) {
        if(length(dat@quantData)<3)
            col <- c("blue","orange")[seq_along(dat@quantData)]
        else {
            col <- brewer.pal(min(length(dat@quantData), 8), "Set2")
        }
    }
    lty <- rep_len(lty, length(dat@quantData))
    col <- rep_len(col, length(dat@quantData))
    for(i in seq_len(length(proplist))){
        x <- fprs[[i]]
        y <- tprs[[i]]
        sey <- sdtprs[[i]] / sqrt(tprs[[i]])
        sex <- sdfprs[[i]] / sqrt(fprs[[i]])
        if(i == 1) {
            plot(x, y, type = 'l', lwd = lwd, col = col[i],
                 lty = lty[i], xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab)
        }else {
            lines(x,y, lwd = lwd, col = col[i], lty = lty[i])
        }
        if(arrow){
            idx <- 2:(length(thresholds) - 1)
            arrows(x[idx], (y-sey)[idx], x[idx], (y+sey)[idx],
                   length = 0.02, angle = 90, code = 3)
            arrows((x-sex)[idx], y[idx], (x+sex)[idx], y[idx],
                   length = 0.02, angle = 90, code = 3)
        }
    }
    abline(a = 0,b = 1,lty = 2)
    legend('topleft', names(dat@quantData), lwd = lwd, col = col,
           lty = lty, cex = 1, bty = "n")
    AUC <- sapply(seq_along(tprs), function(i){
        tpr <- tprs[[i]]
        fpr <- fprs[[i]]
        auc <- 0
        J <- sum(fpr <= xlim[2])
        if(fpr[J] < xlim[2]){
            tpr[J+1] <- tpr[J] +
                (tpr[J+1] - tpr[J]) / (fpr[J+1] - fpr[J]) * (xlim[2] - fpr[J])
            fpr[J+1] <- xlim[2]
            J <- J + 1
        }
        for(j in seq_len(J-1))
            {
            auc <- auc + (fpr[j+1] - fpr[j]) * (tpr[j+1] + tpr[j]) / 2
        }
        auc
    })
    pAUC <- ((AUC - xlim[2]^2 / 2) / (xlim[2] - xlim[2]^2 / 2) + 1) / 2
    names(pAUC) <- names(dat@quantData)
    return(round(pAUC, 3))
}
