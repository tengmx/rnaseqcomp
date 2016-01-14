#' @title Estimate And Plot Express And Non-express Features
#'
#' @description For each cell line, any compared two replicates
#' might have a portion of transcripts that express in one replicate
#' but not the other, depending on what cutoff is used to define
#' non-express. This function estimate and plot the proportion of
#' disagreement using multiple cutoffs. Average is used when multiple
#' two-replicate comparisons included.
#'
#' @param dat A \code{rnaseqcomp} S4 class object.
#' @param steps A numeric vector specifying log-scale cutoffs to be used
#' for calculation and plotting. (default: seq(-0.5, 12, 0.5))
#' @param Ks A numeric vector specifying which cutoffs to be highlighted
#' and to which the reported proportions to be corresponding.
#' @param pchK Plot styles of highlight points corresponding
#' to \code{Ks}. (default: seq_along(Ks) - 1)
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
#' @param xlab Plot label of x-axis.
#' (default: 'Proportion of detrended logSignal below K')
#' @param ylab Plot label of y-axis.
#' (default: 'Proportion of disagreement between replicates')
#' @param ... Other parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{plot}{NE plots of quantification pipelines for
#' selected cell line by \code{plotcell}.}
#' \item{NE}{A list of sub-list with each sub-list corresponding to
#' one cell line. Each element of sub-list is one 2*n matrix for
#' one pipeline, where \code{n=length(Ks)}. The first row of
#' matrix gives the proportion
#' of disagreement and the second row gives the proportion of both
#' replicates under (non-express) correspoinding cutoff \code{Ks}.}
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
#' plotNE(dat,col=c("blue","orange"))

plotNE <- function(dat,  steps = seq(-0.5, 12, 0.5), Ks = 0:3,
                   pchK = seq_along(Ks) - 1, plotcell = 1,
                   lwd = 2, col = NULL, lty = 1, cex.leg = 1,
                   xlim = NULL, ylim = NULL,
                   xlab = "% of detrended logSignal below K",
                   ylab = "% of disagreement between replicates",
                   ...){
    if(!is(dat,'rnaseqcomp'))
        stop('"plotNE" only plots class "rnaseqcomp".')
    cdList <- list()
    for(i in 1:2){
        cdList[[i]] <- lapply(dat@quantData, function(x)
                 log2(x[, dat@condInfo == levels(dat@condInfo)[i], drop=F]))
    }
    steps <- sort(unique(c(steps, Ks)))
    pnelist1 <- lapply(seq_along(dat@quantData), function(j){
        props <- matrix(0, 2, length(steps))
        count <- 0
        for(m in 1:(ncol(cdList[[1]][[j]])-1))
            for(n in (m+1):ncol(cdList[[1]][[j]])){
                props <- props + sapply(steps, function(K){
                    rep1 <- cdList[[1]][[j]][, m] < K
                    rep2 <- cdList[[1]][[j]][, n] < K
                    c(sum((!rep1 & rep2) | (rep1 & !rep2)) / length(rep1),
                      sum(rep1 & rep2)/length(rep1))})
                count <- count + 1
            }
        props/count
    })
    pnelist2 <- lapply(seq_along(dat@quantData),function(j){
        props <- matrix(0, 2, length(steps))
        count <- 0
        for(m in 1:(ncol(cdList[[2]][[j]])-1))
            for(n in (m+1):ncol(cdList[[2]][[j]])){
                props <- props + sapply(steps, function(K){
                    rep1 <- cdList[[2]][[j]][, m] < K
                    rep2 <- cdList[[2]][[j]][, n] < K
                    c(sum((!rep1 & rep2) | (rep1 & !rep2)) / length(rep1),
                      sum(rep1 & rep2)/length(rep1))})
                count <- count + 1
            }
        props/count
    })
    if(is.null(xlim)) xlim <- c(0.3, 1)
    if(is.null(ylim)) ylim <- c(0, 0.2)
    names(pnelist1) <- names(pnelist2) <- names(dat@quantData)
    if(is.null(col))   col <- brewer.pal(min(length(pnelist1), 8), "Set2")
    col <- rep_len(col, length(pnelist1))
    lty <- rep_len(lty, length(pnelist1))
    idx <- match(Ks, steps)
    for(i in seq_along(pnelist1)){
        if(plotcell==1){
             plotx <- pnelist1[[i]][2,]
             ploty <- pnelist1[[i]][1,]
         }else if(plotcell==2){
             plotx <- pnelist2[[i]][2,]
             ploty <- pnelist2[[i]][1,]
         }else{
             plotx <- pnelist1[[i]][2,]
             ploty <- pnelist1[[i]][1,]
             plotx2 <- pnelist2[[i]][2,]
             ploty2 <- pnelist2[[i]][1,]
         }
        if(i == 1) {
            plot(plotx, ploty, type = 'l', lwd = lwd, col = col[i],
                 lty = lty[i], xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab, ...)
        }else {
            lines(plotx, ploty, lwd = lwd, col = col[i], lty = lty[i])
        }
        points(plotx[idx], ploty[idx], pch = pchK, cex = 0.6)
        if(!(plotcell %in% 1:2)){
            lines(plotx2, ploty2, lwd = lwd, col = col[i], lty = lty[i] + 2)
            points(plotx2[idx], ploty2[idx], pch = pchK, cex = 0.6)
        }
    }
    box()
    if(plotcell %in% 1:2){
        legend('topright', names(dat@quantData),
               lwd = lwd, col = col, lty = lty, bty = "n",cex = cex.leg)
    }else{
        cells <- levels(dat@condInfo)
        legend('topright', c(names(dat@quantData), cells),
               lwd = lwd, col = c(col, rep("black", length(cells))),
               lty = c(lty, lty[1], lty[1] + 2), bty = "n", cex = cex.leg)
    }
    legend('bottomleft', paste0("K=", Ks), pch = pchK,
           bty = "n", cex = cex.leg * 0.8)
    cell1 <- lapply(pnelist1, function(x) {
        tmp <- x[, idx]
        colnames(tmp) <- Ks
        rownames(tmp) <- c("y", "x")
        tmp
    })
    cell2 <- lapply(pnelist2,function(x) {
        tmp <- x[, idx]
        colnames(tmp) <- Ks
        rownames(tmp) <- c("y", "x")
        tmp
    })
    NE <- list(cell1 = cell1, cell2 = cell2)
    names(NE) <- levels(dat@condInfo)
    return(NE)
}
