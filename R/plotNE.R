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
#' @param plotcell 1 or 2 indicating which cell line
#' will be plotted. If values other than 1 and 2, both cell
#' lines will be plotted.  This value won't affect estimation for both
#' cell lines. (default: 1)
#' @param ... Parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{plot}{NE plots of quantification pipelines for
#' selected cell line by \code{plotcell}.}
#' \item{NE}{A list of two matrices. The first matrix gives
#' the proportion of disagreement and the second matrix gives the
#' proportion of both replicates under (non-express)
#' correspoinding cutoff \code{Ks}. Values are based on average
#' of two cell lines.}
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
#' plotNE(dat)

plotNE <- function(dat,  steps = seq(-0.5, 12, 0.5), Ks = 0:3,
                   pchK = seq_along(Ks) - 1, plotcell = 1,
                   ...){
    if(!is(dat,'rnaseqcomp'))
        stop('"plotNE" only plots class "rnaseqcomp".')
    para <- list(...)
    if(length(para)!=0 && any(!(names(para) %in%
             c("xlim","ylim","xlab","ylab","lty","lwd","main","col"))))
        stop('... contains non-used arguments.')
    pnelist <- list()
    steps <- sort(unique(c(steps, Ks)))
    for(i in 1:2){
        tmp <- lapply(dat@quantData, function(x)
                log2(x[, dat@condInfo == levels(dat@condInfo)[i], drop=F]))
        pnelist[[i]] <- lapply(seq_along(dat@quantData), function(j){
            props <- matrix(0, 2, length(steps))
            count <- 0
            for(m in 1:(ncol(tmp[[j]])-1))
                for(n in (m+1):ncol(tmp[[j]])){
                    props <- props + sapply(steps, function(K){
                        rep1 <- tmp[[j]][, m] < K
                        rep2 <- tmp[[j]][, n] < K
                        c(sum((!rep1 & rep2) | (rep1 & !rep2)) / length(rep1),
                          sum(rep1 & rep2)/length(rep1))})
                    count <- count + 1
                }
            props/count
        })
    }
    if(!('xlab' %in% names(para)))  xlab <- '% of detrended logSignal below K'
    else xlab <- para$xlab
    if(!('ylab' %in% names(para)))
        ylab <- '% of disagreement between replicates'
    else ylab <- para$ylab
    if(!('xlim' %in% names(para)))  xlim <- c(0.3, 1)
    else xlim <- para$xlim
    if(!('ylim' %in% names(para)))  ylim <- c(0, 0.2)
    else ylim <- para$ylim
    if(!('lty' %in% names(para))) lty <- 1
    else lty <- para$lty
    if(!('lwd' %in% names(para))) lwd <- 2
    else lwd <- para$lwd
    if(!('main' %in% names(para))) main <- "NE plot"
    else main <- para$main
    if(!('col' %in% names(para))) {
        if(length(dat@quantData)<3)
            col <- c("blue","orange")[seq_along(dat@quantData)]
        else {
            col <- brewer.pal(min(length(dat@quantData), 8), "Set2")
        }
    }else col <- para$col
    lty <- rep_len(lty, length(dat@quantData))
    col <- rep_len(col, length(dat@quantData))
    idx <- match(Ks, steps)
    for(i in seq_along(dat@quantData)){
        if(plotcell==1){
             plotx <- pnelist[[1]][[i]][2,]
             ploty <- pnelist[[1]][[i]][1,]
         }else if(plotcell==2){
             plotx <- pnelist[[2]][[i]][2,]
             ploty <- pnelist[[2]][[i]][1,]
         }else{
             plotx <- pnelist[[1]][[i]][2,]
             ploty <- pnelist[[1]][[i]][1,]
             plotx2 <- pnelist[[2]][[i]][2,]
             ploty2 <- pnelist[[2]][[i]][1,]
         }
        if(i == 1) {
            plot(plotx, ploty, type = 'l', lwd = lwd, col = col[i],
                 lty = lty[i], xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab, main = main)
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
               lwd = lwd, col = col, lty = lty, bty = "n", cex = 1)
    }else{
        cells <- levels(dat@condInfo)
        legend('topright', c(names(dat@quantData), cells),
               lwd = lwd, col = c(col, rep("black", length(cells))),
               lty = c(lty, lty[1], lty[1] + 2), bty = "n", cex = 1)
    }
    legend('bottomleft', paste0("K=", Ks), pch = pchK,
           bty = "n", cex = 0.8)
    NEs <-  lapply(pnelist, function(y)
                   lapply(y, function(x) {
                       tmp <- x[, idx]
                       colnames(tmp) <- Ks
                       rownames(tmp) <- c("y", "x")
                       tmp
                   }))
    NE <- (sapply(NEs[[1]],function(x) x[1,]) +
           sapply(NEs[[2]],function(x) x[1,])) / 2
    NN <- (sapply(NEs[[1]],function(x) x[2,]) +
           sapply(NEs[[2]],function(x) x[2,])) / 2
    colnames(NE) <- colnames(NN) <- names(dat@quantData)
    return(list(NE = round(NE, 3),NN = round(NN, 3)))
}
