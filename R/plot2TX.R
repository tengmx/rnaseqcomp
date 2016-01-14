#' @title Estimate And Plot Transcript Proportion Difference
#'
#' @description For any compared two replicates in each cell line,
#' the proportion of one transcript for genes that only include two
#' annotated transcripts can be different even flipped. This function
#' estimates and plots the proportion difference stratefied by detrended
#' logsignal. Means of absolute difference will be reported for three
#' levels of detrened logsignals. Average is used when multiple
#' two-replicate comparisons included.
#'
#' @param dat A \code{rnaseqcomp} S4 class object.
#' @param genes A vector of gene names corresponding to quantified
#' transcripts. Note that \code{length(genes)} should equal to
#' \code{nrow(dat@quantData[[1]])}.
#' @param step A number specifying the resolution on detrended logsignal
#' for calculation and plotting the proportion difference.
#' (default: 0.5)
#' @param thresholds A vector of two numbers define cutoffs for
#' three levels of detreded log signals, where one number summary
#' will be generated. (default: c(1, 6))
#' @param plotcell 1 or 2 indicating which cell line
#' will be plotted. If values other than 1 and 2, both cell
#' lines will be plotted.  This value won't affect estimation for both
#' cell lines. (default: 1)
#' @param ... Parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{plot}{2TX plots of quantification pipelines for
#' selected cell line by \code{plotcell}.}
#' \item{2TX}{A matrix of mean proportion difference. Valuesa are based
#' on averageing two cell lines.}
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
#' plot2TX(dat,genes=simdata$meta$gene)


plot2TX <- function(dat, genes, step = 0.5, thresholds = c(1, 6), plotcell = 1,
                    ...){
    if(!is(dat,'rnaseqcomp'))
        stop('"plot2TX" only plots class "rnaseqcomp".')
    cdList <- list()
    for(i in 1:2){
        cdList[[i]] <- lapply(dat@quantData, function(x)
                        x[, dat@condInfo == levels(dat@condInfo)[i], drop=F])
    }
    tx2idx <- sapply(split(rownames(dat@quantData[[1]]), genes),
                     function(x) length(x) == 2)
    M <- A <- list()
    for(i in seq_along(dat@quantData)){
        prop <- sig <- list()
        for(k in seq_along(cdList)){
            prop[[k]] <- sig[[k]] <- list()
            for(j in seq_len(ncol(cdList[[k]][[i]]))){
                cell.rep <- split(cdList[[k]][[i]][, j], genes)[tx2idx]
                cell.rep.percent <- lapply(cell.rep, function(x){
                    if(sum(x, na.rm = T) == 0) x
                    else x / sum(x, na.rm = T)
                })
                prop[[k]][[j]] <- sapply(cell.rep.percent, function(x) x[1])
                sig[[k]][[j]] <- sapply(cell.rep, function(x) x[1])
            }
        }
        M[[i]] <- A[[i]] <- list()
        for(k in seq_along(cdList)){
            tmp1 <- tmp2 <- c()
            for(m in 1:(length(prop[[k]]) - 1)){
                for(n in (m+1):length(prop[[k]])){
                    tmp1 <- c(tmp1, prop[[k]][[m]] - prop[[k]][[n]])
                    tmp2 <- rbind(tmp2, cbind(sig[[k]][[m]], sig[[k]][[n]]))
                }
            }
            M[[i]][[k]] <- tmp1
            A[[i]][[k]] <- tmp2
        }
    }
    if(!exists('xlab'))  xlab <- 'Detrended logSignal'
    if(!exists('ylab'))  ylab <- 'Mean difference of transcript proportions'
    if(!exists('xlim'))  xlim <- c(0,12)
    if(!exists('ylim'))  ylim <- c(0,1)
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
    steps <- seq(xlim[1], xlim[2], step)
    pnelist <- list()
    for(k in 1:2){
        pnelist[[k]] <- lapply(seq_len(length(M)), function(i){
            sapply(seq_along(steps), function(j){
                if(j==1){
                    idx <- rowMeans(A[[i]][[k]]) <= 2^steps[j] &
                        rowMeans(A[[i]][[k]]) != 0
                }else{
                    idx <- rowMeans(A[[i]][[k]]) <= 2^steps[j] &
                        rowMeans(A[[i]][[k]]) > 2^steps[j-1]
                }
                if(sum(idx)==0) 0
                else mean(abs(M[[i]][[k]][idx]))
            })})
    }
    for(i in seq_along(dat@quantData)){
        if(plotcell == 1){
            ploty <- pnelist[[1]][[i]]
        }else if(plotcell == 2){
            ploty <- pnelist[[2]][[i]]
        }else{
            ploty <- pnelist[[1]][[i]]
            ploty2 <- pnelist[[2]][[i]]
        }
        if(i == 1) {
            plot(steps, ploty, type = 'o', lwd = lwd, col = col[i],
                 lty = lty[i], xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab)
        }else {
            lines(steps, ploty, lwd = lwd, col = col[i], lty = lty[i],
                  type = 'o')
        }
        if(!(plotcell %in% 1:2)){
            lines(steps, ploty2, lwd = lwd, col = col[i], lty = lty[i] + 2,
                  type = 'o')
        }
    }
    box()
    if(plotcell %in% 1:2){
        legend('topright', names(dat@quantData),
               lwd = lwd, col = col, lty = lty, bty = "n",cex = 1)
    }else{
        cells <- levels(dat@condInfo)
        legend('topright', c(names(dat@quantData), cells),
               lwd = lwd, col = c(col, rep("black", length(cells))),
               lty = c(lty, lty[1], lty[1] + 2), bty = "n", cex = 1)
    }
    TX2s <- list()
    for(k in 1:2){
        TX2s[[k]] <- sapply(seq_len(length(M)), function(i){
            idx1 <- rowMeans(A[[i]][[k]]) <= 2^thresholds[1] &
                rowMeans(A[[i]][[k]]) != 0
            idx2 <- rowMeans(A[[i]][[k]]) < 2^thresholds[2] &
                rowMeans(A[[i]][[k]]) > 2^thresholds[1]
            idx3 <- rowMeans(A[[i]][[k]]) >= 2^thresholds[2]
            c(mean(abs(M[[i]][[k]][idx1])), mean(abs(M[[i]][[k]][idx2])),
              mean(abs(M[[i]][[k]][idx3])))
        })
    }
    TX2 <- (TX2s[[1]] + TX2s[[2]])/2
    colnames(TX2) <- names(dat@quantData)
    rownames(TX2) <- c(paste0("A<=", thresholds[1]),
                       paste0(thresholds[1], "<A<", thresholds[2]),
                       paste0("A>=", thresholds[2]))
    return(round(TX2, 3))
}
