#' @title CAT Plots of Differential Express Features
#'
#' @description For each pipeline, 2 biological conditions
#' (cell lines) each with 2 quantification replicates are
#' compared and fold changes of each replicate are calculated.
#' Then, CAT plots between replicates (precision) or between
#' mean of replicates and microarray (accuracy).
#'
#' @param dat1,dat2 \code{rnaseqcomp} S4 class objects for two
#' conditions. \code{dat1} and \code{dat2} should have the same
#' size of slot \code{quantData} and the same \code{repInfo}.
#' @param constant A numeric constant that can be added to
#' quantifications before fold changes calculation
#' (default: NULL).
#' @param infinity A logical indicator that specify if fold change
#' of infinity should be consitered. Functional only if \code{constant}
#' is not a positive number. (default: FALSE)
#' @param microarray A numeric vector of fold change by
#' microarray or other 'gold standard', with each elements
#' corresponding to rows of \code{quantData} slot in \code{dat1}
#' or \code{dat2}. Missing data NA allowed. (default: NULL)
#' @param step Plot steps on x-axis (default: 5).
#' @param type Plot types (default: 'l').
#' @param lwd Plot line weights (default: 2).
#' @param col Plot colors (default: NULL, colors are assigned
#' by package \code{RColorBrewer}).
#' @param lty Plot line styles (default: 1).
#' @param xlim Plot limits of x-axis (default: c(20, 500)).
#' @param ylim Plot limits of y-axis (default: c(0, 1)).
#' @param xlab Plot label of x-axis (default: 'Size of List').
#' @param ylab Plot label of y-axis (default: 'Proportion in Common').
#' @param cex.leg Legend size (default: 0.6).
#' @param ... Other parameters for base function \code{plot}.
#'
#' @import RColorBrewer
#'
#' @return
#' \item{CAT plot}{CAT plots for all the quantification pipelines.}
#' \item{Precision or Accuaracy}{A numeric vector of pipeline
#' precision or accuracy, depending on whether microarray is
#' provided.}
#'
#' @export
#' @examples
#' data(encodeCells)
#' txFIdx <- encodeCells$genemeta$type == "protein_coding"
#' hkIdx <- encodeCells$genemeta$housekeeping
#' unitFIdx <- grepl("Cufflinks",encodeCells$repInfo)
#' dat1 <- matrixFilter(encodeCells$gm12878,encodeCells$repInfo,
#' txFIdx,hkIdx,unitFIdx)
#' dat2 <- matrixFilter(encodeCells$k562,encodeCells$repInfo,
#' txFIdx,hkIdx,unitFIdx)
#'
#' plotCAT(dat1,dat2)
#' plotCAT(dat1,dat2,constant=1)
#'
#' genes <- encodeCells$genemeta[encodeCells$genemeta$type ==
#' "protein_coding", 1]
#' microarray <- encodeCells$arrayFC[match(genes,names(encodeCells$arrayFC))]
#' plotCAT(dat2,dat1,microarray=microarray)
#' plotCAT(dat2,dat1,constant=1,microarray=microarray)

plotCAT <- function(dat1, dat2, constant = NULL, microarray = NULL,
                    infinity = FALSE, step = 5L,
                    type = 'l', lwd = 2, col = NULL, lty = 1,
                    xlim = c(20L, 500L), ylim = c(0, 1), xlab = "Size of List",
                    ylab = "Proportion in Common", cex.leg = 0.6,
                    ...){
    if(!is(dat1, 'rnaseqcomp') || !is(dat2, 'rnaseqcomp'))
        stop('"plotCAT" only plots class "rnaseqcomp".')
    if(!identical(dat1@repInfo, dat2@repInfo))
        stop('"repInfo" are not identical bwteen "dat1" and "dat2".')
    if(nrow(dat1@quantData) != nrow(dat2@quantData))
       stop('Rows are different between "dat1" and "dat2".')
    if(!is.null(constant) && (!is.numeric(constant) || constant < 0))
        stop('"constant" is not correct.')
    if(!is.null(microarray) && (!is.numeric(microarray) ||
                                length(microarray) != nrow(dat1@quantData)))
        stop('"microarray" should be numeric vector or NULL')
    if(!is.numeric(step) || step < 1)
        stop('"step" should be natual number.')
    if(!is.logical(infinity) || length(infinity) != 1)
        stop('"infinity" should be a logical.')
    repInfo <- dat1@repInfo
    if(!is.null(constant)){
        dat1@quantData <- dat1@quantData + constant
        dat2@quantData <- dat2@quantData + constant
    }
    cdList1 <- lapply(levels(repInfo), function(i)
                      dat1@quantData[ ,repInfo == i])
    cdList2 <- lapply(levels(repInfo), function(i)
                      dat2@quantData[ ,repInfo == i])
    # fold change
    if(is.null(microarray)){
        fcList <- lapply(seq_len(length(cdList1)), function(i)
                         log2(cdList1[[i]]) - log2(cdList2[[i]]))
    }else{
        fcList <- lapply(seq_len(length(cdList1)), function(i)
                         cbind(log2(rowMeans(cdList1[[i]])) -
                               log2(rowMeans(cdList2[[i]])),microarray))
    }
    # handling 0s if no constant or constant is 0
    if(is.null(constant) || constant == 0){
        if(!infinity){
            fcList <- lapply(fcList, function(x) {
                x[is.nan(x) | is.infinite(x)] <- 0
                x })
        }else if(is.null(microarray)){
            fcList <- lapply(seq_len(length(fcList)), function(i){
                x <- fcList[[i]]
                x[is.nan(x)] <- 0
                cutinf <- max(abs(x[!is.infinite(x)]),na.rm = TRUE)
                x[is.infinite(x) & x < 0] <-
                    -cutinf - cdList2[[i]][is.infinite(x) & x < 0]
                x[is.infinite(x) & x > 0] <-
                    cutinf + cdList1[[i]][is.infinite(x) & x > 0]
                x
            })
        }else{
            fcList <- lapply(seq_len(length(fcList)), function(i){
                x <- fcList[[i]]
                x[is.nan(x[,1]),1] <- 0
                cutinf <- max(abs(x[!is.infinite(x[,1]),1]), na.rm = TRUE)
                x[is.infinite(x[,1]) & x[,1] < 0, 1] <- -cutinf -
                    rowMeans(cdList2[[i]])[is.infinite(x[,1]) & x[,1] < 0]
                x[is.infinite(x[,1]) & x[,1] > 0, 1] <- cutinf +
                    rowMeans(cdList1[[i]])[is.infinite(x[,1]) & x[,1] > 0]
                x
            })
        }
    }
    # missing data
    fcList <- lapply(fcList, function(x){
        x[is.na(rowSums(x)), ] <- 0
        x
    })
    fcListAbs <- lapply(fcList,abs)
    if(is.null(xlim)) xlim <- c(20L, 500L)
    if(is.null(ylim)) ylim <- c(0, 1)
    # ranks and proportion of concordance
    ranks <- seq(max(xlim[1], 1L), min(xlim[2], nrow(dat1@quantData)), step)
    proplist <- lapply(seq_len(length(fcListAbs)), function(i){
        x <- fcListAbs[[i]]
        y <- fcList[[i]]
        order1 <- order(x[ ,1], decreasing = TRUE)
        order1[y[order1,1] < 0] <- -order1[y[order1,1] < 0]
        order2 <- order(x[ ,2], decreasing = TRUE)
        order2[y[order2,2] < 0] <- -order2[y[order2,2] < 0]
        orders <- cbind(order1,order2)
        sapply(ranks, function(j) sum(table(orders[seq_len(j), ]) == 2) / j)
    })
    names(proplist) <- levels(repInfo)
    if(is.null(xlab))  xlab <- 'Size of List'
    if(is.null(ylab))  ylab <- 'Proportion in Common'
    if(is.null(col))   col <- brewer.pal(min(length(proplist), 9), "Set1")
    col <- rep_len(col, length(proplist))
    type <- rep_len(type, length(proplist))
    lwd <- rep_len(lwd, length(proplist))
    lty <- rep_len(lty, length(proplist))
    for(i in seq_len(length(proplist))){
        if(i == 1) {
            plot(ranks, proplist[[i]], type = type[i], lwd = lwd[i],
                 col = col[i], ylim = ylim, xlim = xlim, xlab = xlab,
                 ylab = ylab, lty = lty[i], ...)
        }else {
            points(ranks, proplist[[i]], type = type[i], lwd = lwd[i],
                   col = col[i], lty = lty[i])
        }
    }
    legend('bottomright', names(proplist), lwd = lwd, col = col,
           lty = lty, cex = cex.leg)
    sapply(proplist, median)
}
