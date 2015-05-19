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
#' @param microarray A numeric vector of fold change by
#' microarray or other 'gold standard', with each elements
#' corresponding to rows of \code{quantData} slot in \code{dat1}
#' or \code{dat2}. Missing data NA allowed. (default: NULL)
#' @param step Plot steps on x-axis (default: 5).
#' @param type Plot types (default: 'l').
#' @param lwd Plot line weights (default: 2).
#' @param col Plot colors (default: NULL, colors are assigned
#' by package \code{RColorBrewer}).
#' @param xlim Plot limits of x-axis (default: c(20, 2000)).
#' @param ylim Plot limits of y-axis (default: c(0, 1)).
#' @param xlab Plot label of x-axis
#' (default: 'Size of List').
#' @param ylab Plot label of y-axis
#' (default: 'Proportion in Common').
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
#' txFIdx <- genemeta$type == "protein_coding"
#' hkIdx <- genemeta$housekeeping
#' unitFIdx <- grepl("Cufflinks",repInfo)
#' dat1 <- matrixFilter(gm12878,repInfo,txFIdx,hkIdx,unitFIdx)
#' dat2 <- matrixFilter(k562,repInfo,txFIdx,hkIdx,unitFIdx)
#' 
#' plotCAT(dat1,dat2)
#' plotCAT(dat1,dat2,constant=1)
#'
#' genes <- genemeta[genemeta$type == "protein_coding", 1]
#' microarray <- arrayFC[match(genes,names(arrayFC))]
#' plotCAT(dat1,dat2,microarray=microarray)
#' plotCAT(dat1,dat2,constant=1,microarray=microarray)

plotCAT <- function(dat1, dat2, constant = NULL, microarray = NULL, step = 5L,
                    type = 'l', lwd = 2, col = NULL, xlim = c(20, 2000),
                    ylim = c(0, 1), xlab = "Size of List",
                    ylab = "Proportion in Common",
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
    fcList <- lapply(seq_len(length(cdList1)), function(i)
                     log2(cdList1[[i]]) - log2(cdList2[[i]]))
    # handling 0s if no constant 
    if(is.null(constant)){
        fcList <- lapply(fcList, function(x) {
                         x[is.nan(x) | is.infinite(x)] <- 0
                         x })
    }
    # mocroarray or not
    if(is.null(microarray)){
        fcList <- lapply(fcList, abs)
    }else{
        fcList <- lapply(fcList, function(x)
                         abs(cbind(rowMeans(x), microarray)))
    }
    # missing data
    fcList <- lapply(fcList, function(x){ 
        x[is.na(rowSums(x)), ] <- 0
        x
    })
    if(is.null(xlim)) xlim <- c(20, 2000)
    if(is.null(ylim)) ylim <- c(0, 1)
    # ranks and proportion of concordance
    ranks <- seq(max(xlim[1], 1), min(xlim[2], nrow(dat1@quantData)), step)
    proplist <- lapply(fcList, function(x){
        orders <- cbind(order(x[ ,1], decreasing = TRUE),
                        order(x[ ,2], decreasing = TRUE))
        sapply(ranks, function(i) sum(table(orders[seq_len(i), ]) == 2) / i)
    })
    names(proplist) <- levels(repInfo)
    if(is.null(xlab))  xlab <- 'Size of List'
    if(is.null(ylab))  ylab <- 'Proportion in Common'
    if(is.null(col))   col <- brewer.pal(min(length(proplist), 9), "Set1")
    col <- rep_len(col, length(proplist))
    type <- rep_len(type, length(proplist))
    lwd <- rep_len(lwd, length(proplist))
    for(i in seq_len(length(proplist))){
        if(i == 1) {
            plot(ranks, proplist[[i]], type = type[i], lwd = lwd[i],
                 col = col[i], ylim = ylim, xlim = xlim, xlab = xlab,
                 ylab = ylab, ...)
        }else {
            points(ranks, proplist[[i]], type = type[i], lwd = lwd[i],
                   col = col[i])
        }
    }
    legend('bottomright', names(proplist), lwd = lwd, col = col)
    sapply(proplist, median)
}
################################################################################
