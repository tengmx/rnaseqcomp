#' @title Quantification Table Filtering
#'
#' @description This is the function to do any pre-filtering or
#' pre-processing analysis for dowstreaming benchmark estimation
#' and visualization. Pre-filtering includes row selection (e.g.
#' protein coding genes) of quantification table ; pre-processing
#' inlcudes calculation on a set of rows as calibration reference
#' (e.g. house keeping genes) across different quantification
#' pipelines, estimation of confidence thresholds on all the
#' pipelins based on given thresholds from selected pipelines.
#'
#' @param quantData A matrix of quantifications with rows by
#' features (genes, transcriptsm, jucntions or exons) and columns
#' by pipelines each with 2 replicates. Missing data is allowed
#' as NA.
#' @param repInfo A factor documenting quantification pipeline
#' names correponding to the columns of \code{quantData}.
#' @param txFIdx A logical vector corresponding to the rows of
#' \code{quantData}, providing which features should be considered
#' for downsteam analyis, e.g. protein coding genes.
#' @param hkIdx A logical vector corresponding to the rows of
#' \code{quantData}, providing which features should be considered
#' as calibration reference, e.g. house keeping genes.
#' @param unitFIdx A logical vector corresponding to the columns
#' of \code{quantData}, providing which columns (i.e. pipelines)
#' should be applied for given threshold \code{unitCut}.
#' @param unitCut A numeric threshold for signal cutoff at given
#' quantification table columns by \code{unitFidx} (default: 1).
#' The effects cutoff are estimated for other columns accordingly.
#'
#' @return A \code{rnaseqcomp} S4 class object
#' \item{quantData}{A filtered table for downsteam analysis.}
#' \item{repInfo}{A factor documenting quantification pipeline
#' names correponding to the columns of \code{quantData}.}
#' \item{hkmed}{A numeric vector giving the log scale medians
#' of calibration reference.}
#' \item{scaler}{A number providing the scales for downstream
#' detrended signal.}
#'
#' @details
#' In the functions \code{plotMAD} and \code{plotNE}, detrended
#' signals with value 0 will be at the same level as value 1 for
#' giving pipelines by \code{unitFIdx}.
#'
#' @importFrom methods show
#' @export
#'
#' @examples
#' data(encodeCells)
#' txFIdx <- encodeCells$genemeta$type == "protein_coding"
#' hkIdx <- encodeCells$genemeta$housekeeping
#' unitFIdx <- grepl("Cufflinks",encodeCells$repInfo)
#' dat <- matrixFilter(encodeCells$gm12878,encodeCells$repInfo,
#' txFIdx,hkIdx,unitFIdx)

matrixFilter <- function(quantData, repInfo, txFIdx, hkIdx,
                                 unitFIdx, unitCut = 1){
    if(!is.factor(repInfo) || sum(summary(repInfo) != 2) > 0 )
        stop('"repInfo" is not factor or not fitting two replicates.')
    if(!is.matrix(quantData) || sum(quantData < 0, na.rm = TRUE) > 0)
        stop('"quantData" is not a non-negative matrix.')
    if(ncol(quantData) != length(repInfo))
        stop('Conformity error between "quantData" and "repInfo".')
    if(length(txFIdx) != nrow(quantData) ||
       length(hkIdx) != nrow(quantData)  ||
       length(unitFIdx) != ncol(quantData))
        stop('Conformity error between "quantData" and
              indices (txFIdx, hkIdx, unitFIdx).')
    if(!is.logical(txFIdx) || !is.logical(hkIdx) || !is.logical(unitFIdx))
        stop('Indices is not logical.')
    if(!is.numeric(unitCut) || length(unitCut) > 1)
        stop('"unitCut" is not single number.')
    # house keeping medians
    hkmed <- log2(apply(quantData[hkIdx, ], 2, median, na.rm = TRUE))
    # cut low signals with equivalent effects as selected units
    if(unitCut > 0){
        min2med <- median(hkmed[unitFIdx] - log2(unitCut), na.rm = TRUE)
        countCut <- 2^(hkmed - min2med)
        quantData[t(t(quantData) < countCut)] <- 0
    }
    # rolling back scaler (detrended 0 equals unit 1)
    scaler <- median(hkmed[unitFIdx])
    dat <- new("rnaseqcomp", quantData = quantData[txFIdx, ],
               repInfo = repInfo, hkmed = hkmed, scaler = scaler)
    return(dat)
}

