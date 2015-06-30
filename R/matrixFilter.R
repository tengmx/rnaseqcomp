#' @title Quantification Table Filtering
#'
#' @description This is the function to do any pre-filtering or
#' pre-processing analysis for downstream benchmark estimation
#' and visualization. Pre-filtering includes row selection (e.g.
#' protein coding genes) of quantification table ; pre-processing
#' includes calculation on a set of rows as calibration reference
#' (e.g. house keeping genes) across different quantification
#' pipelines, estimation of confidence thresholds on all the
#' pipelines based on given thresholds from selected pipelines.
#'
#' @param quantData A matrix of quantifications with rows by
#' features (genes, transcripts, junctions or exons) and columns
#' by pipelines each with 2 replicates. Missing data is allowed
#' as NA.
#' @param repInfo A factor documenting quantification pipeline
#' names corresponding to the columns of \code{quantData}.
#' @param evaluationFeature A logical vector corresponding to the
#' rows of \code{quantData}, providing which features should be
#' considered for downstream evaluation, e.g. protein coding genes.
#' @param calibrationFeature A logical vector corresponding to the
#' rows of \code{quantData}, providing which features should be
#' considered as calibration reference, e.g. house keeping genes.
#' @param unitReference A logical vector corresponding to the columns
#' of \code{quantData}, providing to which columns all units should be
#' unified and given threshold \code{unitCutoff} should be applied.
#' If multiple columns are chosen, it is important that these columns
#' have the same units (e.g. FPKM).
#' @param unitCutoff A numeric threshold for signal cutoff at given
#' quantification table columns by \code{unitReference} (default: 1).
#' Same effects cutoff are estimated for other columns accordingly.
#'
#' @return A \code{rnaseqcomp} S4 class object
#' \item{quantData}{A filtered table for downstream analysis.}
#' \item{repInfo}{A factor documenting quantification pipeline
#' names corresponding to the columns of \code{quantData}.}
#' \item{refMed}{A numeric vector giving the log scale medians
#' of calibration reference.}
#' \item{scaler}{A number providing the scales for downstream
#' detrended signal.}
#'
#' @details
#' In the functions \code{plotMAD} and \code{plotNE}, detrended
#' signals with value 0 will be at the same level as value 1 for
#' giving pipelines by \code{unitReference}.
#'
#' @importFrom methods show
#' @export
#'
#' @examples
#' data(encodeCells)
#' evaluationFeature <- encodeCells$genemeta$type == "protein_coding"
#' calibrationFeature <- encodeCells$genemeta$housekeeping
#' unitReference <- grepl("Cufflinks",encodeCells$repInfo)
#' dat <- matrixFilter(encodeCells$gm12878,encodeCells$repInfo,
#' evaluationFeature,calibrationFeature,unitReference)

matrixFilter <- function(quantData, repInfo, evaluationFeature,
                         calibrationFeature, unitReference, unitCutoff = 1){
    if(!is.factor(repInfo) || sum(summary(repInfo) != 2) > 0 )
        stop('"repInfo" is not factor or not fitting two replicates.')
    if(!is.matrix(quantData) || sum(quantData < 0, na.rm = TRUE) > 0)
        stop('"quantData" is not a non-negative matrix.')
    if(length(repInfo) != ncol(quantData))
        stop(paste0('The length of "repInfo" must be the same as the number',
              'of columns of "quantData" matrix.'))
    if(length(unitReference) != ncol(quantData))
        stop(paste0('The length of "unitReference" must be the same as the',
                    'number of columns of "quantData" matrix.'))
    if(length(evaluationFeature) != nrow(quantData))
        stop(paste0('The length of "evaluationFeature" must be the same as',
                    'the number of rows of "quantData" matrix.'))
    if(length(calibrationFeature) != nrow(quantData))
        stop(paste0('The length of "calibrationFeature" must be the same as',
                    'the number of rows of "quantData" matrix.'))
    if(!is.logical(evaluationFeature))
        stop('"evaluationFeature" must be a logical vector.')
    if(!is.logical(calibrationFeature))
        stop('"calibrationFeature" must be a logical vector.')
    if(!is.logical(unitReference))
        stop('"unitReference" must be a logical vector.')
    if(!is.numeric(unitCutoff) || length(unitCutoff) > 1)
        stop('"unitCutoff" must be a single number.')
    # house keeping medians
    refMed <- log2(apply(quantData[calibrationFeature, ], 2,
                        median, na.rm = TRUE))
    # cut low signals with equivalent effects as selected units
    if(unitCutoff > 0){
        min2med <- median(refMed[unitReference] - log2(unitCutoff),
                          na.rm = TRUE)
        countCut <- 2^(refMed - min2med)
        quantData[t(t(quantData) < countCut)] <- 0
    }
    # rolling back scaler (detrended 0 equals unit 1)
    scaler <- median(refMed[unitReference])
    dat <- new("rnaseqcomp", quantData = quantData[evaluationFeature, ],
               repInfo = repInfo, refMed = refMed, scaler = scaler)
    return(dat)
}

