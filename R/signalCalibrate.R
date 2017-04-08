#' @title Quantification Filtering And Calibration
#'
#' @description This is the function to do any pre-filtering or
#' pre-processing analysis for downstream benchmark estimation
#' and visualization. Pre-filtering includes row selection (e.g.
#' protein coding genes) of quantification table; pre-processing
#' includes calculation on a set of rows as calibration reference
#' (e.g. house keeping genes) across different quantification
#' pipelines, calibration of quantifications across all the
#' pipelines based on given cutoffs from selected pipelines.
#'
#' @param quantData A list of quantification matrices each with rows by
#' features (transcripts, genes, junctions or exons) and columns
#' by samples. Names of the list should be provided. The sizes of each
#' element should be the same. Missing data will be set to 0.
#'
#' @param condInfo A factor documenting condition information of samples,
#' corresponding to the columns of each element in \code{quantData}.
#'
#' @param repInfo A factor documenting replicate information of samples,
#' corresponding to the columns of each element in \code{quantData}.
#'
#' @param evaluationFeature A logical vector corresponding to the
#' rows of each element in \code{quantData}, providing which features
#' should be considered for downstream evaluation, e.g. protein coding
#' genes.
#'
#' @param calibrationFeature A logical vector corresponding to the
#' rows of each element in \code{quantData}, providing which features
#' should be
#' considered as calibration reference, e.g. house keeping genes.
#'
#' @param unitReference A numeric number specifying which pipeline will
#' be selected as reference pipeline, i.e. the index of
#' one element in \code{quantData}.
#'
#' @param unitCutoff A numeric number for signal cutoff on reference
#' pipeline specified by \code{unitReference} (default: 0).
#' Equivalent effects of cutoffs will be applied to other pipelines
#' accordingly.
#'
#' @param calibrationFeature2 A logical vector corresponding to the
#' rows of each element in \code{quantData}, providing which features
#' should be considered as references for calibration across
#' different datasets. Default \code{NULL} means no calibration needed.
#'
#' @param fixMedian A numeric number specifying the median of detrend
#' logsignals for features specified by \code{calibrationFeature2}.
#' When comparing across datasets, those features will be calibrated to
#' have the same median as \code{fixMedian}, while other features
#' calibrated accordingly. The default is 4.776, which was calculated
#' based on one ENCODE dataset used in our web tool.
#'
#' @return A \code{rnaseqcomp} S4 class object
#' \item{quantData}{A filtered and calibrated list of quantifications
#' for downstream analysis.}
#' \item{condInfo}{A factor documenting sample condition information.}
#' \item{repInfo}{A factor documenting sample replicate information.}
#' \item{refMed}{A list of numeric vectors giving the log scale medians
#' of calibration features in different pipelines.}
#' \item{scaler}{A number that was used for scaling quantifications onto
#' reference pipeline.}
#'
#' @details
#' In the functions \code{plotSD} and \code{plot2TX}, detrended
#' signals with value 0 will be at the same level as value 1 for
#' giving pipeline by \code{unitReference}.
#'
#' @importFrom methods show
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

signalCalibrate <- function(quantData, condInfo, repInfo, evaluationFeature,
                         calibrationFeature, unitReference, unitCutoff = 0,
                         calibrationFeature2 = NULL, fixMedian = 4.776){
    if(!is.list(quantData) ||
       length(unique(sapply(quantData,class)))>1 ||
       unique(sapply(quantData,class)) != "matrix" ||
       sum(sapply(quantData,dim) - dim(quantData[[1]]) != 0) > 0)
        stop('"quantData" is not a list of same size matrices.')
    if(sum(sapply(quantData, function(x) sum(x<0, na.rm = TRUE))) > 0)
        stop('matrices in "quantData" have negative values.')
    if(is.null(names(quantData)))
        stop('"quantData" shoud have list names.')
    if(!is.factor(condInfo) || nlevels(condInfo) != 2 ||
       length(condInfo) != ncol(quantData[[1]]))
        stop(paste0('"condInfo" is not factor or not two conditions or',
                    ' length different from columns of "quantData" matrices.'))
    if(!is.factor(repInfo) || length(repInfo) != ncol(quantData[[1]]))
        stop(paste0('"refInfo" is not factor or its length not equal to the ',
                    'column number of "quantData" matrices.'))
    if(length(unitReference)!=1 || !(unitReference %in% seq_along(quantData)))
        stop(paste0('"unitReference" must be one number indicating the index',
                    'for one of "quantData" matrices.'))
    if(length(evaluationFeature) != nrow(quantData[[1]]))
        stop(paste0('The length of "evaluationFeature" must be the same as',
                    'the row number of "quantData" matrices.'))
    if(length(calibrationFeature) != nrow(quantData[[1]]))
        stop(paste0('The length of "calibrationFeature" must be the same as',
                    'the row number of "quantData" matrices.'))
    if(!is.null(calibrationFeature2) &&
       length(calibrationFeature2) != nrow(quantData[[1]]))
        stop(paste0('The length of "calibrationFeature2" must be the same as',
                    'the row number of "quantData" matrices.'))
    if(!is.logical(evaluationFeature))
        stop('"evaluationFeature" must be a logical vector.')
    if(!is.logical(calibrationFeature))
        stop('"calibrationFeature" must be a logical vector.')
    if(!is.null(calibrationFeature2) && !is.logical(calibrationFeature2))
        stop('"calibrationFeature2" must be a logical vector.')
    if(!is.numeric(unitCutoff) || length(unitCutoff) > 1 ||
       unitCutoff < 0)
        stop('"unitCutoff" must be one non-negative number.')
    if(!is.numeric(fixMedian) || length(fixMedian) > 1)
        stop('"fixMedian" must be a single number.')
    quantData <- lapply(quantData, function(x){
        x[is.na(x)] <- 0
        x
    })
    # house keeping medians
    refMed <- lapply(quantData,function(x)
              log2(apply(x[calibrationFeature, ], 2, median, na.rm = TRUE)))
    # rolling back scaler (detrended 0 equals unit 1)
    scaler <- median(unlist(refMed[unitReference]))
    if(!is.null(calibrationFeature2)){
        refMed2 <- lapply(quantData,function(x)
              log2(apply(x[calibrationFeature2, ], 2, median, na.rm = TRUE)))
        scaler2 <- median(unlist(refMed2[unitReference]))
        scaler <- scaler + fixMedian - scaler2
    }
    # calibration
    quantDataC <- lapply(seq_along(quantData), function(i){
        x <- t(t(quantData[[i]]) * 2^(scaler-refMed[[i]]))
        # cut low signals
        if(unitCutoff > 0)
            x[x < unitCutoff] <- 0
        x[evaluationFeature,]
    })
    names(quantDataC) <- names(quantData)
    dat <- new("rnaseqcomp", quantData = quantDataC, condInfo = condInfo,
               repInfo = repInfo, refMed = refMed, scaler = scaler)
    return(dat)
}
