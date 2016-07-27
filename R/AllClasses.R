#' @title Sanity Check of S4 rnaseqcomp Class
#'
#' @description This function always checks if the elements are valid
#' to create a S4 rnaseqcomp object. Specifically, check if \code{quantData}
#' is a list of matrices, if \code{condInfo} has the correct length and
#' levels, etc.
#'
#' @param object A object of S4 rnaseqcomp class
#'
#' @return TRUE, or character if error happens.
#'
check_rnaseqcomp <- function(object) {
    errors <- character()
    if(nlevels(object@condInfo) != 2){
        msg <- "Not two cell lines included."
        errors <- c(errors, msg)
    }
    negcells <- sum(sapply(object@quantData, function(x)
                           sum(x<0, na.rm = TRUE)))
    if(negcells > 0){
        msg <- '"quantData" have negative matrix.'
        errors <- c(errors, msg)
    }
    if(ncol(object@quantData[[1]]) != length(object@repInfo)){
        msg <- '"quantData" column size must equal to "refInfo" length.'
        errors <- c(errors, msg)
    }
    if(ncol(object@quantData[[1]]) != length(object@condInfo)){
        msg <- '"quantData" column size must equal to "condInfo" length.'
        errors <- c(errors, msg)
    }
    if(length(object@quantData)<1){
        msg <- '"quantData" must have at least one matrix.'
        errors <- c(errors, msg)
    }
    if(length(object@refMed) != length(object@quantData)){
        msg <- 'Length must be equivalent for "refMed" and "quantData".'
        errors <- c(errors, msg)
    }
    if(length(object@scaler) > 1){
        msg <- '"scaler" must be one single number.'
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}

#' @title rnaseqcomp
#'
#' @description
#' This is a S4 class to organize data ready for benchmark summarization.
#' There are 5 S3 objects inside this class. \code{quantData} documents a
#' list of data matrices ready for evaluation by functions \code{plotSD},
#' \code{plotNE}, \code{plot2TX} or \code{plotROC}. \code{condInfo} is a
#' factor corresponding to columns of \code{quantData} matrices, indicating
#' to which cell lines each sample belongs. \code{repInfo} is a factor
#' corresponding to columns of \code{quantData} matrices indicating replicate
#' information. \code{repInfo} is a legacy from previous versions, and
#' doesn't have too much meanings in current version. \code{refMed} is the
#' median log2 signal of calibration references. \code{scaler} is a number
#' that point to the median log2 signal of reference pipeline. \code{refMed}
#' and \code{scaler} were used to calibrate and generate \code{quantData}.
#'
#' @exportClass rnaseqcomp
#'
setClass(Class = "rnaseqcomp",
         representation = representation(
         quantData = "list", condInfo = "factor",
         repInfo = "factor", refMed = "list", scaler = "numeric"),
         validity = check_rnaseqcomp)


setMethod("show", "rnaseqcomp", function(object){
    cat("rnaseqcomp: Benchmarks for RNA-seq quantification pipelines\n\n")
    cat("Quantifications pipelins: ", length(object@quantData), "\n")
    cat("Total transcripts: ", nrow(object@quantData[[1]]), "\n")
    cat("Total samples from 2 conditions: ", ncol(object@quantData[[1]]), "\n")
})

