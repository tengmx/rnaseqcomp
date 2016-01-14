#' @name simdata
#' @title Example of Quantifications on Simulation Data
#'
#' @description
#' This dataset include quantifications on 15776 transcripts on
#' two cell lines each with 8 replicates. The true differential
#' expressed transcripts were simulated. Quantifications from
#' two pipelines (RSEM and FluxCapacitor) are included in this
#' dataset at \code{simdata$quant}. Meta information of transcripts
#' is included at \code{simdata$meta}, inlcuding if they belongs
#' to house keeping genes and their true fold change status.
#' Sample information is included at \code{simdata$samp}.
#'
#' @docType data
#' @format A list of objects including list of two 15776*16
#' quantification matrices, one 15776*6 data frame with meta
#' information and one 16*3 data frame with sample information.
#'
#' @rdname simdata
#'
NULL
