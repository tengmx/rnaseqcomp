#' @name encodeCells
#' @title Two Cell Lines from ENCODE by 9 Methods
#'
#' @description
#' This data includes RNA-seq quantifications for 56668 genes
#' in 2 cell lines with 2 replicates (GM12878 and K562), by
#' 9 methods including RSEM, Cufflinks, FluxCapacitor, eXpress
#' and Sailfish etc. Gene meta information is included such as
#' gene name, gene type and if house keeping genes. Also, fold
#' changes for part of protein coding genes from microarray is
#' provided between these cell lines.
#' 
#' @docType data
#' @format A series of objects including two 56668*18
#' quantification matrices (gm12878 & k562), one 56668*3
#' dataframe of gene meta information (genemeta), one factor
#' documenting pipelines of RNA-seq replicates (repInfo) and
#' fold change from microarray (arrayFC).
#'
#' @rdname encodeCells
#'
NULL
