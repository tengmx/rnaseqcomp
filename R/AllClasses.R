#' @title rnaseqcomp
#'
#' @exportClass rnaseqcomp
#'
check_rnaseqcomp <- function(object) {
    errors <- character()
    repsnot2 <- sum(summary(object@repInfo) != 2)
    if(repsnot2 > 0 ){
        msg <- "Replicates of each pipeline should be 2."
        errors <- c(errors, msg)
    }
    negcells <- sum(object@quantData < 0, na.rm = TRUE)
    if(negcells > 0){
        msg <- "'quantData' should be a non-negative matrix."
        errors <- c(errors, msg)
    }
    if(ncol(object@quantData) != length(object@repInfo)){
        msg <- "'quantData' column size should equal to 'refInfo' length."
        errors <- c(errors, msg)
    }
    if(length(object@hkmed) != length(object@repInfo)){
        msg <- "Size should be equivalent for 'refInfo' and 'hkmed'."
        errors <- c(errors, msg)
    }
    if(length(object@scaler) > 1){
        msg <- "'scaler' should be one single number."
        errors <- c(errors, msg)
    }
    if (length(errors) == 0) TRUE else errors
}

setClass(Class = "rnaseqcomp",
         representation = representation(
             quantData = "matrix", repInfo = "factor",
             hkmed = "numeric", scaler = "numeric"),
         validity = check_rnaseqcomp)


setMethod("show", "rnaseqcomp", function(object){
    cat("rnaseqcomp: Benchmark for RNA-seq quantification pipelines\n\n")
    cat("Reps:\n", as.character(object@repInfo), "\n\n")
    cat("House-Keeping log2Median:\n", object@hkmed, "\n\n")
    cat("Detrened signal scaler:\n", object@scaler, "\n\n")
    rown <- nrow(object@quantData)
    coln <- ncol(object@quantData)
    cat("Quantification data has ", rown,
        " rows and ", coln, " columns:\n")
    if(coln > 8 & rown > 10){
        a <- object@quantData[1:4, 1:4]
        b <- object@quantData[1:4, (coln-3):coln]
        c <- object@quantData[(rown-3):rown, 1:4]
        d <- object@quantData[(rown-3):rown, (coln-3):coln]
        text <- rbind(data.frame(a, "...", b), "...",
                      data.frame(c, "...", d))
        colnames(text) <- c(as.character(object@repInfo[1:4]), ".",
                            as.character(object@repInfo[(coln-3):coln]))
        rownames(text)[5] <- "."
    }else if (rown > 10){
        a <- data.frame(object@quantData[1:4, ])
        b <- data.frame(object@quantData[(rown-3):rown, ])
        text <- rbind(a, "...", b)
        colnames(text) <- as.character(object@repInfo)
        rownames(text)[5] <- "."
    }else if (coln > 8){
        a <- data.frame(object@quantData[, 1:4])
        b <- data.frame(object@quantData[, (coln-3):coln])
        text <- cbind(a, "...", b)
        colnames(text) <- c(as.character(object@repInfo[1:4]), ".",
                            as.character(object@repInfo[(coln-3):coln]))
    }else{
        text <- data.frame(object@quantData)
        colnames(text) <- as.character(object@repInfo)
    }
    print(text)
})
          
