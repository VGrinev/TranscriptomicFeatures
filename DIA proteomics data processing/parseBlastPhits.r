#' Parse of blastp hits table
#' @description Parsing of the blastp hits table.
#' @param x name of input file with raw blastp hits.
#' @param pr name of FASTA file with protein sequences of interest.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return an object of class data frame.
#' @author Vasily V. Grinev
#' @examples
#' res <- parseBlastPhits(x="Kasumi-1 proteome, all proteins detected at FDR 5%, validated fragments, all blastp hits.txt",
#'                        pr="Kasumi-1 proteome, all proteins detected at FDR 5%, validated fragments.fasta",
#'                        workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 10, 2023.

parseBlastPhits <- function(x, pr, workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with package Biostrings v.2.56.0.
    suppressMessages(expr=library(package=Biostrings))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of the blastp output as data frame object.
    blastp <- read.table(file=paste(workDir, x, sep = "/"),
                         sep="\t",
                         header=FALSE,
                         quote="\"",
                         as.is=TRUE)
    blastp$length <- 0
    blastp <- blastp[, c(1, 14, 2, 4:11, 3)]
    colnames(x=blastp) <- c("query",
                            "length",
                            "subject",
                            "alignment_length", "mismatches", "gap_opens",
                            "q.start", "q.end",
                            "s.start", "s.end",
                            "e_value",
                            "identity1")
    ### Loading of the in silico predicted proteins as AAStringSet object.
    proteins <- readAAStringSet(filepath=paste(workDir, pr, sep = "/"))
    ### Calculation of overall identity.
    blastp$identity2 <- 0
    blastp$length <- width(x=proteins)[match(blastp$query, names(x=proteins))]
    blastp$identity2 <- (blastp$q.end - blastp$q.start + 1)/blastp$length
    blastp$identity2 <- blastp$identity1 * blastp$identity2
    ### Parsing of blastp hits.
    blastp <- split(x=blastp, f=blastp$query)
    full_set <- list()
    for (i in 1:length(x=blastp)){
        obj <- blastp[[i]]
        obj <- obj[obj$identity2 == max(obj$identity2), ]
        obj <- obj[!duplicated(x=obj), ]
        full_set[[i]] <- obj[1, ]
    }
    full_set <- do.call(what=rbind, args=full_set)
    full_set$nonidentical_length <- round(x=full_set$length - full_set$length *
                                          full_set$identity2/100)
    full_set$N <- full_set$q.start - 1
    full_set$C <- full_set$length - full_set$q.end
    full_set$class <- "I"
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches == 0 & 
                      full_set$gap_opens == 0 &
                      full_set$C == 0, ]$class <- "N", silent=TRUE)
    try(expr=full_set[full_set$N == 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens == 0 &
                      full_set$C == 0, ]$class <- "M", silent=TRUE)
    try(expr=full_set[full_set$N == 0 &
                      full_set$mismatches == 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C == 0, ]$class <- "G", silent=TRUE)
    try(expr=full_set[full_set$N == 0 &
                      full_set$mismatches == 0 & 
                      full_set$gap_opens == 0 &
                      full_set$C > 0, ]$class <- "C", silent=TRUE)
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens == 0 &
                      full_set$C == 0, ]$class <- "N+M", silent=TRUE)
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches == 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C == 0, ]$class <- "N+G", silent=TRUE)
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches == 0 & 
                      full_set$gap_opens == 0 &
                      full_set$C > 0, ]$class <- "N+C", silent=TRUE)
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C == 0, ]$class <- "N+M+G", silent=TRUE)
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens == 0 &
                      full_set$C > 0, ]$class <- "N+M+C", silent=TRUE)
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches == 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C > 0, ]$class <- "N+G+C", silent=TRUE)
    try(expr=full_set[full_set$N > 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C > 0, ]$class <- "N+M+G+C", silent=TRUE)
    try(expr=full_set[full_set$N == 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C == 0, ]$class <- "M+G", silent=TRUE)
    try(expr=full_set[full_set$N == 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens == 0 &
                      full_set$C > 0, ]$class <- "M+C", silent=TRUE)
    try(expr=full_set[full_set$N == 0 &
                      full_set$mismatches > 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C > 0, ]$class <- "M+G+C", silent=TRUE)
    try(expr=full_set[full_set$N == 0 &
                      full_set$mismatches == 0 & 
                      full_set$gap_opens > 0 &
                      full_set$C > 0, ]$class <- "G+C", silent=TRUE)
    ### Returning of the final results.
    return(full_set)
}

setwd(dir="D:/Vasily Grinev")
res <- parseBlastPhits(x="Kasumi-1 proteome, all proteins detected at FDR 5%, validated fragments, all blastp hits.txt",
                       pr="Kasumi-1 proteome, all proteins detected at FDR 5%, validated fragments.fasta",
                       workDir="D:/Vasily Grinev")
write.table(x=res,
            file="Kasumi-1 proteome, all proteins detected at FDR 5%, validated fragments, the best blastp hits.txt",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)
