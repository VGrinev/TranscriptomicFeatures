#' Predict the true ORFs in circular RNA molecules
#' @description Predict the true ORFs in circular RNAs.
#' @param circRNAs character string giving the name of file with the sequences
#'     of circular RNAs. This file must be in the working directory. Allowed
#'     file formats are "fasta" or "fa".
#' @param orf_length_thr integer value giving the threshold for minimal length
#'     of ORF(-s). Default value is 42.
#' @param model character string giving the connection or full path to the file
#'     from which the classification model is read. Use default NULL value to 
#'     use default model from system extdata.
#' @param pr_thr probability threshold for the "winning" class of ORFs.
#'     The default value is NULL, which means that for a given circular RNA
#'     molecule, any ORF with the maximum probability value will be selected.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return The coordinates and sequences of ORFs in circular RNAs of interest.
#' @authors Vasily V. Grinev, Dai Xiaoxuan
#' @examples
#' ORFs <- predictORFsCircRNAs(circRNAs="circRNA_sequences_test.fa",
#'                             orf_length_thr=42,
#'                             model="classRFmodel#1.rds",
#'                             pr_thr=NULL,
#'                             workDir="D:/Vasily Grinev")
#' @export
#' Last updated: April 1, 2023.

predictORFsCircRNAs <- function(circRNAs,
                                orf_length_thr=42,
                                model=NULL,
                                pr_thr=NULL,
                                workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the base package parallel and
    #   packages Biostrings v.2.60.2, Rcpp v.1.0.8.3, randomForest v.4.6-14 and
    #   data.table v.1.14.2.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=parallel))
    suppressMessages(expr=library(package=Rcpp))
    suppressMessages(expr=library(package=randomForest))
    suppressMessages(expr=library(package=data.table))
    source(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/R/codonStartStop.r")
    source(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/R/findORFs.r")
    source(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/R/vectorizeORFs.r")
    sourceCpp(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/src/getCorrelationFactors.cpp")
    sourceCpp(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/src/getBaoMetrics.cpp")
    ### Loading of the circular RNA molecules as a list of character strings.
    ##  Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=circRNAs)
    ##  Validation of file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Full path to the file.
    circRNAs <- paste(workDir, circRNAs, sep="/")
    ##  Loading of the circular RNA molecules.
    circRNAseq <- as.list(x=as.character(x=readDNAStringSet(filepath=circRNAs)))
    ### Identification of all possible ORFs in circular RNAs of interest.
    cl <- makeCluster(spec=detectCores() - 1)
    clusterExport(cl=cl,
               varlist=c("findORFs", "DNAStringSet", "matchPDict", "DNAString",
                         "aggregate", "codonStartStop", "c", "outer", "round",
                         "length", "t", "lower.tri", "cbind", "list", "min",
                         "order", "substring", "do.call", "colnames", "return",
                         "names", "as.character", "sort", "unlist", "start"))
    all_orfs <- parLapply(X=circRNAseq,
                          fun=function(y){orf <- paste(rep(x=y, times=4),
                                                       collapse="")
                                          orf <- findORFs(x=orf)
                                          orf <- orf[as.numeric(x=orf[, 1]) <=
                                                     nchar(x=y) + 2, ]},
                          cl=cl)
    stopCluster(cl=cl)
    all_orfs <- cbind(transcript_id=rep(x=names(x=circRNAseq),
                                        times=lengths(all_orfs)/4),
                      do.call(what=rbind, args=all_orfs))
    all_orfs <- all_orfs[as.numeric(all_orfs[, "length"]) >= orf_length_thr, ]
    ### Classification of the all identified ORFs.
    ##  Annotation of the all identified ORFs with sequence features.
    prob_orfs <- vectorizeORFs(x=DNAStringSet(x=all_orfs[, "orf.sequence"]))
    ##  Import of the classification model.
    if(is.null(x=model)){
        model <- download_model_file()
    } else {
        model <- readRDS(file=file(description=paste(workDir, model, sep="/")))
    }
    ##  Classification.
    prob_orfs <- predict(object=model, newdata=prob_orfs, type="prob")
    prob_orfs <- data.table(cbind(all_orfs[, c("transcript_id",
                                               "start",
                                               "end",
                                               "length")],
                                  prob=prob_orfs[, "3"]),
                                  orf_seq=all_orfs[, c("orf.sequence")])
    ### Data aggregation and filtration.
    true_orfs <- prob_orfs[prob_orfs[, .I[prob == max(x=prob)],
                                     by=transcript_id]$V1]
    true_orfs <- true_orfs[true_orfs[, .I[length == max(x=length)],
                                     by=transcript_id]$V1]
    true_orfs <- true_orfs[true_orfs[, .I[start == min(x=start)],
                                     by=transcript_id]$V1]
    true_orfs <- data.frame(true_orfs)
    if (length(x=circRNAseq) > length(x=true_orfs$transcript_id)){
        mis <- names(x=circRNAseq)[!names(x=circRNAseq) %in% 
                                      true_orfs$transcript_id]
        mat <- cbind(matrix(0L, nrow=length(x=mis), ncol=4), "NA")
        colnames(x=mat) <- c("start", "end", "length", "prob", "orf_seq")
        true_orfs <- rbind(true_orfs, cbind(transcript_id=mis, mat))
    }
    true_orfs <- true_orfs[order(x=true_orfs$transcript_id), ]
    true_orfs[, c("start", "end", "length", "prob")] <-
                      apply(X=true_orfs[, c("start", "end", "length", "prob")],
                            MARGIN=2,
                            FUN=as.numeric)
    length_circRNA <- data.frame(cbind(names(x=circRNAseq),
                                       nchar(x=circRNAseq)))
    colnames(x=length_circRNA) <- c("transcript_id", "length")
    length_circRNA$length <- as.numeric(x=length_circRNA$length)
    length_circRNA <- length_circRNA[order(x=length_circRNA$transcript_id), ]
    rownames(x=length_circRNA) <- NULL
    true_orfs$orf_end <- length_circRNA$length - ((((true_orfs$end -
                         length_circRNA$length) %/%
                         length_circRNA$length) + 1) * length_circRNA$length -
                         (true_orfs$end - length_circRNA$length))
    true_orfs$orf_end <- paste(paste(((true_orfs$end -
                                     length_circRNA$length) %/%
                                     length_circRNA$length) + 1, "r", sep=""),
                               true_orfs$orf_end, sep="+")
    true_orfs$transcript_length <- length_circRNA$length
    true_orfs <- true_orfs[, c(1, 8, 2, 7, 4:6)]
    colnames(x=true_orfs) <- c("transcript_id", "transcript_length",
                               "orf_start", "orf_end", "orf_length",
                               "orf_probability", "orf_sequence")
    if(!is.null(x=pr_thr)){
        true_orfs <- true_orfs[true_orfs[, "prob"] >= pr_thr, ]
    }
    rownames(x=true_orfs) <- NULL
    ### Returning a final object of class data.frame.
    return(true_orfs)
}

download_model_file <- function(){
    path <- system.file("extdata", "cl_model.rds", package="ORFhunteR")
    model <- readRDS(path)
    if(length(model) < 1){
        cat("Loading model file (once)...")
        model <- readRDS(file=url("http://www.sstcenter.com/download/ORFhunteR/classRFmodel_1.rds", "r"))
        saveRDS(object=model, file=path)
        cat("Done.")
    }
    return(model)
}
