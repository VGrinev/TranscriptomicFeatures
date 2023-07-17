#' Predict true open reading frames in RNA molecules
#' @description Prediction of the true open reading frames in RNA molecules.
#' @param tr character string giving the name of file with experimental
#'     transcripts. Allowed file formats are "fasta", "fa", "gtf" or "gff".
#'     Alternatively, it may be a ready-to-use output object of function
#'     findAllORFs() stored in the computer's RAM.
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param codStart character string with type of start codon: "ATG", "GTG",
#'     "TTG" or "CTG". Default value is "ATG".
#' @param prThr probability threshold for the "winning" class of open reading
#'     frames. Default value is NULL.
#' @param model character string giving the connection or full path to the file
#'     from which the classification model is read. Use default NULL value to 
#'     use our default model.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return An object of class data frame containing the following fields:
#'     i) transcript_id  - ID of sequence;
#'     ii) orf_id        - ID of open reading frame;
#'     iii) orf_prob     - probability that given open reading frame is true
#'                         open reading frame;
#'     iv) orf_start     - start coordinate of open reading frame in a
#'                         sequence;
#'     v) orf_end        - end coordinate of open reading frame in a sequence;
#'     vi) prf_length    - length of open reading frame;
#'     vii) orf_sequence - sequence of open reading frame.
#' @author Mikalai M. Yatskou, Vasily V. Grinev.
#' @examples
#' tr <- "Kasumi-1, transcriptome, StringTie, filtered, all transcripts, part 1.fasta"
#' genome <- "BSgenome.Hsapiens.UCSC.hg38"
#' model <- "D:/Vasily Grinev/classRFmodel_1.rds"
#' workDir <- "D:/Vasily Grinev"
#' ORFs <- predictORF(tr=tr,
#'                    genome=genome,
#'                    codStart="ATG",
#'                    model=model,
#'                    prThr=NULL,
#'                    workDir=workDir)
#' @export
#' Last updated: April 19, 2023.

predictORF <- function(tr,
                       genome="BSgenome.Hsapiens.UCSC.hg38",
                       codStart="ATG",
                       model=NULL,
                       prThr=NULL,
                       workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the packages Rcpp v.1.0.8.3,
    #   randomForest v.4.6-14 and data.table v.1.14.2.
    suppressMessages(expr=library(package=Rcpp))
    suppressMessages(expr=library(package=randomForest))
    suppressMessages(expr=library(package=data.table))
    source(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/R/loadTrExper.r")
    source(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/R/findAllORFs.r")
    source(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/R/vectorizeORFs.r")
    sourceCpp(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/src/getCorrelationFactors.cpp")
    sourceCpp(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/src/getBaoMetrics.cpp")
    ### Identification of all possible open reading frames in RNA molecules.
    if (class(x=tr) == "character"){
        ##  Loading of the experimental transcripts as an object of class
        #   DNAStringSet.
        exp_trans <- loadTrExper(tr=tr, genome=genome, workDir=workDir)
        ##  Identification of all possible open reading frames.
        all_orfs <- findAllORFs(tr=exp_trans)
    }else{
        all_orfs <- tr
    }
    ### Classification of the all identified open reading frames.
    ##  Annotation of the all identified open reading frames with features.
    prob_orfs <- vectorizeORFs(x=all_orfs$seqs)
    ##  Import of the classification model.
    if(is.null(x=model)){
        model <- download_model_file()
    }else{
        model <- readRDS(file=file(description=model))
    }
    ##  Classification.
    prob_orfs <- predict(object=model, newdata=prob_orfs, type="prob")
    prob_orfs <- data.table(cbind(all_orfs$ORFs, prob=prob_orfs[, "3"]))
    ### Data aggregation and filtration.
    true_orfs <- prob_orfs[prob_orfs[, .I[prob == max(x=prob)],
                                     by=sequence_id]$V1]
    true_orfs <- true_orfs[true_orfs[, .I[orf_length == max(x=orf_length)],
                                     by=sequence_id]$V1]
    true_orfs <- true_orfs[true_orfs[, .I[orf_start == min(x=orf_start)],
                                     by=sequence_id]$V1]
    true_orfs <- data.frame(true_orfs)
    if (length(x=exp_trans) > length(x=true_orfs$sequence_id)){
        mis <- names(x=exp_trans)[!names(x=exp_trans) %in% 
                                      true_orfs$sequence_id]
        mat <- matrix(NA, nrow=length(x=mis), ncol=6)
        colnames(x=mat) <- c("orf_id",
                             "orf_start", "orf_end", "orf_length",
                             "orf_sequence",
                             "prob")
        true_orfs <- rbind(true_orfs, cbind(sequence_id=mis, mat))
    }
    true_orfs <- true_orfs[order(x=true_orfs$sequence_id), ]
    true_orfs <- true_orfs[, c(1:2, 7, 3:6)]
    colnames(x=true_orfs) <- c("transcript_id",
                               "orf_id",
                               "orf_prob",
                               "orf_start", "orf_end", "orf_length",
                               "orf_sequence")
    if(!is.null(x=prThr)){
        true_orfs <- true_orfs[true_orfs[, "orf_prob"] >= prThr, ]
    }
    rownames(x=true_orfs) <- NULL
    ### Returning a final object of class data frame.
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
