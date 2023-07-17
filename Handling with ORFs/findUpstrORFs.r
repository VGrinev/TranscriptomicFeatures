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

#' @authors Vasily V. Grinev
#' @examples
#' tr <- "Kasumi-1, transcriptome, StringTie, filtered, alternative assemblies.fasta"
#' orfs <- "Kasumi-1, transcriptome, StringTie, filtered, alternative assemblies, ORFs.txt"
#' res <- findUpstrORFs(tr=tr,
#'                      orfs=orfs,
#'                      prThr=0.9,
#'                      workDir="D:/Vasily Grinev")
#' @export
#' Last updated: June 21, 2023.

findUpstrORFs <- function(tr,
                          orfs,
                          prThr=0.9,
                          workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with packages Biostrings v.2.60.2 and
    #   GenomicRanges v.1.40.0.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=GenomicRanges))
    ### Loading of the RNA molecules as an object of class DNAStringSet.
    ##  Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=tr)
    ##  Validation of file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Full path to the file.
    trans <- paste(workDir, tr, sep="/")
    ##  Loading of the RNA molecules.
    trans <- readDNAStringSet(filepath=trans)
    trans <- trans[order(x=names(x=trans)), ]
    ### Loading the coordinates of identified open reading frames as an object
    #   of class data frame.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=orfs)
    ##  Validation of file format.
    if (!frt %in% "txt"){
        stop("Invalid file format")
    }
    ##  Full path to the file.
    ORFs <- paste(workDir, orfs, sep="/")
    ##  Loading of the coordinates.
    coordORFs <- read.table(file=ORFs, header=TRUE, quote="\"", as.is=TRUE)
    coordORFs <- coordORFs[coordORFs$orf_prob >= prThr, ]
    coordORFs <- coordORFs[order(x=coordORFs$transcript_id), ]

    tr_start <- unlist(x=vmatchPattern(pattern="ATG", subject=trans))
    tr_start <- data.frame(tr_start)
    tr_start <- tr_start[, c(4, 1:2)]
    colnames(x=tr_start) <- c("seqnames", "start", "end")
    tr_stop <- unlist(x=c(vmatchPattern(pattern="TAA", subject=trans),
                          vmatchPattern(pattern="TAG", subject=trans),
                          vmatchPattern(pattern="TGA", subject=trans)))
    tr_stop <- data.frame(tr_stop)
    tr_stop <- tr_stop[, c(4, 1:2)]
    colnames(x=tr_stop) <- c("seqnames", "start", "end")
    tr_start <- tr_start[tr_start$seqnames %in% tr_stop$seqnames, ]
    tr_start <- tr_start[order(x=tr_start$seqnames), ]
    rownames(x=tr_start) <- NULL
    tr_stop <- tr_stop[tr_stop$seqnames %in% tr_start$seqnames, ]
    tr_stop <- tr_stop[order(x=tr_stop$seqnames), ]
    rownames(x=tr_stop) <- NULL
    tr_start <- split(x=tr_start, f=tr_start$seqnames)
    tr_stop <- split(x=tr_stop, f=tr_stop$seqnames)
    ### Equalness of loaded data.
    trans <- trans[names(x=trans) %in% coordORFs$transcript_id, ]
    trans <- trans[names(x=trans) %in% names(x=tr_start)]
    tr_start <- tr_start[names(x=tr_start) %in% names(x=trans)]
    tr_stop <- tr_stop[names(x=tr_stop) %in% names(x=trans)]
    ### Identification of upstream open reading frames.
    uORFs <- list()
    for (i in 1:length(x=tr_start)){
        inFrame <- outer(X=tr_stop[[i]]$start,
                         Y=tr_start[[i]]$start,
                         FUN="-")/3
        inFrame <- which(x=round(x=inFrame) == inFrame & inFrame > 0,
                         arr.ind=TRUE)
        if (length(x=inFrame) > 0){
            inFrame[, 1] <- tr_stop[[i]]$start[inFrame[, 1]]
            inFrame[, 2] <- tr_start[[i]]$start[inFrame[, 2]]
            inFrame <- aggregate(x=inFrame[, 1],
                                 by=list(inFrame[, 2]),
                                 FUN=min)
            inFrame <- aggregate(x=inFrame[, 1],
                                 by=list(inFrame[, 2]),
                                 FUN=min)
            inFrame <- inFrame[, c(2, 1)]
            inFrame <- inFrame[order(inFrame[, 1]), ]
            inFrame <- inFrame[inFrame[, 2] - inFrame[, 1] > 24, ]
            if (nrow(x=inFrame) > 0){
                inFrame[, 2] <- inFrame[, 2] + 2
                inFrame <- cbind(names(x=tr_stop[i]),
                                 inFrame,
                                 substring(text=as.character(x=trans[[i]]),
                                           first=inFrame[, "x"],
                                           last=inFrame[, "Group.1"]))
                colnames(x=inFrame) <- c("seqnames",
                                         "start", "end",
                                         "orf.sequence")
                rownames(x=inFrame) <- NULL
                uORFs[[i]] <- inFrame
            }
        }
    }
    uORFs <- do.call(what=rbind, args=uORFs)
    uORFs <- makeGRangesFromDataFrame(df=uORFs, keep.extra.columns=TRUE)
    trans <- trans[names(x=trans) %in% seqnames(x=uORFs)]
    tr_start <- tr_start[names(x=tr_start) %in% names(x=trans)]
    tr_stop <- tr_stop[names(x=tr_stop) %in% names(x=trans)]
    coordORFs <- coordORFs[coordORFs$transcript_id %in% names(x=trans), ]
    coordORFs <- coordORFs[, c(1, 4:5)]
    colnames(x=coordORFs) <- c("seqnames", "start", "end")
    rownames(x=coordORFs) <- NULL
    coordORFs <- makeGRangesFromDataFrame(df=coordORFs)
    hits <- findOverlaps(query=uORFs, subject=coordORFs, type="equal")
    uORFs <- uORFs[-queryHits(x=hits), ]
    end(x=coordORFs) <- start(x=coordORFs)
    start(x=coordORFs) <- 1
    uORFs <- subsetByOverlaps(x=uORFs, ranges=coordORFs, type="any")
    uORFs <- data.frame(uORFs)
    uORFs$seqnames <- as.character(x=uORFs$seqnames)
    uORFs$strand <- as.character(x=uORFs$strand)
    uORFs_id <- unlist(x=sapply(X=table(uORFs$seqnames),
                                FUN=function(y){paste(".uORF",
                                                      formatC(x=seq(from=1,
                                                                    to=y),
                                                              width=3,
                                                              flag="0"),
                                                      sep="")}))
    uORFs$uORFs_id <- paste(uORFs$seqnames, as.vector(x=uORFs_id), sep="")
    uORFs <- uORFs[, c(1, 7, 2:4, 6)]
    colnames(uORFs) <- c("sequence_id", "orf_id",
                         "orf_start", "orf_end", "orf_length",
                         "orf_sequence")
    uORFs <- list(uORFs=uORFs, seqs=DNAStringSet(x=uORFs[, 6]))
    names(x=uORFs$seqs) <- uORFs$uORFs$orf_id
    ### Returning the final object.
    return(uORFs)
}

tr <- "Kasumi-1, transcriptome, StringTie, filtered, all transcripts.fasta"
orfs <- "Kasumi-1, transcriptome, StringTie, filtered, all transcripts, ORFs.txt"
res <- findUpstrORFs(tr=tr,
                     orfs=orfs,
                     prThr=0.9,
                     workDir="D:/Vasily Grinev")
setwd(dir="D:/Vasily Grinev")
write.table(x=res$uORFs,
            file="Kasumi-1, transcriptome, StringTie, filtered, all transcripts, uORFs.txt",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)
prob_orfs <- vectorizeORFs(x=res$seqs)
model <- readRDS(file=file(description="D:/Vasily Grinev/classRFmodel_1.rds"))
prob_orfs <- predict(object=model, newdata=prob_orfs, type="prob")
prob_orfs <- data.table(cbind(res$uORFs, prob=prob_orfs[, "3"]))
true_orfs <- prob_orfs[prob_orfs[, .I[prob == max(x=prob)],
                                 by=sequence_id]$V1]
true_orfs <- true_orfs[true_orfs[, .I[orf_length == max(x=orf_length)],
                                 by=sequence_id]$V1]
true_orfs <- true_orfs[true_orfs[, .I[orf_start == min(x=orf_start)],
                                 by=sequence_id]$V1]
true_orfs <- data.frame(true_orfs)
true_orfs <- true_orfs[order(x=true_orfs$sequence_id), ]
true_orfs1 <- prob_orfs[, c(1:2, 7, 3:6)]
colnames(x=true_orfs1) <- c("transcript_id",
                           "orf_id",
                           "orf_prob",
                           "orf_start", "orf_end", "orf_length",
                           "orf_sequence")
rownames(x=true_orfs1) <- NULL
write.table(x=true_orfs1,
            file="uorfs.txt",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)
