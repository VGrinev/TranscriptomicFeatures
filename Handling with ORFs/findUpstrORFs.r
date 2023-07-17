#' Find all upstream open reading frames in mRNA molecule(-s)
#' @description Finding all variants of the upstream open reading frames in
#'     mRNA molecule(-s).
#' @param tr character string giving the name of file with the sequence(-s) of
#'     the mRNA molecule(-s). This file must be in the working directory.
#'     Allowed file formats are "fasta" or "fa".
#' @param orfs character string giving the name of tab-delimited TXT file with
#'     coordinates of main open reading frame(-s). It is usually the output of
#'     the function findAllORFs(). This file must be in the working directory.
#'     Allowed file format is "txt". This file should include at least four
#'     mandatory fields:
#'     i) transcript_id - ID of sequence;
#'     ii) orf_start    - start coordinate of open reading frame in a sequence;
#'     iii) orf_end     - end coordinate of open reading frame in a sequence;
#'     iv) orf_prob     - probability that given open reading frame is true
#'                        open reading frame.
#' @param prThr probability threshold for the "winning" class of main open
#'     reading frame(-s). The default value is 0.9.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return list containing data frame and object of class DNAStringSet. The
#'     data frame includes six fields:
#'     i) sequence_id   - ID of sequence;
#'     ii) orf_id       - ID of upstream open reading frame;
#'     iii) orf_start   - start coordinate of upstream open reading frame in a
#'                        sequence;
#'     iv) orf_end      - end coordinate of upstream open reading frame in a
#'                        sequence;
#'     v) prf_length    - length of upstream open reading frame;
#'     vi) orf_sequence - sequence of upstream open reading frame.
#'     DNAStringSet object contains the sequence(-s) of identified variant(-s)
#'     of upstream open reading frame(-s).
#' @authors Vasily V. Grinev
#' @examples
#' tr <- "Kasumi-1, transcriptome, StringTie, filtered, alternative assemblies.fasta"
#' orfs <- "Kasumi-1, transcriptome, StringTie, filtered, alternative assemblies, ORFs.txt"
#' res <- findUpstrORFs(tr=tr,
#'                      orfs=orfs,
#'                      prThr=0.9,
#'                      workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 17, 2023.

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
