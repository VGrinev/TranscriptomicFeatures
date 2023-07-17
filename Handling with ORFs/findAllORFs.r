#' Find all variants of open reading frames in nucleotide sequence(-s)
#' @description Identification all possible variants of open reading frames in
#'     nucleotide sequence(-s) of interest.
#' @param tr character string giving the name of file with the sequences of
#'     interest. Usually it is set of sequences of RNA molecules. This file
#'     must be in the working directory. Allowed file formats are "fasta" or
#'     "fa". Alternatively, it may be a ready-to-use object of class
#'     DNAStringSet stored in the computer's RAM.
#' @param codStart character string with type of start codon: "ATG", "GTG",
#'     "TTG" or "CTG". Default value is "ATG".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return list containing data frame and object of class DNAStringSet. The
#'     data frame includes six fields:
#'     i) sequence_id   - ID of sequence;
#'     ii) orf_id       - ID of open reading frame;
#'     iii) orf_start   - start coordinate of open reading frame in a sequence;
#'     iv) orf_end      - end coordinate of open reading frame in a sequence;
#'     v) prf_length    - length of open reading frame;
#'     vi) orf_sequence - sequence of open reading frame.
#'     DNAStringSet object contains the sequence(-s) of identified variant(-s)
#'     of open reading frame(-s).
#' @author Vasily V. Grinev
#' @examples
#' tr <- "Kasumi-1, transcriptome, StringTie, filtered, all transcripts.fasta"
#' orf <- findAllORFs(tr=tr, codStart="ATG", workDir="D:/Vasily Grinev")
#' @export
#' Last updated: June 21, 2023.

findAllORFs <- function(tr, codStart="ATG", workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with packages Biostrings v.2.60.2 and
    #   GenomicRanges v.1.40.0.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=GenomicRanges))
    ### Loading of the nucleotide sequence(-s) of interest as an object of
    #   class DNAStringSet.
    if (class(x=tr) == "character"){
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
        seqs <- paste(workDir, tr, sep="/")
        ##  Loading of the sequence(-s) of interest.
        seqs <- readDNAStringSet(filepath=seqs)
        seqs <- seqs[order(x=names(x=seqs)), ]
    }else{
        seqs <- tr
    }
    ### Identification of all possible locations of start and stop codons.
    tr_start <- unlist(x=vmatchPattern(pattern=codStart, subject=seqs))
    tr_start <- data.frame(tr_start)
    tr_start <- tr_start[, c(4, 1:2)]
    colnames(x=tr_start) <- c("seqnames", "start", "end")
    tr_stop <- unlist(x=c(vmatchPattern(pattern="TAA", subject=seqs),
                          vmatchPattern(pattern="TAG", subject=seqs),
                          vmatchPattern(pattern="TGA", subject=seqs)))
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
    seqs <- seqs[names(x=seqs) %in% names(x=tr_start)]
    tr_start <- tr_start[names(x=tr_start) %in% names(x=seqs)]
    tr_stop <- tr_stop[names(x=tr_stop) %in% names(x=seqs)]
    ### Identification of all possible variants of open reading frames.
    ORFs <- list()
    for (i in 1:length(x=tr_start)){
        inFrame <- outer(X=tr_stop[[i]]$start,
                         Y=tr_start[[i]]$start,
                         FUN="-")/3
        inFrame <- which(x=round(x=inFrame) == inFrame & inFrame > 0,
                         arr.ind=TRUE)
        if (nrow(x=inFrame) > 0){
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
            if (nrow(x=inFrame) > 0){
                inFrame[, 2] <- inFrame[, 2] + 2
                inFrame <- cbind(names(x=tr_stop[i]),
                                 inFrame,
                                 substring(text=as.character(x=seqs[[i]]),
                                           first=inFrame[, "x"],
                                           last=inFrame[, "Group.1"]))
                colnames(x=inFrame) <- c("seqnames",
                                         "start", "end",
                                         "orf.sequence")
                rownames(x=inFrame) <- NULL
                ORFs[[i]] <- inFrame
            }
        }
    }
    ORFs <- do.call(what=rbind, args=ORFs)
    orf_id <- unlist(x=sapply(X=table(ORFs$seqnames),
                              FUN=function(y){paste(".ORF",
                                                    formatC(x=seq(from=1,
                                                                  to=y),
                                                            width=3,
                                                            flag="0"),
                                                    sep="")}))
    ORFs$orf_id <- paste(ORFs$seqnames, as.vector(x=orf_id), sep="")
    ORFs$length <- ORFs$end - ORFs$start + 1
    ORFs <- ORFs[, c(1, 5, 2:3, 6, 4)]
    colnames(ORFs) <- c("sequence_id", "orf_id",
                        "orf_start", "orf_end", "orf_length",
                        "orf_sequence")
    ORFs <- list(ORFs=ORFs, seqs=DNAStringSet(x=ORFs[, 6]))
    names(x=ORFs$seqs) <- ORFs$ORFs$orf_id
    ### Returning the final object.
    return(ORFs)
}
