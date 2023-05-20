#' Parse of blastp hits for proteins encoded by circRNAs
#' @description Parsing of the blastp hits for proteins encoded by circRNAs.
#' @param msProts character string giving the name of input file with
#'     sequences of the mass spectrometry (MS) approved circRNAs encoded
#'     proteins. This file must be in the working directory. Valid file format
#'     can be either "fasta" or "fa".
#' @param blastpHits character string giving the name of input file with raw
#'     blastp hits. This file must be in the working directory. Valid file
#'     format is "txt". It is typically hit table of blastp.
#' @param thr1 integer value giving the threshold for identification of protein
#'     stretch(-s) with low identity. Default value is 90.
#' @param msPepts character string giving the name of file with the annotations
#'     of MS detected circRNAs encoded peptides. This file must be in the
#'     working directory. Valid file format is "txt". It is typically output
#'     file of function identifyBSJPeptidesFromMS().
#' @param thr2 integer value giving the minimal number of MS detected
#'     peptide(-s) supporting the stretch(-s) with low identity. Default value
#'     is 1.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return The data frame containing the following fields:
#'     i) protein_id              - protein ID;
#'     ii) fragment_id            - ID of new peptide (fragment) encoded by
#'                                  circRNA;
#'     iii) fragment_start        - start coordinate (in protein) new peptide
#'                                  (fragment) encoded by circRNA;
#'     iv) fragment_end           - end coordinate (in protein) new peptide
#'                                  (fragment) encoded by circRNA;
#'     v) fragment_seq            - the sequence of new peptide (fragment)
#'                                  encoded by circRNA;
#'     vi) fragment_peptides_n    - number of MS detected peptide(-s)
#'                                  supporting of new peptide (fragment)
#'                                  encoded by circRNA;
#'     vii) fragment_peptides_seq - the sequence(-s) of MS detected peptide(-s)
#'                                  supporting the new peptide (fragment)
#'                                  encoded by circRNA.
#' @author Vasily V. Grinev
#' @examples
#' msProts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, sequences.fasta"
#' blastpHits="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, blastp hits.txt"
#' thr1=90
#' msPepts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, MS peptides, annotations.txt"
#' thr2=1
#' workDir="D:/Vasily Grinev"
#' res <- parseBlastPhitsCircRNAproteins(msProts=msProts, blastpHits=blastpHits, thr1=thr1, msPepts=msPepts, thr2=thr2, workDir=workDir)
#' @export
#' Last updated: May 20, 2023.

parseBlastPhitsCircRNAproteins <- function(msProts,
                                           blastpHits,
                                           thr1=90,
                                           msPepts,
                                           thr2=1,
                                           workDir=NULL){
    ### Loading of the required packages.
    #   This code was successfully tested with packages Biostrings v.2.60.2 and
    #   GenomicRanges v.1.46.1.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=GenomicRanges))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of the in silico predicted proteins as an AAStringSet object.
    proteins <- readAAStringSet(filepath=paste(workDir, msProts, sep = "/"))
    ### Loading of the blastp output as an object of class list.
    hits <- read.table(file=paste(workDir, blastpHits, sep = "/"),
                       sep="\t",
                       header=FALSE,
                       quote="\"",
                       as.is=TRUE)
    colnames(x=hits) <- c("query",
                          "subject",
                          "identity",
                          "alignment_length", "mismatches", "gap_opens",
                          "q.start", "q.end",
                          "s.start", "s.end",
                          "e_value",
                          "bit_score",
                          "positives")
    hits <- split(x=hits, f=hits$query)
    ### Parsing of blastp hits.
    identity <- list()
    for (i in 1:length(x=hits)){
        protein <- proteins[names(x=proteins) == names(x=hits[i])]
        query <- GRanges(seqnames=names(x=protein),
                         ranges=IRanges(start=1:width(x=protein),
                                        end=1:width(x=protein)))
        query$identity <- 0
        subject <- GRanges(seqnames=names(x=protein),
                           ranges=IRanges(start=hits[[i]]$q.start,
                                          end=hits[[i]]$q.end),
                           mcols=DataFrame(hits[[i]]))
        subject <- subject[order(start(x=subject), end(x=subject)), ]
        overlaps <- findOverlaps(query=query, subject=subject, type="any")
        overlaps <- cbind(data.frame(query[queryHits(x=overlaps), ])[, 1:2],
                          data.frame(subject[subjectHits(x=overlaps), 3])[, 6])
        overlaps <- aggregate(x=overlaps[, 3], by=list(overlaps[, 2]), FUN=max)
        query <- data.frame(query)[, c(1:2, 6)]
        query[overlaps$Group.1, ]$identity <- overlaps$x
        if (length(x=unique(x=overlaps[, 2])) > 1){
            ident <- which(x=query[, 3] <= thr1)
                if (length(x=ident) > 0){
                    ident <- tapply(X=ident,
                                    INDEX=cumsum(x=c(1,
                                                 diff(x=ident)) - 1 >= 8),
                                    FUN=range)
                    ident <- cbind(names(x=protein),
                                   data.frame(do.call(what=rbind,
                                                      args=ident)))
                    colnames(x=ident) <- c("protein_id", "start", "end")
                    ident <- ident[(ident$end - ident$start + 1) >= 8, ]
                    identity[[i]] <- ident
                }
        }
    }
    identity <- do.call(what=rbind, args=identity)
    rownames(x=identity) <- NULL
    ### Loading of the MS detected peptides as an object of class data frame.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=msPepts)
    ##  Validation of the file format.
    if (!frt %in% c("txt")){
        stop("Invalid file format")
    }
    ##  Loading of the peptides.
    pe <- read.table(file=paste(workDir, msPepts, sep = "/"),
                     sep="\t",
                     header=TRUE,
                     quote="\"",
                     as.is=TRUE)[, c(1, 4)]
    ### Developing of peptide map.
    peptideMap <- list()
    for (i in 1:nrow(x=identity)){
        protein <- proteins[names(x=proteins) %in% identity$protein_id[i]][[1]]
        pepts <- pe[pe$protein_id %in% identity$protein_id[i], ]$peptides
        pepts <- AAStringSet(x=strsplit(x=pepts, split=",")[[1]])
        query <- GRanges(seqnames=identity$protein_id[i],
                         ranges=IRanges(start=identity$start[i],
                                        end=identity$end[i]))
        subject <- sort(x=unlist(x=matchPDict(pdict=pepts,
                                              subject=protein)))
        subject <- GRanges(seqnames=identity$protein_id[i],
                           ranges=subject)
        overlaps <- findOverlaps(query=query,
                                 subject=subject,
                                 type="any",
                                 minoverlap=3L)
        peptideMap[[i]] <- overlaps
    }
    names(x=peptideMap) <- identity$protein_id
    peptideMap <- peptideMap[lengths(x=peptideMap) >= thr2]
    ### Developing a final object with annotations and subsequences.
    subseqs <- identity[identity$protein_id %in% names(x=peptideMap), ]
    colnames(x=subseqs) <- c("protein_id", "fragment_start", "fragment_end")
    subseqs$fragment_peptides_n <- 0
    subseqs$fragment_peptides_seq <- ""
    subseqs$fragment_id <- paste(paste(paste(subseqs$protein_id,
                                             "_fr.",
                                             sep=""),
                                       subseqs$fragment_start,
                                       sep=""),
                                 subseqs$fragment_end,
                                 sep="-")
    subseqs$fragment_seq <- ""
    subseqs <- subseqs[, c(1, 6, 2:3, 7, 4:5)]
    for (i in 1:nrow(x=subseqs)){
        protein <- proteins[names(x=proteins) %in% subseqs$protein_id[i]][[1]]
        pepts <- pe[pe$protein_id %in% gsub(pattern="_fr.+",
                                            replacement="",
                                            x=subseqs$protein_id[i]), ]$peptides
        pepts <- AAStringSet(x=strsplit(x=pe[pe$protein_id %in%
                                             subseqs$protein_id[i], 2],
                                        split=",")[[1]])
        query <- GRanges(seqnames=subseqs$protein_id[i],
                         ranges=IRanges(start=subseqs$fragment_start[i],
                                        end=subseqs$fragment_end[i]))
        subject <- sort(x=unlist(x=matchPDict(pdict=pepts,
                                              subject=protein)))
        subject <- GRanges(seqnames=subseqs$protein_id[i],
                           ranges=subject)
        overlaps <- findOverlaps(query=query,
                                 subject=subject,
                                 type="any",
                                 minoverlap=3L)
        if (length(x=overlaps) > 0){
            pepts <- substring(text=as.character(x=protein),
                               first=start(x=subject[subjectHits(x=overlaps), ]),
                               last=end(x=subject[subjectHits(x=overlaps), ]))
            pepts <- sort(x=unique(x=pepts))
            subseqs$fragment_peptides_n[i] <- length(x=pepts)
            subseqs$fragment_peptides_seq[i] <- paste(pepts, collapse=",")
            sub_str <- subseq(x=protein,
                              start=subseqs$fragment_start[i],
                              end=subseqs$fragment_end[i])
            subseqs$fragment_seq[i] <- as.character(x=sub_str)
        }
    }
    subseqs <- subseqs[subseqs$fragment_peptides_n > 0, ]
    anno <- subseqs
    subseqs <- AAStringSet(x=subseqs$fragment_seq)
    names(x=subseqs) <- anno$fragment_id
    ### Returning of the final results.
    return(list(annotations=anno, "fragment sequences"=subseqs))
}
