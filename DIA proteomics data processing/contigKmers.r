#' Merge of unique k-mer into contig(-s)
#' @description Merging of unique k-mer into contig(-s).
#' @param x character string giving the name of file with the sequence(-s) of
#'     origin of k-mer. This file must be in the working directory. Allowed
#'     file formats are "fasta" or "fa".
#' @param y character string giving the name of file with the sequences of
#'     k-mer. This file must be in the working directory. Allowed file formats
#'     are "fasta" or "fa".
#' @param type type of sequence(-s): "DNA", "RNA" or "AA".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return An object of class IRanges.
#' @authors Vasily V. Grinev
#' @examples
#' x <- "Kasumi-1 proteome, June 23, 2023.fasta"
#' y <- "Kasumi-1 proteome, June 23, 2023, all specific 8-mers.fasta"
#' res <- contigKmers(x=x,
#'                    y=y,
#'                    type="AA",
#'                    workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 5, 2023.

contigKmers <- function(x,
                        y,
                        type="AA",
                        workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with package Biostrings v.2.60.2.
    suppressMessages(expr=library(package=Biostrings))
    ### Loading of the set of sequences.
    ##  Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=x)
    ##  Validation of file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Full path to the file.
    seqs <- paste(workDir, x, sep="/")
    ##  Loading of the sequences.
    if (type == "DNA"){
        seqs <- readDNAStringSet(filepath=seqs)
    }
    if (type == "RNA"){
        seqs <- readRNAStringSet(filepath=seqs)
    }
    if (type == "AA"){
        seqs <- readAAStringSet(filepath=seqs)
    }
    names(x=seqs) <- substring(text=names(x=seqs), first=4, last=16)
    ### Loading of the k-mer.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=y)
    ##  Validation of file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Full path to the file.
    mers <- paste(workDir, y, sep="/")
    ##  Loading of the sequences.
    if (type == "DNA"){
        mers <- readDNAStringSet(filepath=mers)
    }
    if (type == "RNA"){
        mers <- readRNAStringSet(filepath=mers)
    }
    if (type == "AA"){
        mers <- readAAStringSet(filepath=mers)
    }
    names(x=mers) <- substring(text=names(x=mers), first=1, last=13)
    ### Development of contigs.
    seqs <- seqs[names(seqs) %in% names(x=mers)]
    z <- function(y){mer=mers[names(x=mers) == names(x=seqs)[[y]]]
                     contig <- unlist(x=matchPDict(pdict=mer,
                                                   subject=seqs[[y]],
                                                   algorithm="naive-exact"))
                     contig <- intersect(range(contig), contig)}
    contigs <- lapply(X=seq_along(seqs), FUN=z)
    names(x=contigs) <- names(x=seqs)
    contigs <- unlist(x=as(object=contigs, Class="SimpleIRangesList"))
    ### Returning the final object.
    return(contigs)
}

setwd(dir="D:/Vasily Grinev")
x <- "Kasumi-1 proteome, June 23, 2023.fasta"
y <- "Kasumi-1 proteome, June 23, 2023, all specific 8-mers.fasta"
res <- contigKmers(x=x,
                   y=y,
                   type="AA",
                   workDir="D:/Vasily Grinev")
write.table(x=res,
            file="Kasumi-1 proteome, June 23, 2023, all specific 8-mers, contigs.txt",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)
