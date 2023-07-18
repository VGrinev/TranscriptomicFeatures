#' Extract all unique k-mer from biological string(-s)
#' @description Extraction of all unique k-mer from biological string(-s).
#' @param x character string giving the name of file with the sequence(-s) of
#'     interest. This file must be in the working directory. Allowed file
#'     formats are "fasta" or "fa".
#' @param type type of sequence(-s): "DNA", "RNA" or "AA".
#' @param k width of k-mer.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return An object of class XStringSet.
#' @authors Vasily V. Grinev
#' @examples
#' x <- "AML proteome, Kasumi-1, June 23, 2023.fasta"
#' res <- extractKmer(x=x,
#'                    type="AA",
#'                    k=8,
#'                    workDir="D:/Vasily Grinev")
#' @export
#' Last updated: June 24, 2023.

extractKmer <- function(x,
                        type,
                        k,
                        workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with packages Biostrings v.2.60.2 and
    #   stringr v.1.4.0.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=stringr))
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
    ### Extraction of k-mer.
    mers <- lapply(X=seqs,
                   FUN=function(y){
                            y=str_sub(string=y,
                                      start=seq(from=1, to=nchar(x=y) - k + 1),
                                      end=seq(from=k, to=nchar(x=y)))
                            y=unique(x=y)
                        }
                   )
    if (type == "DNA"){
        mers <- DNAStringSet(x=unlist(x=mers))
    }
    if (type == "RNA"){
        mers <- RNAStringSet(x=unlist(x=mers))
    }
    if (type == "AA"){
        mers <- AAStringSet(x=unlist(x=mers))
    }
    mers <- mers[order(substring(text=names(x=mers), first=1, last=14)), ]
    IDs <- unlist(x=apply(X=data.frame(table(substring(text=names(x=mers),
                                                                  first=1,
                                                                  last=14)
                                             )
                                       ),
                           MARGIN=1,
                           FUN=function(y){paste(y[1],
                                                 formatC(x=seq(from=1,
                                                               to=y[2]),
                                                         width=4,
                                                         flag="0"),
                                                 sep=".")}))
    names(x=mers) <- IDs
    ### Returning the final object.
    return(mers)
}
