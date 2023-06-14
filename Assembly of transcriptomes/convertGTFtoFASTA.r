#' Conversion of GTF/GFF annotations into FASTA sequences
#' @description Convert the GTF/GFF annotations into FASTA sequences.
#' @param gtf character string giving the name of input file in GTF/GFF format.
#' @param format format of the input file (typically "gtf" or "gff").
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return DNAStringSet object with FASTA sequences.
#' @author Vasily V. Grinev
#' @examples
#' res <- convertGTFtoFASTA(gtf="Kasumi-1, transcriptome, StringTie, filtered.gtf",
#'                          format="gtf",
#'                          genome="BSgenome.Hsapiens.UCSC.hg38", 
#'                          workDir="D:/Vasily Grinev")
#' @export
#' Last updated: June 14, 2023.

convertGTFtoFASTA <- function(gtf,
                              format="gtf",
                              genome="BSgenome.Hsapiens.UCSC.hg38", 
                              workDir=NULL){
    ### Loading of the required packages.
    #   This code was successfully tested with packages rtracklayer v.1.44.3
    #   and BSgenome.Hsapiens.UCSC.hg38 v.1.4.3.
    suppressMessages(expr=library(package=rtracklayer))
    suppressMessages(expr=library(package=BSgenome.Hsapiens.UCSC.hg38))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of GTF/GFF data as a GRanges object.
    tr <- sort(x=import(con=paste(workDir, gtf, sep="/"), format=format))
    tr <- tr[tr$type == "exon", ]
    ### Converting the GRanges object into transcript-splitted GRangesList
    #   object.
    tr_seq <- split(x=tr, f=tr$transcript_id)
    ### Retrieving the sequence(-s).
    tr_seq <- lapply(X=tr_seq,
                     FUN=function(y){if(as.character(strand(y)@values) == "-"){
                                        unlist(x=rev(x=getSeq(x=get(x=genome),
                                                              names=y)))
                                     }else{
                                        unlist(x=getSeq(x=get(x=genome),
                                                              names=y))
                                     }
                         }
                     )
    tr_seq <- DNAStringSet(x=tr_seq)
    ### Returning of the final object of class DNAStringSet.
    return(tr_seq)
}
