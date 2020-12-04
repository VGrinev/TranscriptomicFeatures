#' Load a set of transcripts
#' @description Load a set of experimental transcripts.
#' @param tr character string giving the name of file with experimental
#'     transcripts. Allowed file formats are "fasta", "fa", "gtf" or "gff".
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return List of loaded transcript sequences.
#' @author Vasily V. Grinev
#' @examples
#' trans <- system.file("extdata",
#'                      "Set.trans_sequences.fasta",
#'                      package="ORFhunteR")
#' trans_seq <- loadTrExper(tr=trans)
#' @export

loadTrExper <- function(tr,
                        genome="BSgenome.Hsapiens.UCSC.hg38",
                        workDir=NULL){
### Validation of file format.
##  Retrieving the file extension.
frt <- tools::file_ext(x=tr)
##  Validation of file format.
if (!frt %in% c("fasta", "fa", "gtf", "gff")){
    stop("Invalid file format")
}
### Full path to the file with experimental transcripts.
if (!is.null(workDir)){
    tr <- paste(workDir, tr, sep="/")
}
### Loading of experimental transcripts in FASTA/FA format.
if (frt == "fasta" || frt == "fa"){
##  Loading of experimental transcripts as a list of character strings.
    trSeq <- as.list(x=as.character(x=readDNAStringSet(filepath=tr,
                                                       format=frt)))
}
### Loading of experimental transcripts in GTF/GFF format.
if (frt == "gtf" || frt == "gff"){
#   Loading of the transcript(-s) as an object of class GRanges.
    trSeq <- sort(x=import(con=tr, format=frt))
    trSeq <- trSeq[trSeq$type == "exon", ]
#   Converting a GRanges object into transcript-splitted GRangesList object.
    trSeq <- split(x=trSeq, f=trSeq$transcript_id)
#   Retrieving the sequence(-s).
    trSeq <- lapply(X=trSeq,
                    FUN=function(y){if (as.character(strand(y)@values) == "-"){
                          as.character(x=unlist(x=rev(x=getSeq(x=get(x=genome),
                                                               names=y))))
                                   }else{
                          as.character(x=unlist(x=getSeq(x=get(x=genome),
                                                         names=y)))
                                   }
                    }
             )
}
return(trSeq)
}
