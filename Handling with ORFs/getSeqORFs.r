#' Extract the sequences of identified open reading frames
#' @description Extract the sequences of identified open reading frames.
#' @param orfs character string giving the name of tab-delimited TXT file with
#'     coordinates of open reading frame(-s). This file should include three
#'     mandatory fields:
#'     i) transcript_id - transcript ID;
#'     ii) start        - start coordinate of open reading frame
#'                        in a transcript;
#'     iii) end         - end coordinate of open reading frame in a transcript.
#' @param tr character string giving the name of file with transcript(-s)
#'     of interest. This file must include the transcript(-s) for which the
#'     open reading frame(-s) was/were identified and listed in above orfs
#'     file. Valid formats are "gtf", "gff", "fasta" and "fa".
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return DNAStringSet object with sequences of extracted open reading frames.
#' @author Vasily V. Grinev
#' @examples
#' orfs_path <- system.file("extdata",
#'                          "Set.trans_ORFs.coordinates.txt",
#'                          package="ORFhunteR")
#' tr_path <- system.file("extdata",
#'                        "Set.trans_sequences.fasta",
#'                        package="ORFhunteR")
#' seq_orfs <- getSeqORFs(orfs=orf_path,
#'                        tr=trans_path,
#'                        genome="BSgenome.Hsapiens.UCSC.hg38",
#'                        workDir=NULL)
#' @export

getSeqORFs <- function(orfs,
                       tr,
                       genome="BSgenome.Hsapiens.UCSC.hg38", 
                       workDir=NULL){
### Loading of the coordinates of identified open reading frames as an object
#   of class data.frame.
##  Full path to the file.
if (!is.null(workDir)){
    orfs <- paste(workDir, orfs, sep="/")
}
##  Loading of the coordinates.
coordORFs <- read.table(file=orfs, header = TRUE, quote = "\"", as.is = TRUE)
coordORFs <- coordORFs[order(coordORFs$transcript_id), ]
### Loading of the transcript(-s) as an object of class DNAStringSet.
##  Retrieving the file extension.
frt <- tools::file_ext(x=tr)
##  Validation of file format.
if (!frt %in% c("fasta", "fa", "gtf", "gff")){
    stop("Invalid file format")
}
##  Full path to the file with the transcript(-s).
if (!is.null(workDir)){
    tr <- paste(workDir, tr, sep="/")
}
##  Loading of the transcript(-s) in FASTA/FA format.
if (frt == "fasta" || frt == "fa"){
    trSeq <- readDNAStringSet(filepath=tr, format=frt)
    trSeq <- trSeq[order(names(x=trSeq))]
}
##  Loading of the transcript(-s) in GTF/GFF format.
if (frt == "gtf" || frt == "gff"){
#   Loading of the transcript(-s) as an object of class GRanges.
    trSeq <- sort(x=import(con=tr, format=frt))
    trSeq <- trSeq[trSeq$type == "exon", ]
#   Converting a GRanges object into transcript-splitted GRangesList object.
    trSeq <- split(x=trSeq, f=trSeq$transcript_id)
#   Retrieving the sequence(-s).
    trSeq <- lapply(X=trSeq,
                    FUN=function(y){if (as.character(strand(y)@values) == "-"){
                                        unlist(x=rev(x=getSeq(x=get(x=genome),
                                                              names=y)))
                                   }else{
                                        unlist(x=getSeq(x=get(x=genome),
                                                        names=y))
                                   }
                    }
             )
    trSeq <- DNAStringSet(x=trSeq)
}
### Extraction the sequences of identified open reading frames.
seqORFs <- DNAStringSet(x=substr(x=trSeq,
                                 start=coordORFs$start,
                                 stop=coordORFs$end),
                        use.names=TRUE)
### Return the final object of class DNAStringSet.
return(seqORFs)
}
