#' Load a set of RNA molecules
#' @description Loading of the set of RNA molecules.
#' @param tr character string giving the name of file with RNA molecules.
#'     Allowed file formats are "fasta", "fa", "gtf" or "gff".
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return An object of class DNAStringSet.
#' @author Vasily V. Grinev
#' @examples
#' trans <- system.file("extdata",
#'                      "Set.trans_sequences.fasta",
#'                      package="ORFhunteR")
#' trans_seq <- loadTrExper(tr=trans)
#' @export
#' Last updated: April 19, 2023.

loadTrExper <- function(tr,
                        genome="BSgenome.Hsapiens.UCSC.hg38",
                        workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with the package Biostrings v.2.60.2.
    suppressMessages(expr=library(package=Biostrings))
    source(file="D:/Software/Data Analytics Pipelines/Identification and Annotation of the Human ORFs/R/convertGTFtoFASTA.r")
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Retrieving the file extension.
    frt <- tools::file_ext(x=tr)
    ### Validation of file format.
    if (!frt %in% c("fasta", "fa", "gtf", "gff")){
        stop("Invalid file format")
    }
    ### Loading of RNA molecules in FASTA/FA format.
    if (frt == "fasta" || frt == "fa"){
        trSeq <- readDNAStringSet(filepath=paste(workDir, tr, sep="/"),
                                  format=frt)
    }
    ### Loading of RNA molecules in GTF/GFF format.
    if (frt == "gtf" || frt == "gff"){
        trSeq <- convertGTFtoFASTA(gtf=tr,
                                   format=frt,
                                   genome=genome,
                                   workDir=workDir)
    }
    return(trSeq)
}
