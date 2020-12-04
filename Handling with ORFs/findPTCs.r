#' Identify the premature termination codons in nucleotide sequences
#' @description Identify the premature termination codons in nucleotide
#'     sequences of interest.
#' @param orfs character string giving the name of tab-delimited TXT file with
#'     coordinates of open reading frame(-s). This file should include four
#'     mandatory fields:
#'     i) transcript_id - transcript ID;
#'     ii) start        - start coordinate of open reading frame
#'                        in a transcript;
#'     iii) end         - end coordinate of open reading frame in a transcript;
#'     iv) length       - length of open reading frame.
#' @param gtf character string giving the name of GTF/GFF file with annotated
#'     transcripts of interest. Valid format is "gtf" or "gff".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return data.frame object with start and stop positions, length and stop
#'     status of codons for each transcript ID.
#' @author Vasily V. Grinev
#' @examples
#' orfs_path <- system.file("extdata",
#'                          "Set.trans_ORFs.coordinates.txt",
#'                          package="ORFhunteR")
#' gtf_path <- system.file("extdata",
#'                         "Set.trans_sequences.gtf",
#'                         package="ORFhunteR")
#' ptcs <- findPTCs(orfs=orfs_path,
#'                  gtf=gtf_path,
#'                  workDir=NULL)
#' @export

findPTCs <- function(orfs, gtf, workDir=NULL){
### Loading of the coordinates of identified open reading frames as an object
#   of class data.frame.
##  Full path to the file.
if (!is.null(x=workDir)){
    orfs <- paste(workDir, orfs, sep="/")
}
##  Loading of the coordinates.
coordORFs <- read.table(file=orfs, header=TRUE, quote="\"", as.is=TRUE)
coordORFs$tr_length <- 0
coordORFs$l.ex_width <- 0
coordORFs <- as.matrix(x=coordORFs)
### Loading of experimental transcripts in GTF/GFF format.
##  Validation of file format.
#   Retrieving the file extension.
frt <- tools::file_ext(x=gtf)
#   Validation of file format.
if (!frt %in% c("gtf", "gff")){
    stop("Invalid file format")
}
##  Loading of experimental transcripts as a list of character strings.
#   Full path to the file.
if (!is.null(x=workDir)){
    gtf <- paste(workDir, gtf, sep="/")
}
expTrans <- rtracklayer::import(con=gtf, format=frt)
expTrans <- sort(x=expTrans[expTrans$type == "exon", ])
expTrans <- split(x=expTrans, f=expTrans$transcript_id)
### Identification of PTC(-s).
cl <- makeCluster(spec=detectCores() - 1)
clusterExport(cl=cl,
              varlist=c("as.character", "c", "length", "return", 
                        "strand", "sum", "width"))
lengths <- parLapply(X=expTrans,
                     fun=function(y){
                       if (as.character(x=rtracklayer::strand(x=y))[1] == "-"){
                           return(c(sum(x=width(x=y)),
                                  width(x=y[1])))
                      }else{
                           return(c(sum(x=width(x=y)),
                                  width(x=y[length(x=y)])))
                      }
                     },
                     cl=cl)
stopCluster(cl=cl)
lengths <- do.call(what=rbind, args=lengths)
idx <- coordORFs[, 1] %in% rownames(x=lengths)
coordORFs[idx, 6:7] = lengths[coordORFs[idx, 1], 1:2]
coordORFs <- data.frame(coordORFs, stringsAsFactors=FALSE)
coordORFs[, 2:7] <- apply(X=coordORFs[, 2:7], MARGIN=2, FUN=as.numeric)
coordORFs$stop.status <- "mature"
ptc_ctrl <- coordORFs[, 3] - (coordORFs[, 6] - coordORFs[, 7] + 1) <= -50
if (any(ptc_ctrl)){
    coordORFs[ptc_ctrl, ]$stop.status <- "premature"
} 
coordORFs <- coordORFs[, c(1:5, 8)]
### Returning a final object.
return(coordORFs)
}
