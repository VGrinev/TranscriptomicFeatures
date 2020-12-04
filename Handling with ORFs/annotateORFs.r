#' Annotate open reading frames
#' @description Annotate the open reading frames identified in nucleotide
#'     sequences of interest.
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
#'     file. Valid format is "fasta" or "fa".
#' @param gtf character string giving the name of GTF/GFF file with annotated
#'     transcripts of interest. Default value is NULL.
#' @param prts character string giving the name of FASTA file with sequences of
#'     in silico translated proteins encoded by identified open reading frames.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return data.frame object with columns:
#'  \item{transcript_id}{transcript ID.}
#'  \item{f_utr.length}{length of 5'UTRs.}
#'  \item{start.codon}{type of start codon.}
#'  \item{orf.start}{start coordinate of open reading frames.}
#'  \item{orf.stop}{stop coordinate of open reading frames.}
#'  \item{stop.codon}{type of stop codon.}
#'  \item{stop.status}{PTC status of stop codon.}
#'  \item{orf.length}{length of open reading frames.}
#'  \item{t_utr.length}{length of 3'UTRs.}
#'  \item{MW}{molecular weight.}
#'  \item{pI}{isoelectic point of a protein sequence.}
#'  \item{indexPPI}{potential protein interaction index.}
#' @author Vasily V. Grinev
#' @examples
#' orfs_path <- system.file("extdata",
#'                          "Set.trans_ORFs.coordinates.txt",
#'                          package="ORFhunteR")
#' tr_path <- system.file("extdata",
#'                        "Set.trans_sequences.fasta",
#'                        package="ORFhunteR")
#' gtf_path <- system.file("extdata",
#'                         "Set.trans_sequences.gtf",
#'                         package="ORFhunteR")
#' prts_path <- system.file("extdata",
#'                          "Set.trans_proteins.sequences.fasta",
#'                          package="ORFhunteR")
#' anno_orfs <- annotateORFs(orfs=orfs_path,
#'                           tr=tr_path,
#'                           gtf=gtf_path,
#'                           prts=prts_path,
#'                           workDir=NULL)
#' @export

annotateORFs <- function(orfs,
                         tr,
                         gtf=NULL,
                         prts,
                         workDir=NULL){
### Loading of the coordinates of identified open reading frames as an object
#   of class data.frame.
##  Full path to the file.
if (!is.null(x=workDir)){
    orfs_path <- paste(workDir, orfs, sep="/")
}
##  Loading of the coordinates.
coordORFs <- read.table(file=orfs_path, header=TRUE, quote="\"", as.is=TRUE)
coordORFs <- coordORFs[order(coordORFs$transcript_id), ]
### Creation of an object with annotations.
##  Length of 5'UTRs.
annoORFs <- coordORFs[, 1:2]
annoORFs$start <- annoORFs$start - 1
colnames(x=annoORFs) <- c("transcript_id", "f_utr.length")
##  Type of start codons.
#   Loading of nucleotide sequences of interest as an object of class 
#   DNAStringSet.
if (!is.null(x=workDir)){
    tr <- paste(workDir, tr, sep="/")
}
seqTrans <- readDNAStringSet(filepath=tr)
seqTrans <- seqTrans[order(x=names(x=seqTrans))]
#   Extraction of start codons.
annoORFs$start.codon <- substr(x=seqTrans,
                               start=coordORFs$start,
                               stop=coordORFs$start + 2)
##  Extraction of coordinates of open reading frames.
annoORFs$orf.start <- coordORFs$start
annoORFs$orf.stop <- coordORFs$end
##  Type of stop codon.
annoORFs$stop.codon <- substr(x=seqTrans,
                              start=annoORFs$orf.stop - 2, 
                              stop=annoORFs$orf.stop)
##  PTC status of stop codon.
if (!is.null(x=gtf)){
    ptc <- findPTCs(orfs=orfs, gtf=gtf, workDir=workDir)[, c(1, 6)]
    ptc <- ptc[order(x=ptc$transcript_id), ]
    annoORFs$stop.status <- ptc$stop.status
    ptc = NULL
}else{
    annoORFs$stop.status <- "ND"
}
### Length of open reading frames.
annoORFs$orf.length <- coordORFs$length
### Length of 3'UTRs.
annoORFs$t_utr.length <- width(x=seqTrans) - coordORFs$end
### Calculation the molecular weight of a protein sequence.
##  Loading of amino acid sequences of interest as an object of class
#   AAStringSet.
if (!is.null(x=workDir)){
    prts <- paste(workDir, prts, sep="/")
}
seqPrts <- readAAStringSet(filepath=prts)
##  Molecular weight calculation.
cl <- makeCluster(spec=detectCores() - 1)
clusterExport(cl=cl, varlist=c("mw", "pI", "boman"))
annoORFs$MW <- unlist(x=parLapply(X=seqPrts,
                                  fun=function(y){mw(seq=y)},
                                  cl=cl))/1000
annoORFs$MW <- round(x=as.vector(x=annoORFs$MW), digits=2)
### Calculation the isoelectic point of a protein sequence.
annoORFs$pI <- round(x=unlist(x=parLapply(X=seqPrts, 
                                          fun=function(y){pI(seq=y)},
                                          cl=cl)),
                     digits=2)
### Calculation the potential protein interaction index.
annoORFs$indexPPI <- unlist(x=parLapply(X=seqPrts, 
                                        fun=function(y){boman(seq=y)},
                                        cl=cl))
annoORFs$indexPPI <- round(x=annoORFs$indexPPI, digits=2)
stopCluster(cl=cl)
### Returning a final object.
return(annoORFs)
}
