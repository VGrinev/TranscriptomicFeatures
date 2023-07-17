#' Filtration of the Cufflinks/StringTie assembled transcripts
#' @description Filtrate of the Cufflinks/StringTie assembled transcripts.
#' @param gtf name of the input file in GTF/GFF format. It is typically output
#'     file of Cufflinks or StringTie.
#' @param format format of the output file (typically "gtf" or "gff").
#' @param unstr a logical argument indicating that the unstranded transcripts
#'     should be removed. Default value is TRUE.
#' @param difStr a logical argument specifying that the transcripts that match
#'     two different strands should be removed. Default value is TRUE.
#' @param canChr a logical argument specifying that the records from
#'     non-canonical chromosomes/scaffolds should be removed. Default value is
#'     TRUE.
#' @param oneEx a logical argument indicating that the one-exon transcripts
#'     should be removed. Default value is FALSE.
#' @param exL an integer argument. It is a threshold for minimal length of
#'     exon(-s) in transcripts. Default Ensembl-based value is 25.
#' @param intrL an integer argument. It is a threshold for minimal length of
#'     intron(-s) in transcripts. Default Ensembl-based value is 50.
#' @param trL an integer argument. It is a threshold for minimal length of
#'     transcript. Default value is 300.
#' @param samples named list of samples to be sub-setted during filtration
#'     against transcripts with too low abundance.
#' @param subsetting character string giving the way for subsetting of data
#'     matrix during determination of the expressed transcripts. Default value
#'     is "by_groups". An alternative way is "by_individual", that is, the
#'     determination of the expressed transcripts for each sample separately.
#' @param fpkmData name of the file with FPKM data on RNA isoforms expression.
#' @param fpkm an integer argument. It is a threshold for minimal transcripts
#'     abundance (in FPKM). Default value is 1.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return an object of class GRanges containing filtered transcripts.
#' @author Vasily V. Grinev
#' @examples
#' fpkmData="Kasumi-1, transcriptome, StringTie, primary set of transcripts.txt"
#' res <- filtrateTrans(gtf="Kasumi-1, transcriptome, StringTie, primary.gtf",
#'                      format="gtf",
#'                      unstr=TRUE,
#'                      difStr=TRUE,
#'                      canChr=TRUE,
#'                      oneEx=FALSE,
#'                      exL=25,
#'                      intrL=50,
#'                      trL=0,
#'                      samples=list(ctrl=c("GSM1316401",
#'                                          "GSM1316402",
#'                                          "GSM1316403"),
#'                                   expr=c("GSM1316404",
#'                                          "GSM1316405",
#'                                          "GSM1316406")),
#'                      subsetting="by_groups",
#'                      fpkmData=fpkmData,
#'                      fpkm=0.1,
#'                      workDir="D:/Vasily Grinev")
#' @export
#' Last updated: June 14, 2023.

filtrateTrans <- function(gtf,
                          format="gtf",
                          unstr=TRUE,
                          difStr=TRUE,
                          canChr=TRUE,
                          oneEx=FALSE,
                          exL=25,
                          intrL=50,
                          trL=300,
                          samples,
                          subsetting="by_groups",
                          fpkmData,
                          fpkm=1,
                          workDir=NULL){
    ### Loading of the required package.
    #   This code was successfully tested with package rtracklayer v.1.44.3.
    suppressMessages(expr=library(package=rtracklayer))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of primary data as a GRanges object.
    tr <- import(con=paste(workDir, gtf, sep="/"), format=format)
    tr <- tr[tr$type == "exon", ]
    ### Removing of unstranded transcripts.
    if (unstr == TRUE){
        tr <- tr[!strand(tr) == "*", ]
    }
    ### Removing of transcripts that match two different strands.
    if (difStr == TRUE){
        d_str <- cbind(as.character(x=strand(x=tr)), tr$transcript_id)
        d_str <- table(d_str[!duplicated(x=d_str), ][, 2])
        tr <- tr[!tr$transcript_id %in% attr(x=d_str[d_str > 1],
                                             which="dimnames")[[1]], ]
    }
    ### Only canonical chromosomes.
    if (canChr == TRUE){
        tr <- as.data.frame(tr)
        chr <- c(paste("chr", c(1:22, "X", "Y"), sep=""))
        tr <- tr[tr$seqnames %in% chr, ][, -8:-9]
        rownames(x=tr) <- NULL
        tr <- makeGRangesFromDataFrame(df=tr, keep.extra.columns=TRUE)
    }
    ### Removing of one-exon transcripts.
    if (oneEx == TRUE){
        tr <- tr[tr$transcript_id %in%
                 attr(x=table(x=tr$transcript_id)[table(x=tr$transcript_id) >
                                                  length(x=unique(x=tr$type))],
                      which="dimnames")[[1]], ]
    }
    ### Removing of transcripts with too short exon(-s).
    if (exL > 0){
        tr <- tr[!tr$transcript_id %in%
                 unique(x=tr[width(x=tr) < exL]$transcript_id), ]
    }
    ### Removing of transcripts with too short intron(-s).
    if (intrL > 0){
        introns <- split(x=tr, f=tr$transcript_id)
        introns <- unlist(x=setdiff(x=range(introns), y=introns))
        tr <- tr[!tr$transcript_id %in%
                 unique(x=names(x=introns[width(x=introns) < intrL, ])), ]
    }
    ### Removing too short transcripts.
    if (trL > 0){
        tr_length <- aggregate(x=width(x=tr),
                               by=list(tr$transcript_id),
                               FUN=sum)
        tr=tr[tr$transcript_id %in% tr_length[tr_length$x >= trL, ]$Group.1, ]
    }
    ### Removing of transcripts with too low abundance.
    fpkm_data <- read.table(file=paste(workDir, fpkmData, sep="/"),
                            sep="\t",
                            header=TRUE,
                            quote="\"",
                            as.is=TRUE)
    fpkm_data <- fpkm_data[, c(as.numeric(x=colnames(x=fpkm_data) ==
                                                              "transcript_id"),
                               grep(pattern="FPKM.+",
                                    x=colnames(x=fpkm_data),
                                    ignore.case=TRUE))]
    colnames(x=fpkm_data) <- gsub(pattern="FPKM.",
                                  replacement="",
                                  colnames(x=fpkm_data))
    rownames(x=fpkm_data) <- fpkm_data[, 1]
    fpkm_data <- fpkm_data[, -1]
    if (subsetting =="by_groups"){
        groups <- lapply(X=samples,
                         FUN=function(y){fpkm_data[, y]})
        groups <- lapply(X=groups,
                         FUN=function(y){y[rowSums(x=y >= fpkm) == ncol(x=y), ]})
        groups <- lapply(X=groups,
                         FUN=function(y){sort(x=unique(x=rownames(x=y)))})
        IDs <- sort(x=unique(x=as.character(x=do.call(what=c, args=groups))))
        tr <- tr[tr$transcript_id %in% IDs, ]
    }else{
        IDs <- apply(X=fpkm_data,
                     MARGIN=2,
                     FUN=function(y){rownames(x=fpkm_data[y >= fpkm, ])})
        IDs <- sort(x=unique(x=as.character(x=do.call(what=c, args=IDs))))
        tr <- tr[tr$transcript_id %in% IDs, ]
    }
    ### Returning of the final object of class GRanges.
    return(tr)
}
