#' Detect cis-splicing events between adjacent and/or overlapped genes
#' @description Detection of cis-splicing events between adjacent and/or
#'     overlapped genes based on subjunc's outputs.
#' @param anno character vector giving the name of tab-delimited TXT file with
#'     annotated set of genes. This file should include five mandatory fields:
#'     i) gene_id   - IDs of gene;
#'     ii) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     iii) start   - start genomic coordinate of the gene;
#'     iv) end      - end genomic coordinate of the gene;
#'     v) strand    - strand information about gene location.
#' @param HRs character vector giving the name of subjunc's output file in
#'     format BED (postfix .junction.bed) or VCF (postfix .breakpoints.vcf) 
#'     with inferred canonical and/or non-canonical splicing junctions.
#' @param thr_type character vector giving the type of threshold for minimal
#'     coverage (sequencing depth) for junction to be considered detected.
#'     "CPM" by default. Alternative value is "RAW" that means direct usage of
#'     raw reads (counts).
#' @param thr_value an integer argument. It is a threshold value for minimal
#'     coverage (sequencing depth) for junction to be considered detected.
#'     Default value is 1 but it can be specified by the user.
#' @param canChr a logical argument specifying that the records from
#'     non-canonical chromosomes/scaffolds should be removed. Default value
#'     is TRUE.
#' @param minLength an integer argument giving the minimal length (distance) of
#'     splicing event. Default value is 50.
#' @param maxLength an integer argument giving the maximal length (distance) of
#'     splicing event. Default value is 123461. This default value corresponds
#'     to the 99th quantile of the intron length distribution from Ensembl
#'     annotation of human genome/transcriptome (Ensembl release 114 based on
#'     GRCh38.p14 human reference genome assembly).
#' @param parlgs character vector giving the name of tab-delimited TXT file
#'     with pairs of paralogous genes. This file should include two fields:
#'     i) gene1  - identifier of the first gene in the pair;
#'     ii) gene2 - Identifier of the second gene in the pair.
#'     The parlgs argument defaults to NULL, meaning that gene pairs will not
#'     be annotated by paralogous information. 
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return data frame with inferred cis-SAGe.
#' @author Vasily V. Grinev.
#' @examples
#' anno="Ensembl release 114, GRCh38.p14, annotations of genes.txt"
#' res <- detectCisSAGe(anno=anno,
#'                      HRs="GSM1316401.bed",
#'                      thr_type="CPM",
#'                      thr_value=1,
#'                      canChr=TRUE,
#'                      minLength=50,
#'                      maxLength=123461,
#'                      parlgs=NULL,
#'                      workDir="D:/Vasily Grinev")
#' @export
#' Last updated: June 27, 2025.

detectCisSAGe <- function(anno,
                          HRs,
                          thr_type="CPM",
                          thr_value=1,
                          canChr=TRUE,
                          minLength=50,
                          maxLength=123461,
                          parlgs=NULL,
                          workDir=NULL){
    ### Loading of the required packages.
    #   This code was successfully tested with packages GenomicRanges v.1.40.0,
    #   stringr v.1.5.1, edgeR v.3.30.3 and data.table v.1.15.4.
    suppressMessages(expr=library(package=GenomicRanges))
    suppressMessages(expr=library(package=tools))
    suppressMessages(expr=library(package=stringr))
    suppressMessages(expr=library(package=edgeR))
    suppressMessages(expr=library(package=data.table))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of the annotation data as an object of class GRanges.
    GRs <- read.table(file=paste(workDir, anno, sep="/"),
                      header=TRUE,
                      sep="\t",
                      quote="\"",
                      as.is=TRUE)
    GRs <- makeGRangesFromDataFrame(df=GRs,
                                    keep.extra.columns=TRUE)
    ### Loading of the hybrid RNA data as an object of class data frame.
    #   Handling with VCF file.
    if (file_ext(x=HRs) == "vcf"){
        hEEJs <- read.table(file=paste(workDir, HRs, sep="/"),
                            header=FALSE,
                            sep="\t",
                            quote="\"",
                            as.is=TRUE)
        hEEJs <- hEEJs[seq(from=1, to=nrow(x=hEEJs) - 1, by=2), ]
        hEEJs$V3 <- gsub(pattern=":.+",
                         replacement="",
                         x=gsub(pattern=".+chr",
                                replacement="chr",
                                x=hEEJs$V5))
        hEEJs$V4 <- str_extract_all(string=gsub(pattern=".+:",
                                                replacement="",
                                                x=hEEJs$V5),
                                    pattern="\\d+")
        hEEJs$V4 <- as.numeric(x=unlist(x=hEEJs$V4))
        hEEJs$V5 <- as.numeric(x=gsub(pattern=".+SR=",
                                      replacement="",
                                      x=hEEJs$V8))
        hEEJs <- hEEJs[, 1:5]
        if (thr_type == "CPM"){
            hEEJs$V5 <- as.vector(x=cpm(y=hEEJs$V5, log=TRUE, prior.count=0))
        }
        hEEJs <- hEEJs[hEEJs$V5 >= thr_value, ]
        if (canChr == TRUE){
            hEEJs <- hEEJs[hEEJs$V1 %in%
                                  c(paste("chr", c(1:22, "X", "Y"), sep="")), ]
        }
        hEEJs <- hEEJs[hEEJs$V1 == hEEJs$V3, ]
        hEEJs <- hEEJs[hEEJs$V4 - hEEJs$V2 >= minLength, ]
        if (!is.null(x=maxLength)){
            hEEJs <- hEEJs[hEEJs$V4 - hEEJs$V2 <= maxLength, ]
        }
        rownames(x=hEEJs) <- NULL
        hEEJs$id <- paste(paste(hEEJs$V1, hEEJs$V2, sep=":"), hEEJs$V4, sep="-")
        hEEJs <- hEEJs[, c(6, 1, 2, 4, 5)]
        colnames(x=hEEJs) <- c("id", "seqnames", "start", "end", "sr")
    }
    #   Handling with BED file.
    if (file_ext(x=HRs) == "bed"){
        hEEJs <- read.table(file=paste(workDir, HRs, sep="/"),
                            header=FALSE,
                            sep="\t",
                            quote="\"",
                            as.is=TRUE)
        hEEJs <- data.frame(cbind(hEEJs$V1,
           hEEJs$V2 + as.numeric(x=do.call(what=rbind, 
                                           args=strsplit(x=hEEJs$V11,
                                                         split=","))[, 1]),
           hEEJs$V3 - as.numeric(x=do.call(what=rbind,
                                           args=strsplit(x=hEEJs$V11,
                                                         split=","))[, 2]) + 1,
                        hEEJs$V5), stringsAsFactors=FALSE)
        hEEJs$id <- paste(paste(hEEJs$X1,
                                hEEJs$X2,
                                sep=":"),
                          hEEJs$X3,
                          sep="-")
        hEEJs <- hEEJs[, c(5, 1:4)]
        hEEJs[, 3:5] <- apply(X=hEEJs[, 3:5], MARGIN=2, FUN=as.numeric)
        colnames(hEEJs) <- c("id", "seqnames", "start", "end", "sr")
        if (length(x=unique(x=hEEJs$id)) < length(x=hEEJs$id)){
            sr <- data.table(hEEJs[, c(1, 5)])
            sr <- data.frame(sr[, list(sum(x=sr)), by="id"])
            sr <- sr[order(x=sr$id), ]
            hEEJs <- hEEJs[, -5]
            hEEJs <- hEEJs[!duplicated(x=hEEJs), ]
            hEEJs <- hEEJs[order(x=hEEJs$id), ]
            hEEJs$sr <- sr$V1
        }
        if (thr_type == "CPM"){
            hEEJs$sr <- as.vector(x=cpm(y=hEEJs$sr, log=TRUE, prior.count=0))
        }
        hEEJs <- hEEJs[hEEJs$sr >= thr_value, ]
        if (canChr == TRUE){
            hEEJs <- hEEJs[hEEJs$seqnames %in%
                                  c(paste("chr", c(1:22, "X", "Y"), sep="")), ]
        }
        hEEJs <- hEEJs[(hEEJs$end - hEEJs$start - 1) >= minLength, ]
        if (!is.null(x=maxLength)){
            hEEJs <- hEEJs[(hEEJs$end - hEEJs$start - 1) <= maxLength, ]
        }
        rownames(x=hEEJs) <- NULL
    }
    ### Inferring of cis-SAGe.
    start_plus <- hEEJs
    start_plus$strand <- "+"
    start_plus <- makeGRangesFromDataFrame(df=start_plus,
                                           keep.extra.columns=TRUE)
    hits <- findOverlaps(query=start_plus, subject=GRs,
                         type="within",
                         select="all",
                         ignore.strand=FALSE)
    start_plus <- start_plus[-queryHits(x=hits), ]
    end_plus <- start_plus
    end(x=start_plus) <- start(x=start_plus)
    start(x=end_plus) <- end(x=end_plus)
    hits <- findOverlaps(query=start_plus, subject=GRs,
                         type="any",
                         select="all",
                         ignore.strand=FALSE)
    start_plus <- start_plus[queryHits(x=hits), ]
    start_plus$gene1 <- GRs[subjectHits(x=hits), 1]
    hits <- findOverlaps(query=end_plus, subject=GRs,
                         type="any",
                         select="all",
                         ignore.strand=FALSE)
    end_plus <- end_plus[queryHits(x=hits), ]
    end_plus$gene2 <- GRs[subjectHits(x=hits), 1]
    start_plus <- start_plus[start_plus$id %in% end_plus$id, ]
    end_plus <- end_plus[end_plus$id %in% start_plus$id, ]
    end(x=start_plus) <- as.numeric(x=gsub(pattern=".+-",
                                           replacement="",
                                           x=start_plus$id))
    start(x=end_plus) <- as.numeric(x=gsub(pattern="c.+:",
                                           replacement="",
                                           x=gsub(pattern="-.+",
                                                  replacement="",
                                                  x=end_plus$id)))
    hits <- findOverlaps(query=start_plus, subject=end_plus,
                         type="equal",
                         select="all",
                         ignore.strand=FALSE)
    start_plus <- start_plus[queryHits(x=hits), ]
    end_plus <- end_plus[subjectHits(x=hits), ]
    plus <- cbind(data.frame(start_plus), data.frame(end_plus))
    plus <- plus[, c(6, 1:2, 16, 4:5, 7, 13, 8:10, 26, 21:23)]
    plus[, c(2, 6, 9, 13)] <- apply(X=plus[, c(2, 6, 9, 13)],
                                    MARGIN=2,
                                    FUN=as.character)
    start_minus <- hEEJs
    start_minus$strand <- "-"
    start_minus <- makeGRangesFromDataFrame(df=start_minus,
                                            keep.extra.columns=TRUE)
    hits <- findOverlaps(query=start_minus, subject=GRs,
                         type="within",
                         select="all",
                         ignore.strand=FALSE)
    start_minus <- start_minus[-queryHits(x=hits), ]
    end_minus <- start_minus
    end(x=start_minus) <- start(x=start_minus)
    start(x=end_minus) <- end(x=end_minus)
    hits <- findOverlaps(query=start_minus, subject=GRs,
                         type="any",
                         select="all",
                         ignore.strand=FALSE)
    start_minus <- start_minus[queryHits(x=hits), ]
    start_minus$gene1 <- GRs[subjectHits(x=hits), 1]
    hits <- findOverlaps(query=end_minus, subject=GRs,
                         type="any",
                         select="all",
                         ignore.strand=FALSE)
    end_minus <- end_minus[queryHits(x=hits), ]
    end_minus$gene2 <- GRs[subjectHits(x=hits), 1]
    start_minus <- start_minus[start_minus$id %in% end_minus$id, ]
    end_minus <- end_minus[end_minus$id %in% start_minus$id, ]
    end(x=start_minus) <- as.numeric(x=gsub(pattern=".+-",
                                            replacement="",
                                            x=start_minus$id))
    start(x=end_minus) <- as.numeric(x=gsub(pattern="c.+:",
                                            replacement="",
                                            x=gsub(pattern="-.+",
                                                   replacement="",
                                                   x=end_minus$id)))
    hits <- findOverlaps(query=start_minus, subject=end_minus,
                         type="equal",
                         select="all",
                         ignore.strand=FALSE)
    start_minus <- start_minus[queryHits(x=hits), ]
    end_minus <- end_minus[subjectHits(x=hits), ]
    minus <- cbind(data.frame(start_minus), data.frame(end_minus))
    minus <- minus[, c(6, 1:2, 16, 4:5, 7, 13, 8:10, 26, 21:23)]
    minus[, c(2, 6, 9, 13)] <- apply(X=minus[, c(2, 6, 9, 13)],
                                     MARGIN=2,
                                     FUN=as.character)
    cisSAGe <- rbind(plus, minus)
    cisSAGe <- cisSAGe[order(x=cisSAGe$id), ]
    rownames(x=cisSAGe) <- NULL
    ### Annotation of cis-SAGe using pairs of paralogous genes.
    if (!is.null(x=parlgs)){
        PGs <- read.table(file=paste(workDir, parlgs, sep="/"),
                          header=TRUE,
                          sep="\t",
                          quote="\"",
                          as.is=TRUE)
        cisSAGe$paralogous <- "no"
        cisSAGe[paste(cisSAGe$gene1.gene_id, cisSAGe$gene2.gene_id) %in%
                paste(PGs$gene1, PGs$gene2), ]$paralogous <- "yes"
    }
    ### Returning a final object of class data frame.
    return(cisSAGe)
}
