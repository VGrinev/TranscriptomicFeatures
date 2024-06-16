#' Extract the nucleotide sequence(-s) flanking breakpoint(-s)
#' @description Extraction of the nucleotide sequence(-s) flanking breakpoint(-s).
#' @param vcfFiles character vector giving the name(-s) of VCF file(-s) to be
#'     processed.
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param canChr a logical argument specifying that the records from
#'     non-canonical chromosomes/scaffolds should be removed. Default value is
#'     FALSE.
#' @param thr an integer argument. It is a threshold for minimal coverage
#'     (sequencing depth) of breakpoint (in raw reads). Default value is 5.
#' @param group1 character vector of samples for group 1. Default value is NULL.
#' @param group2 character vector of samples for group 2. Default value is NULL.
#' @param flankL an integer argument giving the length for flanking sequences.
#'     Default value is 100.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return DNAStringSet object of breakpoint flanking sequence(-s).
#' @author Vasily V. Grinev, Katerina V. Huzava
#' @examples
#' res <- getFlankSeqBPs(vcfFiles=NULL,
#'                       genome="BSgenome.Hsapiens.UCSC.hg38",
#'                       canChr=FALSE,
#'                       thr=5,
#'                       group1=NULL,
#'                       group2=NULL,
#'                       flankL=100,
#'                       workDir="D:/Files_VCF")
#' @export
#' Last updated: June 16, 2024.

getFlankSeqBPs <- function(vcfFiles=NULL,
                           genome="BSgenome.Hsapiens.UCSC.hg38",
                           canChr=FALSE,
                           thr=5,
                           group1=NULL,
                           group2=NULL,
                           flankL=100,
                           workDir=NULL){
    ### Loading of the required packages.
    #   This code was successfully tested with packages
    #   StructuralVariantAnnotation v.1.18.0
    #   and BSgenome.Hsapiens.UCSC.hg38 v.1.4.5.
    suppressMessages(expr=library(package=StructuralVariantAnnotation))
    suppressMessages(expr=library(package=BSgenome.Hsapiens.UCSC.hg38))
    ### Development a list of VCF files.
    ##  Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ##  Development the list of VCF files.
    if (is.null(x=vcfFiles)){
        vcf <- list.files(path=workDir, pattern="[.]breakpoints.vcf$")
    vcf <- paste(workDir, vcf, sep="/") 
    }else{
        vcf <- paste(workDir, vcfFiles, sep="/")
    }
    ### Loading of the breakpoints data as an object of class list of GRanges.
    GRs <- list()
    IDs <- list()
    for (i in 1:length(x=vcf)){
        bnd <- readVcf(file=vcf[i])
        gr <- breakpointRanges(x=bnd)[, c(2, 3, 6:9)]
        gr$sr <- info(x=bnd)$SR
        if (canChr == TRUE){
            gr <- gr[seqnames(x=gr) %in%
                                    c(paste("chr", c(1:22, "X", "Y"), sep=""))]
        }
        l <- grep(pattern="^[A-Z]", x=gr$ALT)
        r <- grep(pattern="[A-Z]$", x=gr$ALT)
        gr$id <- NA
        gr$id[l] <- paste0(paste(seqnames(x=gr)[l],
                                 start(x=gr)[l],
                                 sep=":"),
                           gsub(pattern="^[A-Z]", replacement="", x=gr$ALT[l]))
        gr$id[r] <- paste0(gsub(pattern="[A-Z]$", replacement="", x=gr$ALT[r]),
                           paste(seqnames(x=gr)[r],
                                 start(x=gr)[r],
                                 sep=":"))
        d <- names(x=gr[duplicated(gr$id), ])
        gr[d, ]$id <- paste0(gr[d, ]$id, ".1")
        GRs[[i]] <- gr
        names(x=GRs)[i] <- sub(pattern=".breakpoints.vcf",
                               replacement="",
                               x=basename(path=vcf[i]))
        IDs[[i]] <- cbind(names(x=GRs)[[i]], gr[gr$sr >= thr]$id)
    }
    ### Developing a list of unique breakpoint identifiers.
    if (!is.null(x=group1)){
        IDs <- do.call(what=rbind, args=IDs)
        IDs_group1 <- table(IDs[IDs[, 1] %in% group1, 2])
        IDs_group1 <- IDs_group1[IDs_group1 == length(x=group1)]
        IDs_group2 <- table(IDs[IDs[, 1] %in% group2, 2])
        IDs_group2 <- IDs_group2[IDs_group2 == length(x=group2)]
        IDs <- sort(x=unique(x=c(attr(IDs_group1, "dimnames")[[1]],
                                 attr(IDs_group2, "dimnames")[[1]])))
    }else{
        id <- IDs[[1]]
        for (i in 2:length(x=IDs)){
            id <- id[id %in% IDs[[i]]]
        }
        IDs <- sort(x=unique(x=id))
    }
    ### Extraction of the breakpoint sequence(-s).
    SEQs <- list()
    for (i in 1:length(x=GRs)){
        gr <- GRs[[i]]
        gr <- gr[gr$id %in% IDs, ]
        seqs <- extractBreakpointSequence(gr=gr,
                                          ref=get(x=genome),
                                          anchoredBases=flankL,
                                          remoteBases=flankL)
        seqs <- DNAStringSet(x=seqs)
        names(x=seqs) <- gr$id
        SEQs[[i]] <- seqs
    }
    SEQs <- do.call(what=c, args=SEQs)
    SEQs <- aggregate(x=names(x=SEQs),
                      by=list(as.character(x=SEQs)),
                      FUN=function(y){paste(sort(x=unique(x=y)),
                                            collapse="; ")})
    SEQs_seq <- DNAStringSet(x=SEQs$Group.1)
    names(x=SEQs_seq) <- SEQs$x
    SEQs <- SEQs_seq
    SEQs_seq <- NULL
    ### Returning a final object of class DNAStringSet.
    return(SEQs)
}
