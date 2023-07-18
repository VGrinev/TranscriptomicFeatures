#' Identify of region specific peptides from mass spectrometer data
#' @description Identification of the region specific peptides from mass
#'     spectrometer data.
#' @param x character string giving the name of first input file with proteins
#'     (and related peptides) detected by mass spectrometer. It is usually
#'     output file of function parseDIANNoutput(). This file must be in the
#'     working directory. Allowed file format is "txt".
#' @param y character string giving the name of second input file with
#'     sequences of proteins detected by mass spectrometer. This file must be
#'     in the working directory. Allowed file formats are "fasta" or "fa".
#' @param z character string giving the name of third input file with
#'     protein coordinates of specific regions. It is usually output file of
#'     function contigKmers(). This file must be in the working directory.
#'     Allowed file format is "txt".
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return The list containing sequences (in the form of objects of class
#'     AAStringSet) of validated proteins, validated specific regions, region
#'     specific piptides and object of class data frame containing the
#'     following four fields:
#'     i) protein_id     - protein ID;
#'     ii) region_coords - protein coordinates of specific region;
#'     iii) region_seqs  - sequence of specific region;
#'     iv) region_pepts  - region specific piptide(-s).
#' @author Vasily V. Grinev
#' @examples
#' res <- identifyRSPeptidesFromMS(x="Kasumi-1 proteome, all proteins detected at FDR 5%.txt",
#'                                 y="Kasumi-1 proteome, June 23, 2023.fasta",
#'                                 z="Kasumi-1 proteome, June 23, 2023, all specific 8-mers, contigs.txt",
#'                                 workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 10, 2023.

identifyRSPeptidesFromMS <- function(x,
                                     y,
                                     z,
                                     workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with packages Biostrings v.2.60.2 and
    #   GenomicRanges v.1.40.0.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=GenomicRanges))
    ### Loading of the mass spectrometer detected proteins as an object of
    #   class data frame.
    ##  Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=x)
    ##  Validation of the file format.
    if (!frt %in% c("txt")){
        stop("Invalid file format")
    }
    ##  Loading of the proteins.
    pr_peptides <- read.table(file=paste(workDir, x, sep = "/"),
                              sep="\t",
                              header=TRUE,
                              quote="\"",
                              as.is=TRUE)[, c(1, 3)]
    ### Loading of the mass spectrometer detected proteins as an object of
    #   class AAStringSet.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=y)
    ##  Validation of file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Loading of the proteins sequences.
    pr_seq <- readAAStringSet(filepath=paste(workDir, y, sep = "/"))
    names(x=pr_seq) <- substring(text=names(x=pr_seq), first=4, last=16)
    pr_seq <- pr_seq[names(x=pr_seq) %in% pr_peptides$protein_id, ]
    ### Loading of the specific protein regions as an object of class GRanges.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=z)
    ##  Validation of the file format.
    if (!frt %in% c("txt")){
        stop("Invalid file format")
    }
    ##  Loading of the regions.
    regions <- read.table(file=paste(workDir, z, sep = "/"),
                          header=TRUE,
                          quote="\"",
                          as.is=TRUE)
    regions <- makeGRangesFromDataFrame(df=regions, keep.extra.columns=TRUE)
    regions <- regions[seqnames(x=regions) %in% pr_peptides$protein_id, ]
    ### Subsetting of all data sets.
    IDs <- pr_peptides$protein_id[pr_peptides$protein_id %in%
                                  seqnames(x=regions)]
    ### Extraction of the mass spectrometer peptides as an object of class
    #   AAStringSet.
    peps_old <- strsplit(x=pr_peptides[!pr_peptides$protein_id %in%
                                       IDs, ]$support_peptides_seq,
                         split=";")
    names(x=peps_old) <- pr_peptides[!pr_peptides$protein_id %in%
                                     IDs, ]$protein_id
    peps_old <- AAStringSet(x=unlist(x=peps_old))
    names(x=peps_old) <- substring(text=names(x=peps_old), first=1, last=13)
    peps <- strsplit(x=pr_peptides[pr_peptides$protein_id %in%
                                   IDs, ]$support_peptides_seq,
                     split=";")
    names(x=peps) <- pr_peptides[pr_peptides$protein_id %in%
                                 IDs, ]$protein_id
    peps <- AAStringSet(x=unlist(x=peps))
    names(x=peps) <- substring(text=names(x=peps), first=1, last=13)
    pr_seq <- pr_seq[names(x=pr_seq) %in% names(x=peps), ]
    ### Development of peptide map.
    pMap <- list()
    for (i in 1:length(x=pr_seq)){
        pdict <- peps[names(x=peps) == names(x=pr_seq)[i], ]
        pMap[[i]] <- GRanges(seqnames=names(x=pr_seq[i]),
                             ranges=unlist(x=matchPDict(pdict=pdict,
                                                        subject=pr_seq[[i]])))
    }
    pMap <- sort(x=suppressWarnings(expr=unlist(x=as(object=pMap,
                                                     Class="GRangesList"))))
    names(x=pMap) <- NULL
    ### Development of complex report.
    start(regions) <- start(regions) + 1
    end(regions) <- end(regions) - 1
    hits <- findOverlaps(query=pMap,
                         subject=regions,
                         minoverlap=7L,
                         type="any")
    regions <- regions[sort(x=unique(x=subjectHits(x=hits))), ]
    pMap <- pMap[sort(x=unique(x=queryHits(x=hits))), ]
    start(regions) <- start(regions) - 1
    end(regions) <- end(regions) + 1
    pr_seqs <- pr_seq[names(x=pr_seq) %in% seqnames(x=regions), ]
    sub_seqs <- list()
    for (i in 1:length(x=pr_seqs)){
        region <- regions[seqnames(x=regions) %in% names(x=pr_seqs[i]), ]
        sub_seq <- AAStringSet(x=substring(text=pr_seqs[i],
                                           first=start(x=region),
                                           last=end(x=region)))
        names(x=sub_seq) <- rep(x=names(x=pr_seqs[i]), length(x=(sub_seq))) 
        sub_seqs[[i]] <- sub_seq
    }
    pr_frags <- do.call(what=c, args=sub_seqs)
    sub_seqs <- list()
    for (i in 1:length(x=pr_seqs)){
        region <- pMap[seqnames(x=pMap) %in% names(x=pr_seqs[i]), ]
        sub_seq <- AAStringSet(x=substring(text=pr_seqs[i],
                                           first=start(x=region),
                                           last=end(x=region)))
        names(x=sub_seq) <- rep(x=names(x=pr_seqs[i]), length(x=(sub_seq)))
        sub_seqs[[i]] <- sub_seq
    }
    pr_pepts <- do.call(what=c, args=sub_seqs)
    anno <- regions
    anno$pr_frags <- pr_frags
    anno$region <- paste(start(x=anno), end(x=anno), sep="-")
    anno <- data.frame(anno)[, c(1, 6:7)]
    anno$seqnames <- as.character(x=anno$seqnames)
    f <- function(y){paste(unique(x=y), collapse=";")}
    anno <- aggregate(x=anno[, -1], by=list(anno[, 1]), FUN=f)
    anno <- cbind(anno[, c(1, 3:2)],
                  aggregate(x=as.character(x=pr_pepts),
                  by=list(names(x=as.character(x=pr_pepts))),
                  FUN=f)[, 2])
    colnames(x=anno) <- c("protein_id",
                          "region_coords",
                          "region_seqs",
                          "region_pepts")
    ### Returning the final object of class data frame.
    return(list(proteins=pr_seqs,
                fragments=pr_frags,
                peptides=pr_pepts,
                annotations=anno))
}

setwd(dir="D:/Vasily Grinev")
res <- identifyRSPeptidesFromMS(x="Kasumi-1 proteome, all proteins detected at FDR 5%.txt",
                                y="Kasumi-1 proteome, June 23, 2023.fasta",
                                z="Kasumi-1 proteome, June 23, 2023, all specific 8-mers, contigs.txt",
                                workDir="D:/Vasily Grinev")
suppressMessages(expr=library(package=Biostrings))
writeXStringSet(x=res$fragments, filepath="fragments.fasta")
