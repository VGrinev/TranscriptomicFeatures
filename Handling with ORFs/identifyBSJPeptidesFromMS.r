#' Identification of peptides spanning the amino acid sequence encoded by the
#'     back-splicing junction site region of human circular RNA molecules
#' @description Identify the BSJ encoded peptides of circRNAs.
#' @param circRNAs character string giving the name of file with the sequences
#'     of circRNAs. This file must be in the working directory. Allowed file
#'     formats are "fasta" or "fa".
#' @param circORFs character string giving the name of file with the sequences
#'     of circRNA open reading frames. This file must be in the working
#'     directory. Allowed file formats are "fasta" or "fa".
#' @param proteins character string giving the name of file with the sequences
#'     of circRNAs encoded proteins. This file must be in the working directory.
#'     Allowed file formats are "fasta" or "fa".
#' @param MSproteins character string giving the name of file with the
#'     annotations of mass spectrometer detected circRNAs encoded proteins. In
#'     the current version of function, only thirteen-field CSV files generated
#'     by PEAKS(R) Studio 11 software are accepted.
#' @param MSpeptides character string giving the name of file with the
#'     annotations of mass spectrometer detected circRNAs encoded peptides. In
#'     the current version of function, only seventeen-field CSV files generated
#'     by PEAKS(R) Studio 11 software are accepted.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return The data frame containing the following fields:
#'     i) protein_id       - protein ID;
#'     ii) pr_length       - length of protein;
#'     iii) junSites       - location of BSJ site(-s) in protein;
#'     iv) peptides        - character string giving the sequence(-s) of MS
#'                           detected peptide(-s);
#'     v) BSJpeptides      - character string giving the sequence(-s) of MS
#'                           detected BSJ peptide(-s);
#'     vi) neg_log10_p     - negative value of log10-transformed adjusted
#'                           p-value;
#'     vii) coverage1      - protein sequence coverage (%) by all MS detected
#'                           peptides;
#'     viii) coverage2     - protein sequence coverage (%) by peptides detected
#'                           in given sample;
#'     ix) area            - the area under the curve of the peptide feature
#'                           found at the same m/z and retention time as the
#'                           MS/MS scan;
#'     x) peptides_total   - total number of MS detected peptides supporting
#'                           given protein;
#'     xi) peptides_unique - number of unique MS detected peptides supporting
#'                           given protein;
#'     xii) uniqueness     - uniqueness score.
#' @authors Vasily V. Grinev, Dai Xiaoxuan
#' @examples
#' pepBSJ <- identifyBSJPeptidesFromMS(circRNAs="TransCirc, circRNAs, sequences.fa",
#'                                     circORFs="TransCirc, circRNAs, ORF sequences.fa",
#'                                     proteins="TransCirc, circRNAs, ORF translations.fasta",
#'                                     MSproteins="proteins.csv",
#'                                     MSpeptides="peptides.csv",
#'                                     workDir="D:/Vasily Grinev")
#' @export
#' Last updated: May 8, 2023.

identifyBSJPeptidesFromMS <- function(circRNAs,
                                      circORFs,
                                      proteins,
                                      MSproteins,
                                      MSpeptides,
                                      workDir=NULL){
    ### Loading of the required package.
    #   This code was successfully tested with the package Biostrings v.2.60.2.
    suppressMessages(expr=library(package=Biostrings))
    ### Loading of the circRNA sequences as an object of class DNAStringSet.
    ##  Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=circRNAs)
    ##  Validation of file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Full path to the file.
    circRNA_seq <- paste(workDir, circRNAs, sep="/")
    ##  Loading of the circRNA sequences.
    circRNA_seq <- readDNAStringSet(filepath=circRNA_seq)
    ### Loading of the sequences of circRNA ORFs as an object of class
    #   DNAStringSet.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=circORFs)
    ##  Validation of file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Full path to the file.
    orf_seq <- paste(workDir, circORFs, sep="/")
    ##  Loading of the sequences of circRNA ORFs.
    orf_seq <- readDNAStringSet(filepath=orf_seq)
    orf_seq <- orf_seq[width(x=orf_seq) < 20000]
    orf_seq <- orf_seq[width(x=orf_seq) >= 153]
    ### Loading of the MS detected peptide annotations as an object of class
    #   data frame.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=MSpeptides)
    ##  Validation of the file format.
    if (!frt %in% c("csv")){
        stop("Invalid file format")
    }
    ##  Loading of the peptide annotations.
    pe <- read.csv(file=paste(workDir, MSpeptides, sep = "/"),
                   comment.char="!")[, c(1, 14)]
    colnames(x=pe) <- c("peptide", "protein_id")
    ##  Clearing of peptide sequences.
    pe$peptide <- gsub(pattern="[(]sub ", replacement="", x=pe$peptide)
    pe$peptide <- gsub(pattern="[-+.0-9()]", replacement="", x=pe$peptide)
    ##  Filtering against empty protein records.
    pe <- pe[!nchar(x=pe$protein_id) == 0, ]
    ##  Development of final data frame object.
    pe$protein_id <- strsplit(x=pe$protein_id, split="[:]")
    res <- list()
    for (i in 1:nrow(x=pe)){
        res[[i]] <- cbind(pe$protein_id[i][[1]], pe$peptide[i])
    }
    pe <- do.call(what=rbind, args=res)
    pe <- aggregate(x=pe[, 2],
                    by=list(pe[, 1]),
                    FUN=function(y){y=paste(sort(x=unique(x=y)),
                                            collapse=",")})
    colnames(x=pe) <- c("protein_id", "peptides")
    pe$protein_id <- gsub(pattern=" ", replacement=",", x=pe$protein_id)
    ### Loading of the MS detected protein annotations as an object of class
    #   data frame.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=MSproteins)
    ##  Validation of the file format.
    if (!frt %in% c("csv")){
        stop("Invalid file format")
    }
    ##  Loading of the protein annotations.
    pr <- read.csv(file=paste(workDir, MSproteins, sep = "/"),
                   comment.char="!")[, 3:10]
    colnames(x=pr) <- c("protein_id",
                        "neg_log10_p",
                        "coverage1", "coverage2", "area",
                        "peptides_total", "peptides_unique",
                        "uniqueness")
    pr$protein_id <- gsub(pattern=" ", replacement=",", x=pr$protein_id)
    ##  Filtration of the primary data frame.
    pr <- pr[pr[, 2] >= 40, ]
    pr <- pr[pr[, 6] >= 2, ]
    pr <- pr[pr[, 7] >= 1, ]
    ##  Annotation of proteins with peptides.
    rownames(x=pe) <- pe[, 1]
    pe <- as.matrix(pe)
    pr$peptides <- ""
    pr <- as.matrix(pr)
    idx <- pr[, 1] %in% pe[, 1]
    pr[idx, 9] <- pe[pr[idx, 1], 2]
    pr <- data.frame(pr)
    pr[, 2:8] <- apply(X=pr[, 2:8], MARGIN=2, FUN=as.numeric)
    ### Loading of the MS detected proteins as an object of class DNAStringSet.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=proteins)
    ##  Validation of the file format.
    if (!frt %in% c("fasta", "fa")){
        stop("Invalid file format")
    }
    ##  Loading of the protein sequences.
    pr_seq <- readAAStringSet(filepath=paste(workDir, proteins, sep = "/"))
    pr_seq <- pr_seq[width(x=pr_seq) >= 50]
    pr <- pr[pr$protein_id %in% names(x=pr_seq), ]
    pr_seq <- pr_seq[names(x=pr_seq) %in% pr$protein_id]
    ### Mapping of back-splicing junction sites in circRNA ORFs.
    orf_seq <- orf_seq[names(x=orf_seq) %in%
                       unlist(x=strsplit(x=pr$protein_id, split=","))]
    circRNA_seq <- circRNA_seq[names(x=circRNA_seq) %in%
                               gsub(pattern="_ORF.+",
                                    replacement="",
                                    x=names(x=orf_seq))]
    BSJsites <- list()
    for (i in 1:length(x=orf_seq)){
        ORF <- as.character(x=orf_seq[[i]])
        name <- sub(pattern="_ORF.+", replacement="", x=names(x=orf_seq[i]))
        circRNA <- as.character(x=circRNA_seq[names(circRNA_seq) == name][[1]])
        circRNA <- paste(rep(x=circRNA, times=4), collapse="")
        coord <- vmatchPattern(pattern=ORF, subject=circRNA)[[1]]
        coord <- coord[start(x=coord) == min(x=start(x=coord)), ]
        length_circRNA <- width(x=circRNA_seq[names(circRNA_seq) == name])
        orf_end <- length_circRNA - ((((end(x=coord) - length_circRNA) %/%
                   length_circRNA) + 1) * length_circRNA - (end(x=coord) -
                   length_circRNA))
        orf_end <- paste(paste(((end(x=coord) - length_circRNA) %/%
                                 length_circRNA) + 1, "r", sep=""),
                         orf_end, sep="+")
        BSJsites[[i]] <- c(name,
                           names(x=orf_seq[i]),
                           length_circRNA,
                           start(x=coord),
                           orf_end,
                           nchar(x=ORF))
    }
    BSJsites <- data.frame(do.call(what=rbind, args=BSJsites))
    colnames(x=BSJsites) <- c("transcript_id",
                              "orf_id",
                              "transcript_length",
                              "orf_start",
                              "orf_end",
                              "orf_length")
    BSJsites$transcript_length <- as.numeric(x=BSJsites$transcript_length)
    BSJsites$orf_start <- as.numeric(x=BSJsites$orf_start)
    BSJsites$orf_length <- as.numeric(x=BSJsites$orf_length)
    BSJsites <- BSJsites[as.numeric(x=gsub(pattern="r.+",
                                    replacement="",
                                    x=BSJsites$orf_end)) > 0, ]
    junSites <- list()
    for (i in 1:nrow(x=BSJsites)){
        sites <- c(BSJsites[i, 3] - BSJsites[i, 4] + 1,
                   if (as.numeric(x=gsub(pattern="r.+",
                                         replacement="",
                                         x=BSJsites[i, 5])) > 1){
                       1:(as.numeric(x=gsub(pattern="r.+",
                                            replacement="",
                                            x=BSJsites[i, 5])) - 1) *
                       BSJsites[i, 3] + (BSJsites[i, 3] - BSJsites[i, 4] + 1)})
        sites <- round(x=sites/3, digits=0)
        junSites[[i]] <- IRanges(start=sites, end=sites)
    }                      
    names(x=junSites) <- BSJsites$orf_id
    ### Development of peptide map.
    peptideMap <- list()
    for (i in 1:nrow(x=pr)){
        peptides <- AAStringSet(x=strsplit(x=pr$peptides[i], split=",")[[1]])
        protein <- pr_seq[names(x=pr_seq) == pr$protein_id[i]][[1]]
        map <- sort(x=unlist(x=matchPDict(pdict=peptides, subject=protein)))
        peptideMap[[i]] <- map
    }
    names(x=peptideMap) <- gsub(pattern=",.+", replacement="", x=pr$protein_id)
    ### Identification of peptides overlapping BSJ site(-s).
    res <- list()
    for (i in 1:length(x=peptideMap)){
        if (names(x=peptideMap[i]) %in% names(x=junSites)){
            subject <- junSites[names(x=junSites) %in%
                                names(x=peptideMap[i])][[1]]
            res[[i]] <- findOverlaps(query=peptideMap[[i]], subject=subject)
        }
    }
    names(x=res) <- names(x=peptideMap)
    res <- res[lengths(x=res) > 0]
    overlaps <- list()
    for (i in 1:length(x=res)){
        query <- peptideMap[names(x=peptideMap) == names(x=res)[i]][[1]]
        query <- query[queryHits(x=res[[i]])]
        subject <- junSites[names(x=junSites) == names(x=res)[i]][[1]]
        subject <- subject[subjectHits(x=res[[i]])]
        if (max(x=end(x=query) - end(x=subject)) > 2){
            overlaps[[i]] <- end(x=query) - end(x=subject)
        }else{
            overlaps[[i]] <- ""
        }
    }
    names(x=overlaps) <- names(x=res)
    overlaps <- overlaps[!overlaps == ""]
    res <- res[names(x=res) %in% names(x=overlaps)]
    overlaps <- list()
    for (i in 1:length(x=res)){
        query <- peptideMap[names(x=peptideMap) == names(x=res)[i]][[1]]
        query <- query[queryHits(x=res[[i]])]
        subject <- junSites[names(x=junSites) == names(x=res)[i]][[1]]
        subject <- subject[subjectHits(x=res[[i]])]
        if (max(x=start(x=subject) - start(x=query) + 1) > 2){
            overlaps[[i]] <- start(x=subject) - start(x=query) + 1
        }else{
            overlaps[[i]] <- ""
        }
    }
    names(x=overlaps) <- names(x=res)
    overlaps <- overlaps[!overlaps == ""]
    res <- res[names(x=res) %in% names(x=overlaps)]
    BSJpeptides <- list()
    for (i in 1:length(x=res)){
        peptides <- pr[gsub(pattern=",.+", replacement="", x=pr$protein_id) ==
                       names(x=res)[i], ]$peptides
        peptides <- AAStringSet(x=strsplit(x=peptides, split=",")[[1]])
        protein <- pr_seq[gsub(pattern=",.+",
                               replacement="",
                               x=names(x=pr_seq)) == names(x=res)[i]][[1]]
        peptides <- peptides[lengths(x=matchPDict(pdict=peptides,
                                                  subject=protein)) > 0]
        query <- matchPDict(pdict=peptides, subject=protein)
        subject <- junSites[names(x=junSites) %in%
                            names(x=res)[i]][[1]]
        overlap <- lapply(X=query, FUN=function(y){findOverlaps(query=y, subject=subject)})
        peptides <- peptides[lengths(x=overlap) > 0]
        peptide <- list()
        for (j in 1:length(x=peptides)){
            query <- unlist(x=matchPDict(pdict=peptides[j], subject=protein))
            overlap <- findOverlaps(query=query, subject=subject)
            query <- query[queryHits(x=overlap)]
            subj <- subject[subjectHits(x=overlap)]
            if (sum(x=as.numeric(x=(((start(x=subj) - start(x=query) + 1) > 2) &
                                    ((end(x=query) - end(x=subj)) > 2)))) > 0){
                peptide[[j]] <- as.character(x=peptides[j][[1]])
            }
        }
        peptide <- do.call(what=c, args=peptide)
        peptide <- paste(sort(x=unique(x=peptide)), collapse=",")
        BSJpeptides[[i]] <- c(names(x=res)[i], peptide)
    }
    BSJpeptides <- do.call(what=rbind, args=BSJpeptides)
    ### Creating the final object.
    junSites <- unlist(x=lapply(X=junSites,
                                FUN=function(y){paste(start(x=y),
                                                      collapse=",")}))
    junSites <- cbind(attr(x=junSites, which="names"), junSites)
    rownames(x=junSites) <- NULL
    junSites <- junSites[junSites[, 1] %in% gsub(pattern=",.+",
                                                 replacement="",
                                                 x=pr$protein_id), ]
    junSites <- junSites[!duplicated(x=junSites), ]
    pr$name <- gsub(pattern=",.+", replacement="", x=pr$protein_id)
    pr <- pr[, c(10, 1:9)]
    pr$pr_length <- 0
    pr$junSites <- ""
    pr$BSJpeptides <- ""
    pr <- as.matrix(pr)
    pr_length <- cbind(gsub(pattern=",.+", replacement="", x=names(x=pr_seq)),
                       width(x=pr_seq))
    rownames(x=pr_length) <- pr_length[, 1]
    pr_length <- as.matrix(pr_length)
    idx <- pr[, 1] %in% pr_length[, 1]
    pr[idx, 11] <- pr_length[pr[idx, 1], 2]
    rownames(x=junSites) <- junSites[, 1]
    idx <- pr[, 1] %in% junSites[, 1]
    pr[idx, 12] <- junSites[pr[idx, 1], 2]
    rownames(x=BSJpeptides) <- BSJpeptides[, 1]
    idx <- pr[, 1] %in% BSJpeptides[, 1]
    pr[idx, 13] <- BSJpeptides[pr[idx, 1], 2]
    pr <- data.frame(pr)
    pr <- pr[, c(2, 11:12, 10, 13, 3:9)]
    pr[, c(2, 6:12)] <- apply(X=pr[, c(2, 6:12)], MARGIN=2, FUN=as.numeric)
    pr <- pr[!pr$junSites == "", ]
    ### Returning the final object of class data frame.
    return(pr)
}
