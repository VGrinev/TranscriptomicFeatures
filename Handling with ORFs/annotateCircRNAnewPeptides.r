#' Annotate new peptides encoded by circRNAs
#' @description Annotation of the new peptides encoded by circRNAs.
#' @param newPepts character string giving the name of input file with
#'     sequences of the new peptides encoded by circRNAs. This file must be in
#'     the working directory. Valid file format can be either "fasta" or "fa".
#' @param msProts character string giving the name of input file with
#'     sequences of the mass spectrometry (MS) approved circRNAs encoded
#'     proteins. This file must be in the working directory. Valid file format
#'     can be either "fasta" or "fa".
#' @param msPepts character string giving the name of file with the annotations
#'     of MS detected circRNAs encoded peptides. This file must be in the
#'     working directory. Valid file format is "txt". It is typically output
#'     file of function identifyBSJPeptidesFromMS().
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return The data frame containing the following fields:
#'     i) gene_symbol                    - gene symbol;
#'     ii) protein_id                    - protein ID;
#'     iii) protein_length               - length of protein;
#'     iv) support_peptides              - character string giving the
#'                                         sequence(-s) of MS detected
#'                                         peptide(-s);
#'     v) neg_log10_p                    - negative value of log10-transformed
#'                                         adjusted p-value;
#'     vi) coverage1                     - protein sequence coverage (%) by all
#'                                         MS detected peptides;
#'     vii) coverage2                    - protein sequence coverage (%) by
#'                                         peptides detected in given sample;
#'     viii) area                        - the area under the curve of the
#'                                         peptide feature found at the same
#'                                         m/z and retention time as the MS/MS
#'                                         scan;
#'     ix) support_peptides_total        - total number of MS detected peptides
#'                                         supporting given protein;
#'     x) support_peptides_unique        - number of unique MS detected
#'                                         peptides supporting given protein;
#'     xi) uniqueness                    - uniqueness score.
#'     xii) BSJsite(-s)                  - location of BSJ site(-s) in protein;
#'     xiii) BSJpeptide(-s)              - character string giving the
#'                                         sequence(-s) of MS detected BSJ
#'                                         peptide(-s);
#'     xiv) new_peptide                  - the sequence of new peptide encoded
#'                                         by circRNA;
#'     xv) new_peptide_location          - protein location of new peptide
#'                                         encoded by circRNA;
#'     xvi) new_peptide_support          - MS detected peptide(-s) supporting
#'                                         of new peptide encoded by circRNA;
#'     xvii) new_peptide_support_n       - number of MS detected peptide(-s)
#'                                         supporting of new peptide encoded
#'                                         by circRNA;
#'     xviii) new_peptide_support_length - length (in a.a.) of MS detected
#'                                         peptide(-s) supporting of new
#'                                         peptide encoded by circRNA.
#' @author Vasily V. Grinev
#' @examples
#' newPepts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, new peptides, sequences.fasta"
#' msProts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, sequences.fasta"
#' msPepts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, MS peptides, annotations.txt"
#' workDir="D:/Vasily Grinev"
#' res <- annotateCircRNAnewPeptides(newPepts=newPepts, msProts=msProts, msPepts=msPepts, workDir=workDir)
#' @export
#' Last updated: May 20, 2023.

annotateCircRNAnewPeptides <- function(newPepts,
                                       msProts,
                                       msPepts,
                                       workDir=NULL){
    ### Loading of the required packages.
    #   This code was successfully tested with packages Biostrings v.2.60.2 and
    #   GenomicRanges v.1.46.1.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=GenomicRanges))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of the new peptides encoded by circRNAs as an AAStringSet
    #   object.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=newPepts)
    ##  Validation of the file format.
    if (!frt %in% c("fasta", "fa")){
        stop("Invalid file format")
    }
    ##  Loading of the new peptides encoded by circRNAs.
    new_pepts <- readAAStringSet(filepath=paste(workDir, newPepts, sep = "/"))
    ### Loading of the circRNAs encoded MS detected proteins as an AAStringSet
    #   object.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=msProts)
    ##  Validation of the file format.
    if (!frt %in% c("fasta", "fa")){
        stop("Invalid file format")
    }
    ##  Loading of the circRNAs encoded MS detected proteins.
    ms_prots <- readAAStringSet(filepath=paste(workDir, msProts, sep = "/"))
    ms_prots <- ms_prots[names(x=ms_prots) %in% names(x=new_pepts), ]
    ### Loading of the MS detected peptides as an object of class data frame.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=msPepts)
    ##  Validation of the file format.
    if (!frt %in% c("txt")){
        stop("Invalid file format")
    }
    ##  Loading of the peptides.
    ms_pepts <- read.table(file=paste(workDir, msPepts, sep = "/"),
                           sep="\t",
                           header=TRUE,
                           quote="\"",
                           as.is=TRUE)
    ms_pepts <- ms_pepts[ms_pepts$protein_id %in% names(x=new_pepts), ]
    ms_pepts$new_peptide <- ""
    ms_pepts$new_peptide_location <- ""
    ms_pepts$new_peptide_support <- ""
    ms_pepts$new_peptide_support_n <- 0
    ms_pepts$new_peptide_support_length <- ""
    ### Annotation of the circRNAs encoding new peptides.
    for (i in 1:nrow(x=ms_pepts)){
        protein <- ms_prots[names(x=ms_prots) %in% ms_pepts$protein_id[i]][[1]]
        pepts <- AAStringSet(x=strsplit(x=ms_pepts[i, ]$peptides,
                                        split=",")[[1]])
        names(x=pepts) <- pepts
        pepts <- GRanges(seqnames=ms_pepts$protein_id[i],
                         ranges=unlist(x=matchPDict(pdict=pepts,
                         subject=protein)))
        fr <- new_pepts[names(x=new_pepts) %in% ms_pepts$protein_id[i]][[1]]
        subject <- matchPattern(pattern=fr, subject=protein)
        subject <- GRanges(seqnames=ms_pepts$protein_id[i],
                           ranges=IRanges(start=start(x=subject),
                                          end=end(x=subject)))
        pepts <- subsetByOverlaps(x=pepts,
                                  ranges=subject,
                                  minoverlap=3L,
                                  type="any")
        if (length(x=pepts) > 0){
            ms_pepts$new_peptide[i] <- as.character(x=fr)
            ms_pepts$new_peptide_location[i] <-
              paste(sort(x=paste(start(x=subject), end(x=subject), sep="...")),
                                   collapse=",")
            ms_pepts$new_peptide_support[i] <-
                          paste(sort(x=unique(x=names(x=pepts))), collapse=",")
            ms_pepts$new_peptide_support_n[i] <-
                                     length(x=sort(x=unique(x=names(x=pepts))))
            ms_pepts$new_peptide_support_length[i] <-
                 paste(nchar(x=sort(x=unique(x=names(x=pepts)))), collapse=",")
        }
    }
    ms_pepts <- ms_pepts[!ms_pepts$new_peptide == "", ]
    rownames(x=ms_pepts) = NULL
    ms_pepts$gene_symbol <- gsub(pattern="_.+",
                                 replacement="",
                                 x=gsub(pattern="TC-hsa-",
                                        replacement="",
                                        x=ms_pepts$protein_id))
    ms_pepts <- ms_pepts[, c(18, 1:2, 4, 6:12, 3, 5, 13:17)]
    colnames(x=ms_pepts) <- c("gene_symbol",
                              "protein_id",
                              "protein_length",
                              "support_peptides",
                              "neg_log10_p",
                              "coverage1",
                              "coverage2",
                              "area",
                              "support_peptides_total",
                              "support_peptides_unique",
                              "uniqueness",
                              "BSJsite(-s)",
                              "BSJpeptide(-s)",
                              colnames(x=ms_pepts)[14:18])
    ### Returning of the final results.
    return(ms_pepts)
}

setwd(dir="D:/Vasily Grinev")
newPepts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, new peptides, sequences.fasta"
msProts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, sequences.fasta"
msPepts="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, MS peptides, annotations.txt"
workDir="D:/Vasily Grinev"
res <- annotateCircRNAnewPeptides(newPepts=newPepts, msProts=msProts, msPepts=msPepts, workDir=workDir)
write.table(x=res,
            file="AML proteome, Kasumi-1, siMM for 3 days, #1, circRNA proteins, new peptides, annotations.txt",
            sep="\t",
            quote=FALSE,
            col.names=TRUE,
            row.names=FALSE)
