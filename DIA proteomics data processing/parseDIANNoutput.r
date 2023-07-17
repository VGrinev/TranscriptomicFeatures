#' Parse of DIA-NN output data
#' @description Parsing of the DIA-NN output data.
#' @param x character string giving the name of main output DIA-NN matrix. This
#'     file must be in the working directory. Valid file format can be either
#'     "txt" or "tsv".
#' @param groups list of experimental groups of samples. NULL by default.
#' @param y character string giving the name of DIA-NN pr matrix. This file
#'     must be in the working directory. Valid file format can be either "txt"
#'     or "tsv".
#' @param fdr FDR threshold for filtration of precursors and proteins.
#'     Default value is 0.05.
#' @param pep PEP threshold for filtration of precursors and proteins.
#'     Default value is 0.05.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return List containing general filtering statistics and the main data frame
#'     with the following fields:
#'     i) protein_id             - protein ID;
#'     ii) gene_symbol           - gene symbol;
#'     iii) support_peptides_seq - sequence of supporting peptides;
#'     iv) support_peptides_n    - quantity of supporting peptides;
#'     v-...) ...                - one or several sample fields with intensity
#'                                 values.
#' @author Vasily V. Grinev
#' @examples
#' res <- parseDIANNoutput(x="kasumi_6d_fdr5.tsv",
#'                         groups=list(siMM=c("Kasumi_siMM_6d_1",
#'                                            "Kasumi_siMM_6d_2",
#'                                            "Kasumi_siMM_6d_3",
#'                                            "Kasumi_siMM_6d_4",
#'                                            "Kasumi_siMM_6d_5"),
#'                                     siRR=c("Kasumi_siRR_6d_1",
#'                                            "Kasumi_siRR_6d_2",
#'                                            "Kasumi_siRR_6d_3",
#'                                            "Kasumi_siRR_6d_4",
#'                                            "Kasumi_siRR_6d_5")),
#'                         y="kasumi_6d_fdr5.pr_matrix.tsv",
#'                         fdr=0.05,
#'                         pep=0.05,
#'                         workDir="D:/Vasily Grinev")
#' @export
#' Last updated: July 9, 2023.

parseDIANNoutput <- function(x, groups=NULL, y, fdr=0.05, pep, workDir=NULL){
    ### Loading the required package.
    #   This code was successfully tested with package stringr v.1.4.0.
    suppressMessages(expr=library(package=stringr))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading the main output DIA-NN matrix as an object of class data frame.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=x)
    ##  Validation of the file format.
    if (!frt %in% c("txt", "tsv")){
        stop("Invalid file format")
    }
    ##  Loading of the data.
    runs <- read.table(file=paste(workDir, x, sep = "/"),
                       sep="\t",
                       header=TRUE,
                       quote="\"",
                       as.is=TRUE)
    ### Loading the  DIA-NN pr_matrix as an object of class data frame.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=y)
    ##  Validation of the file format.
    if (!frt %in% c("txt", "tsv")){
        stop("Invalid file format")
    }
    ##  Loading of the data.
    prots <- read.table(file=paste(workDir, y, sep = "/"),
                        sep="\t",
                        header=TRUE,
                        quote="\"",
                        as.is=TRUE)
    prots <- prots[, c(2, 4, 7, 11:ncol(x=prots))]
    colnames(x=prots) <- c("protein_id",
                           "gene_symbol", 
                           "peptide_stripped",
                           colnames(x=prots)[-1:-3])
    if (!is.null(x=groups)){
        prots <- prots[, c(colnames(x=prots)[1:3], unlist(x=groups))]
    }
    ### Collection of general statistics.
    stats <- list()
    genes <- length(x=unique(x=unlist(x=strsplit(x=prots$gene_symbol,
                                                 split=";"))))
    proteins <- length(x=unique(x=unlist(x=strsplit(x=prots$protein_id,
                                                    split=";"))))
    peptides <- length(x=unique(x=unlist(x=strsplit(x=prots$peptide_stripped,
                                                    split=";"))))
    stats[[1]] <- c(genes=genes, proteins=proteins, peptides=peptides)
    ### Filtration.
    ##  First filter: run-specific precursor q-value.
    filter1 <- sort(x=unique(x=runs[runs$Q.Value <= fdr, ]$Stripped.Sequence))
    prots <- prots[prots$peptide_stripped %in% filter1, ]
    ##  Collection of step specific statistics.
    genes <- length(x=unique(x=unlist(x=strsplit(x=prots$gene_symbol,
                                                 split=";"))))
    proteins <- length(x=unique(x=unlist(x=strsplit(x=prots$protein_id,
                                                    split=";"))))
    peptides <- length(x=unique(x=unlist(x=strsplit(x=prots$peptide_stripped,
                                                    split=";"))))
    stats[[2]] <- c(genes=genes, proteins=proteins, peptides=peptides)
    ##  Second filter: the posterior error probability for the precursor
    #   identification.
    filter2 <- sort(x=unique(x=runs[runs$PEP <= pep, ]$Stripped.Sequence))
    prots <- prots[prots$peptide_stripped %in% filter2, ]
    ##  Collection of step specific statistics.
    genes <- length(x=unique(x=unlist(x=strsplit(x=prots$gene_symbol,
                                                 split=";"))))
    proteins <- length(x=unique(x=unlist(x=strsplit(x=prots$protein_id,
                                                    split=";"))))
    peptides <- length(x=unique(x=unlist(x=strsplit(x=prots$peptide_stripped,
                                                    split=";"))))
    stats[[3]] <- c(genes=genes, proteins=proteins, peptides=peptides)
    ##  Third filter: run and experimental group specific precursor quantity.
    if (!is.null(x=groups)){
        filter3 <- list()
        for (i in 1:length(x=groups)){
            intens <- prots[, c("peptide_stripped", groups[[i]])]
            intens <- intens[rowSums(x=(intens[, -1] > 0),
                                     na.rm=TRUE) >= ncol(x=intens[, -1])/2,
                            ]$peptide_stripped
            filter3[[i]] <- intens
        }
        filter3 <- sort(x=unique(x=unlist(x=filter3)))
    }else{
        filter3 <- prots[rowSums(x=(prots[, -1:-3] > 0),
                                     na.rm=TRUE) > ncol(x=prots[, -1:-3])/2,
                            ]$peptide_stripped
        filter3 <- sort(x=unique(x=filter3))
    }
    prots <- prots[prots$peptide_stripped %in% filter3, ]
    ##  Collection of step specific statistics.
    genes <- length(x=unique(x=unlist(x=strsplit(x=prots$gene_symbol,
                                                 split=";"))))
    proteins <- length(x=unique(x=unlist(x=strsplit(x=prots$protein_id,
                                                    split=";"))))
    peptides <- length(x=unique(x=unlist(x=strsplit(x=prots$peptide_stripped,
                                                    split=";"))))
    stats[[4]] <- c(genes=genes, proteins=proteins, peptides=peptides)
    ##  Fourth filter: run-specific q-value for the unique protein, that is
    #   protein identified with proteotypic (specific to it) peptide(-s).
    filter4 <- sort(x=unique(x=runs[runs$Protein.Q.Value <= fdr, ]$Protein.Ids))
    prots <- prots[prots$protein_id %in% filter4, ]
    ##  Collection of step specific statistics.
    genes <- length(x=unique(x=unlist(x=strsplit(x=prots$gene_symbol,
                                                 split=";"))))
    proteins <- length(x=unique(x=unlist(x=strsplit(x=prots$protein_id,
                                                    split=";"))))
    peptides <- length(x=unique(x=unlist(x=strsplit(x=prots$peptide_stripped,
                                                    split=";"))))
    stats[[5]] <- c(genes=genes, proteins=proteins, peptides=peptides)
    ### Collecting of the proteins and supporting peptides.
    pr <- list()
    for (i in 1:nrow(x=prots)){
        suppressWarnings(expr=pr[[i]] <- cbind(strsplit(x=prots$protein_id[i],
                                                        split=";")[[1]],
                                               prots[i, -1]))
    }
    prots <- do.call(what=rbind, args=pr)
    prots <- cbind(aggregate(x=prots[, 2:3],
                             by=list(prots[, 1]),
                             FUN=function(y){paste(sort(x=unique(x=y)),
                                                   collapse=";")}),
                   aggregate(x=prots[, 4:ncol(x=prots)],
                             by=list(prots[, 1]),
                             FUN=function(y){mean(y, na.rm=TRUE)})[, -1])
    z <- function(y){
                     paste(sort(x=unique(x=strsplit(x=y, split=";")[[1]])),
                           collapse=";")
         }
    prots$gene_symbol <- unlist(x=lapply(X=prots$gene_symbol, FUN=z))
    prots$peptide_stripped <- unlist(x=lapply(X=prots$peptide_stripped, FUN=z))
    prots$support_peptides_n <- str_count(string=prots$peptide_stripped,
                                          pattern=";") + 1
    prots <- prots[, c(1:3, ncol(x=prots), 4:(ncol(x=prots) - 1))]
    colnames(x=prots) <- c("protein_id",
                           "gene_symbol",
                           "support_peptides_seq",
                           "support_peptides_n",
                           colnames(x=prots)[-1:-4])
    ### The last filtration: number of support peptides.
    prots <- prots[prots$support_peptides_n > 1, ]
    ##  Collection of step specific statistics.
    genes <- length(x=unique(x=unlist(x=strsplit(x=prots$gene_symbol,
                                                 split=";"))))
    proteins <- length(x=unique(x=unlist(x=strsplit(x=prots$protein_id,
                                                    split=";"))))
    peptides <- length(x=unique(x=unlist(x=strsplit(x=prots$support_peptides_seq,
                                                    split=";"))))
    stats[[6]] <- c(genes=genes, proteins=proteins, peptides=peptides)
    names(x=stats) <- c("input", paste("step", 1:5, sep=""))
    ### Returning of the final results.
    return(list(proteins=prots, statistics=stats))
}
