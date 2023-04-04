#' Calculation of recall/sensitivity and precision/accuracy of a predictor of
#'     open reading frames for circular RNAs
#' @description Calculate the recall/sensitivity and precision/accuracy of the
#'     ORF predictor for circRNAs.
#' @param x character string giving the name of file with the nucleotide
#'     sequences of the reference circRNA ORFs. This file must be in the
#'     working directory. Allowed file formats are "fasta" or "fa".
#' @param codonTable character string giving the name of TXT file with genetic
#'     code in standard TAB-delimited format.
#' @param aaSymbol type of amino acid symbols (one- or three-letter coding).
#'     Default value is 1.
#' @param y character string giving the name of file with the nucleotide
#'     sequences of the empirical/predicted circRNA ORFs. This file must be in
#'     the working directory. Allowed file format is TAB-delimited TXT format
#'     with two obligatory fields "transcript_id" (identificator of circRNA)
#'     and "orf_sequence" (nucleotide sequence of empirical/predicted ORF).
#' @param size integer value giving the size for random sub-sampling of y.
#'     NULL by default (all events of y will be treated).
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @return An object of class data.frame containing three fields
#'     "emp_transcript_id", "recall" and "precision".
#' @authors Vasily V. Grinev, Dai Xiaoxuan
#' @examples
#' rec_prec <- calculateRecPrecORFpredictoR(x="TransCirc, circRNAs, ORF sequences.fa",
#'                                          codonTable="codon.table.txt",
#'                                          aaSymbol=1,
#'                                          y="TransCirc, circRNAs, ORFs by predictORFsCircRNAs.txt",
#'                                          size=1000,
#'                                          workDir="D:/Vasily Grinev")
#' @export
#' Last updated: April 4, 2023.

calculateRecPrecORFpredictoR <- function(x,
                                         codonTable,
                                         aaSymbol=1,
                                         y,
                                         size=NULL,
                                         workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with the base package parallel and
    #   packages Biostrings v.2.60.2.
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=parallel))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading of the codon table.
    codon_table <- read.table(file=paste(workDir, codonTable, sep="/"),
                              sep="\t",
                              header=TRUE,
                              quote="\"",
                              as.is=TRUE)
    ### Loading of the reference data.
    ##  Retrieving the file extension.
    frt <- tools::file_ext(x=x)
    ##  Validation of the file format.
    if (!frt %in% c("fa", "fasta")){
        stop("Invalid file format")
    }
    ##  Loading of the data.
    if (frt %in% c("fa", "fasta")){
        ref <- readDNAStringSet(filepath=paste(workDir, x, sep="/"))
    }
    ##  Subsetting of the reference data.
    ref <- ref[subseq(x=ref, start=1, end=3) == "ATG", ]
    ref <- ref[width(x=ref) < 60000, ]
    ref <- ref[sort(x=names(x=ref)), ]
    ref <- unique(x=ref)
    ### Reference open reading frames to proteins translation.
    ref_prts <- lapply(X=ref,
                       FUN=function(a){a <- as.character(x=a)
                                       prt <- sapply(X=seq(from=1,
                                                           to=nchar(x=a) - 3,
                                                           by=3),
                                                     FUN=function(b){substr(x=a,
                                                                            start=b,
                                                                            stop=b + 2)
                                                     })
                       prt <- paste(codon_table[match(x=prt,
                                                      table=codon_table[, 1]),
                                                if (aaSymbol == 1){
                                                    4
                                                }else{
                                                    3
                                                }],
                                    collapse="")
                       })
    ref_prts <- AAStringSet(x=unlist(x=ref_prts), use.names=TRUE)
    ### Loading of the empirical data.
    emp <- read.table(file=paste(workDir, y, sep="/"),
                      sep="\t",
                      header=TRUE,
                      quote="\"",
                      as.is=TRUE)
    emp <- emp[!is.na(x=emp$orf_sequence), ]
    names <- emp$transcript_id
    emp <- DNAStringSet(x=emp$orf_sequence)
    names(x=emp) <- names
    emp <- emp[width(x=emp) < 60000, ]
    emp <- emp[sort(x=names(x=emp)), ]
    emp <- unique(x=emp)
    if (!is.null(x=size)){
        emp <- emp[sample(x=1:length(x=emp), size=size), ]
    }
    ### Empirical open reading frames to proteins translation.
    emp_prts <- lapply(X=emp,
                       FUN=function(a){a <- as.character(x=a)
                                       prt <- sapply(X=seq(from=1,
                                                           to=nchar(x=a) - 3,
                                                           by=3),
                                                     FUN=function(b){substr(x=a,
                                                                            start=b,
                                                                            stop=b + 2)
                                                     })
                       prt <- paste(codon_table[match(x=prt,
                                                      table=codon_table[, 1]),
                                                if (aaSymbol == 1){
                                                    4
                                                }else{
                                                    3
                                                }],
                                    collapse="")
                       })
    emp_prts <- AAStringSet(x=unlist(x=emp_prts), use.names=TRUE)
    ### Counting of matchs.
    cl <- makeCluster(spec=detectCores() - 1)
    clusterExport(cl=cl,
                  varlist=c("vcountPattern", "length", "unique", "names",
                            "ref_prts", "paste", "sort", "c", "return"))
    matchs <- parLapply(X=emp_prts,
                        fun=function(d){recall <- vcountPattern(pattern=d,
                                                                subject=ref_prts)
                                        if (length(x=unique(x=recall)) == 1){
                                            precision <- "FALSE"
                                        }else{
                                            precision <- "TRUE"
                                        }
                                        recall <- names(x=ref_prts)[recall > 0]
                                        if (length(x=recall) == 0){
                                            recall <- ""
                                        }else{
                                            recall <- paste(sort(x=unique(x=recall)),
                                                            collapse=", ")
                                        }
                                        return(c(recall, precision))},
                        cl=cl)
    stopCluster(cl=cl)
    matchs <- cbind(names(x=matchs), do.call(what=rbind, args=matchs))
    rownames(x=matchs) <- NULL
    colnames(x=matchs) <- c("emp_transcript_id", "recall", "precision")
    matchs <- data.frame(matchs)
    ### Returning a final object of class data.frame.
    return(matchs)
}
