#' Translate retained introns to peptide sequences
#' @description Translation of the retained introns into peptide sequences.
#' @param RIs character string giving the name of tab-delimited TXT file with
#'     genomic coordinates of retained introns. This file should include four
#'     mandatory fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of retained intron;
#'     iii) end    - end coordinate of retained intron;
#'     iv) strand  - strand information about retained intron.
#' @param genome character string giving the name of BSgenome data package with
#'     full genome sequences. Default value is "BSgenome.Hsapiens.UCSC.hg38".
#' @param codones location of tab-delimited TXT file with genetic code. Default
#'     value is "internal" for system extdata. Alternatively, the name of
#'     tab-delimited TXT file located in work directory should be specified.
#' @param aaSymbol type of amino acid symbols (one- or three-letter coding).
#'     Default value is 1.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return list of three AAStringSet objects with peptide sequences.
#' @author Vasily V. Grinev
#' @examples
#' RIs <- "AML proteome project, Kasumi-1, RIs.txt"
#' codones <- "codon.table.txt"
#' workDir <- "D:/Vasily Grinev"
#' prot_seq_RIs <- translateRIs(RIs=RIs, codones=codones, workDir=workDir)
#' @export
#' Last updated: July 14, 2023.

translateRIs <- function(RIs,
                         genome="BSgenome.Hsapiens.UCSC.hg38",
                         codones="internal",
                         aaSymbol=1, 
                         workDir=NULL){
    ### Loading of required packages.
    #   This code was successfully tested with packages GenomicRanges v.1.40.0,
    #   Biostrings v.2.56.0 and BSgenome.Hsapiens.UCSC.hg38 v.1.4.3.
    suppressMessages(expr=library(package=GenomicRanges))
    suppressMessages(expr=library(package=Biostrings))
    suppressMessages(expr=library(package=BSgenome.Hsapiens.UCSC.hg38))
    ### Setting the work directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Loading the gemonic coordinates of RIs as an object of class GRanges.
    gcRIs <- read.table(file=paste(workDir, RIs, sep = "/"),
                        sep="\t",
                        header=TRUE,
                        quote="\"",
                        as.is=TRUE)
    gcRIs <- makeGRangesFromDataFrame(df=gcRIs,
                                      keep.extra.columns=TRUE)
    ### Retrieving the sequences of RIs.
    seqRIs0 <- getSeq(x=get(x=genome), names=gcRIs)
    names(x=seqRIs0) <- paste(paste(seqnames(x=gcRIs),
                                    ranges(x=gcRIs),
                                    sep=":"),
                              strand(x=gcRIs),
                              sep="_str")
    seqRIs1 <- subseq(x=seqRIs0, start=2, end=width(x=seqRIs0))
    seqRIs2 <- subseq(x=seqRIs0, start=3, end=width(x=seqRIs0))
    ### Calling the codon table.
    if (codones == "internal"){
        codon_table_path <- system.file("extdata",
                                        "codon.table.txt",
                                        package="ORFhunteR")
    }else{
        codon_table_path <- paste(workDir, codones, sep = "/")
    }
    codon_table <- read.table(file=codon_table_path,
                              sep="\t",
                              header=TRUE,
                              quote="\"",
                              as.is=TRUE)
    ### RIs to proteins translation.
    aaRIs0 <- lapply(X=seqRIs0,
                     FUN=function(y){y <- as.character(x=y)
                                     prt <- sapply(X=seq(from=1,
                                                         to=nchar(x=y) - 3,
                                                         by=3),
                                            FUN=function(z){substr(x=y,
                                                            start=z,
                                                            stop=z + 2)
                                            }
                                     )
                                     prt <- paste(codon_table[match(x=prt,
                                                       table=codon_table[, 1]),
                                                  if (aaSymbol == 1){
                                                      4
                                                  }else{
                                                      3
                                                  }],
                                                  collapse="")
                                     }
    )
    aaRIs0 <- AAStringSet(x=unlist(x=aaRIs0), use.names=TRUE)
    aaRIs0 <- AAStringSet(x=gsub(pattern="\\*.+", replacement="", x=aaRIs0))
    aaRIs0 <- aaRIs0[width(x=aaRIs0) >= 20, ]
    aaRIs1 <- lapply(X=seqRIs1,
                     FUN=function(y){y <- as.character(x=y)
                                     prt <- sapply(X=seq(from=1,
                                                         to=nchar(x=y) - 3,
                                                         by=3),
                                                   FUN=function(z){substr(x=y,
                                                                   start=z,
                                                                   stop=z + 2)
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
    aaRIs1 <- AAStringSet(x=unlist(x=aaRIs1), use.names=TRUE)
    aaRIs1 <- AAStringSet(x=gsub(pattern="\\*.+", replacement="", x=aaRIs1))
    aaRIs1 <- aaRIs1[width(x=aaRIs1) >= 20, ]
    aaRIs2 <- lapply(X=seqRIs2,
                     FUN=function(y){y <- as.character(x=y)
                                     prt <- sapply(X=seq(from=1,
                                                         to=nchar(x=y) - 3,
                                                         by=3),
                                                   FUN=function(z){substr(x=y,
                                                                   start=z,
                                                                   stop=z + 2)
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
    aaRIs2 <- AAStringSet(x=unlist(x=aaRIs2), use.names=TRUE)
    aaRIs2 <- AAStringSet(x=gsub(pattern="\\*.+", replacement="", x=aaRIs2))
    aaRIs2 <- aaRIs2[width(x=aaRIs2) >= 20, ]
    ### Returning a final object.
    return(list(aaRIs0=aaRIs0, aaRIs1=aaRIs1, aaRIs2=aaRIs2))
}
