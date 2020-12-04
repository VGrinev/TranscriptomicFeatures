#' Translate open reading frames to proteins
#' @description Translate the identified open reading frames to proteins.
#' @param seqORFs character string giving the name of FASTA/FA file with
#'     sequences of identified open reading frames.
#' @param aaSymbol type of amino acid symbols (one- or three-letter coding).
#'     Default value is 1.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return AAStringSet object with protein sequences.
#' @author Vasily V. Grinev
#' @examples
#' seq_orf_path <- system.file("extdata",
#'                             "Set.trans_ORFs.sequences.fasta",
#'                             package="ORFhunteR")
#' prot_seqs <- translateORFs(seqORFs=seq_orf_path)
#' @export

translateORFs <- function(seqORFs, aaSymbol=1, workDir=NULL){
### Internal calling of codon table.
codon_table_path <- system.file("extdata",
                                "codon.table.txt",
                                package="ORFhunteR")
codon_table <- read.table(file = codon_table_path,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          as.is = TRUE)
### Loading of the open reading frame sequences
#   as an object of class DNAStringSet.
if (!is.null(workDir)){
    seq_orfs <- paste(workDir, seqORFs, sep = "/")
}
seq_orfs <- readDNAStringSet(filepath=seq_orfs)
### Open reading frames to proteins translation.
seq_prts <- lapply(X=seq_orfs,
                   FUN = function(y){y <- as.character(x=y)
                                     prt <- sapply(X=seq(from=1,
                                                         to=nchar(x=y) - 3,
                                                         by=3),
                                                   FUN=function(z){substr(x=y,
                                                                   start = z,
                                                                   stop=z + 2)
                                                   }
                                     )
                                     prt <- paste(codon_table[
                                                  match(x=prt,
                                                        table=codon_table[, 1]),
                                                  if (aaSymbol == 1){
                                                      4
                                                  }else{
                                                      3
                                                  }],
                                                  collapse="")
                   }
)
seq_prts <- AAStringSet(x=unlist(x=seq_prts), use.names=TRUE)
### Returning a final object.
return(seq_prts)
}
