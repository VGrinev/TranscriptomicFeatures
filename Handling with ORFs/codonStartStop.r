#' Identify all potential start and stop codons in a nucleotide sequence
#' @description This function scans a nucleotide sequence of interest
#'     in the search of canonical start codon ATG or non-canonical start
#'     codons GTG, TTG and CTG as well as stop codons TAA, TAG and TGA.
#' @param x character string giving the nucleotide sequence.
#' @return list of potential start and stop codons with their coordinates.
#' @author Vasily V. Grinev
#' @examples
#' codons <- codonStartStop(x="AAAATGGCATGGTAAGTC")
#' @export

codonStartStop <- function(x){
### Calculation of codon positions.
codons <- DNAStringSet(x=c("ATG", "GTG", "TTG", "CTG", "TAA", "TAG", "TGA"))
names(x = codons) <- as.character(x=codons)
codonPositions <- sort(x = unlist(x=matchPDict(pdict=codons,
                                               subject=DNAString(x=x))))
codonPositions <- list(start(x=codonPositions), names(x=codonPositions))
### Returning a final object of class list.
return(codonPositions)
}
