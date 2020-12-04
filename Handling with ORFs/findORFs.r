#' Identify all variants of open reading frames in a nucleotide sequence
#' @description Identify all possible variants of open reading frames
#'     in a nucleotide sequence of interest.
#' @param x character string giving the nucleotide sequence of interest.
#' @param codStart character string with type of start codon: "ATG", "GTG",
#'      "TTG" or "CTG". Default value is "ATG".
#' @return matrix with start and stop positions, length and sequence of
#'     identified variants of open reading frames.
#' @author Vasily V. Grinev
#' @examples
#' x <- "AAAATGGCTGCGTAATGCAAAATGGCTGCGAATGCAAAATGGCTGCGAATGCCGGCACGTTGCTACGT"
#' orf <- findORFs(x=x, codStart="ATG")
#' @export

findORFs <- function(x, codStart="ATG"){
### Calculation coordinates of defined start codon and all stop codons 
#   in the sequence of interest.
allPos <- codonStartStop(x=x)
allStarts <- allPos[[1]][allPos[[2]] == codStart]
allStops <- allPos[[1]][allPos[[2]] %in% c("TAA", "TAG", "TGA")]
### Development a matrix of in-frame calculated coordinates.
inFrame <- outer(X=allStops, Y=allStarts, FUN="-")
L <- inFrame[lower.tri(x=inFrame, diag=TRUE)]/3
L <- round(x=L) == L & L > 0
if (length(x=L[L == "TRUE"]) > 0){
    ### Calculation of the in-frame start codons.
    inFrameStarts <- t(x=inFrame)
    inFrameStarts[] <- allStarts
    inFrameStarts <- t(x=inFrameStarts)
    inFrameStarts <- inFrameStarts[lower.tri(x=inFrameStarts, diag=TRUE)]
    inFrameStarts <- inFrameStarts[L]
    ### Calculation of the in-frame stop codons.
    inFrameStops <- inFrame
    inFrameStops[] <- allStops
    inFrameStops <- inFrameStops[lower.tri(x=inFrameStops, diag=TRUE)]
    inFrameStops <- inFrameStops[L]
    ### Calculation the coordinates of open reading frames.
    orfs <- cbind(inFrameStarts, inFrameStops)
    orfs <- aggregate(x=orfs[, 2], by=list(orfs[, 1]), FUN=min)
    orfs <- aggregate(x=orfs[, 1], by=list(orfs[, 2]), FUN=min)
    orfs <- orfs[order(orfs[, 2]), ]
    seqorfs <- substring(x, first=orfs[, 2], last=orfs[, 1] + 2)
    ### Final object.
    orfs <- list(orfs[, 2], orfs[, 1] + 2, orfs[, 1] - orfs[, 2] + 3, seqorfs)
    orfs <- do.call(what=cbind, args=orfs)
    colnames(orfs) <- c("start", "end", "length", "orf.sequence")
    return(orfs)
    }
}
