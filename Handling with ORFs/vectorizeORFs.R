#' Vectorize the ORF(-s) sequence(-s) to sequence features
#' @description Vectorize the ORF(-s) sequence(-s) to sequence features.
#' @param x DNAStringSet object with ORF(-s) sequence(-s).
#' @return The object of class data.frame
#'     with ORF(-s) vectorized into sequence features.
#' @author Mikalai M. Yatskou
#' @examples
#' \dontrun{
#' feats <- vectorizeORFs(x=DNAStringSet(x="ATGGGCCTCA"))
#' }
#' @export

vectorizeORFs = function(x){
### A vector of dinucleotides for calculation of correlation factors.
diNucls = c("AT", "AG", "AC", "TG", "TC", "GC")
### Calculation of correlation factors.
corFcs <- lapply(X=diNucls,
                 FUN=function(y){do.call(what=c,
                                         args=sapply(X=as.character(x),
                                                     FUN=getCorrelationFactors,
                                                     pattern=y,
                                                     blockSize=10,
                                                     simplify=FALSE,
                                                     USE.NAMES=FALSE))})
corFcs <- data.frame(do.call(what=cbind, args=corFcs))
### Calculation of Bao metrics.
bao <- do.call(what=rbind, args=sapply(X=as.character(x),
                                       FUN=getBaoMetrics,
                                       indexes=length(x=x),
                                       types=length(x=x),
                                       simplify=FALSE,
                                       USE.NAMES=FALSE))[, -1:-2]
### Calculation of oligonucleotides frequency (1-mer to 3-mer).
mo <- data.frame(oligonucleotideFrequency(x=x, width=1, as.prob=TRUE))
di <- data.frame(oligonucleotideFrequency(x=x, width=2, as.prob=TRUE))
tr <- data.frame(oligonucleotideFrequency(x=x, width=3, as.prob=TRUE))
### Consolidation of all ORF features in a single object of class data.frame.
features <- cbind(mo, di, tr, corFcs, nchar(x=x), log10(x=nchar(x=x)), bao)
rownames(x=features) <- NULL
colnames(x=features) <- c(paste("freq", colnames(x=features)[1:84], sep=""),
                           paste("corrFact", diNucls, sep=""),
                           "length", "logLength",
                           colnames(x=features)[93:104])
### Returning finnal results.
return(features)
}
