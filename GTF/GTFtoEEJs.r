###################################################################################################
##  High-level R function for retrieving genomic coordinates of exon-exon junctions              ##
##  from GTF/GFF infrastructure.                                                                 ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  x - path to folder and name of the GTF/GFF file with gene annotations.
##  with.metadata - logical argument. Default value is FALSE. If TRUE, all exon-exon junctions
#   will be annotated with metadata from the GTF/GFF file.
GTFtoEEJs = function(x, with.metadata = FALSE){
##  Loading of required auxiliary library.
#   This code was successfully tested with library rtracklayer v.1.38.3 (or higher).
if (!suppressMessages(require(rtracklayer))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
}
##  Loading of input GTF/GFF file.
in.gtf = import(con = x, format = "gtf")
in.gtf = in.gtf[in.gtf$type == "exon", ]
##  Development a list of gene transcripts as a GRangesList object.
tr.list = makeGRangesListFromDataFrame(df = as.data.frame(x = in.gtf),
                                       split.field = "transcript_id",
                                       keep.extra.columns = TRUE)
#   Development a list of introns as a GRangesList object.
intr = setdiff(range(tr.list), tr.list)
#   Development a list of exon-exon junctions as a GRanges object.
EEJs = unlist(intr)
start(EEJs) = start(EEJs) - 1
end(EEJs) = end(EEJs) + 1
EEJs$eej_id = paste(paste(paste(seqnames(EEJs),
                                start(EEJs),
                                sep = ":"),
                          end(EEJs),
                          sep = "-"),
                    strand(EEJs),
                    sep = "_str")
if (with.metadata == FALSE){
EEJs = sort(EEJs[!duplicated(EEJs), ])
names(EEJs) = NULL
}else{
EEJs$transcript_id = names(EEJs)
names(EEJs) = NULL
EEJs = as.matrix(as.data.frame(EEJs))
EEJs = cbind(EEJs, matrix(nrow = nrow(EEJs), ncol = length(elementMetadata(in.gtf)) - 1))
tr.anno = as.data.frame(in.gtf)[, -1:-5]
tr.anno = tr.anno[!duplicated(tr.anno), ]
tr.anno = tr.anno[, c(which(colnames(tr.anno) == "transcript_id"),
                      which(colnames(tr.anno)!= "transcript_id"))]
rownames(tr.anno) = tr.anno[, 1]
tr.anno = as.matrix(tr.anno)
idx = EEJs[, 7] %in% tr.anno[, 1]
EEJs[idx, 8:dim(EEJs)[[2]]] = tr.anno[EEJs[idx, 7], 2:dim(tr.anno)[[2]]]
EEJs = data.frame(EEJs, stringsAsFactors = FALSE)
colnames(EEJs) = c(colnames(EEJs)[1:7], attr(tr.anno, "dimnames")[[2]][-1])
EEJs.metadata = aggregate(x = EEJs[, 7:ncol(EEJs)],
                          by = list(EEJs$eej_id),
                          function(x){paste(x, collapse = ", ")})
EEJs = EEJs[, 1:6]
EEJs = EEJs[, c(which(colnames(EEJs) == "eej_id"),
                which(colnames(EEJs)!= "eej_id"))]
EEJs = EEJs[!duplicated(EEJs), ]
EEJs = EEJs[order(EEJs$eej_id), ] 
EEJs = cbind(EEJs, EEJs.metadata[, -1])
EEJs = makeGRangesFromDataFrame(df = EEJs, keep.extra.columns = TRUE)
}
return(EEJs)
}
            col.names = TRUE,
            row.names = FALSE)
