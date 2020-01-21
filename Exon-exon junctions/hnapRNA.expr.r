###################################################################################################
##  High-level R function for calculation of hypothetical "non-alternative"                      ##
##  precursor RNA based on experimentally determined expression level of exons.                  ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2020. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### It is experimental function without a guarantee of proper operation.
### Arguments of function:
##  d.work - character string giving the path to and name of work directory.
##  f.expr - character string giving the name of file in GTF/GFF format containing expression data
#            on transcripts for gene(-s) of interest. It is typically Cufflinks output file.
#   format - format of file with expression data. Default value is "gtf".
##  expr   - name of column (in GTF/GFF file) containing quantitative expression data.
hnapRNA.expr = function(d.work, f.expr, format = "gtf", expr){
##  Loading of required auxiliary library.
#   This code was successfully tested with library rtracklayer v.1.38.3.
suppressMessages(require(rtracklayer))
##  Loading and sub-setting of input data.
in.gtf = import(con = paste(d.work, f.expr, sep = "/"), format = format)
in.gtf = in.gtf[in.gtf$type == "exon", ]
##  Transforming of a GRanges object into GRangesList object.
gene.list = makeGRangesListFromDataFrame(df = as.data.frame(x = in.gtf),
                                         split.field = "gene_id",
                                         keep.extra.columns = TRUE)
hnapRNAs = lapply(X = gene.list, FUN = function(x){
##  Development a list of gene transcripts as a GRangesList object.
tr.list = sort(makeGRangesListFromDataFrame(df = as.data.frame(x = x),
                                            split.field = "transcript_id",
                                            keep.extra.columns = TRUE))
##  First exons.
#   Development a list of first exons as a GRanges object.
if (unique(as.character(strand(unlist(tr.list)))) == "+"){
    ex.first = do.call(what = getMethod(c, "GenomicRanges"),
                       args = lapply(X = tr.list,
                                     FUN = function(x){x[1, ]}))
}else{
    ex.first = do.call(what = getMethod(c, "GenomicRanges"),
                       args = lapply(X = tr.list,
                                     FUN = function(x){x[length(x), ]}))
}
#   Development a list of clusters of first exons as a GRanges object.
ex.first.cluster = intersect(range(ex.first), ex.first)
ex.first.cluster$cluster_id = paste("cl", c(1:length(ex.first.cluster)), sep = "")
#   Summarization of expression of first exons.
ex.first.expr = as.data.frame(x = ex.first[, names(elementMetadata(ex.first)) == expr])
ex.first.expr[, expr] = as.numeric(x = ex.first.expr[, expr])
ex.first.expr = aggregate(x = ex.first.expr[, expr],
                          by = list(ex.first.expr$seqnames,
                                    ex.first.expr$start,
                                    ex.first.expr$end,
                                    ex.first.expr$strand),
                          FUN = sum)
colnames(ex.first.expr) = c("seqnames", "start", "end", "strand", "expr")
ex.first.expr = makeGRangesFromDataFrame(df = ex.first.expr, keep.extra.columns = TRUE)
#   Development a final list of first exons as a GRanges object.
hits = findOverlaps(query = ex.first.expr, subject = ex.first.cluster, type = "within")
ex.first.expr$cluster_id = ""
ex.first.expr[queryHits(hits), ]$cluster_id = ex.first.cluster[subjectHits(hits), ]$cluster_id
ex.first.expr = sort(do.call(what = getMethod(c, "GenomicRanges"),
                             args = unlist(lapply(X = split(x = ex.first.expr,
                                                            f = ex.first.expr$cluster_id),
                                                  FUN = function(x){x[x$expr == max(x$expr), ]}))))
ex.first.expr = ex.first.expr[, 0]
if (unique(as.character(strand(unlist(ex.first.expr)))) == "+"){
    ex.first = ex.first.expr[1, ]
    ex.first$type = "first_exon"
    ex.alternative.first = ex.first.expr[-1, ]
    ex.alternative.first$type = "alternative_first_exon"
}else{
    ex.first = ex.first.expr[length(ex.first.expr), ]
    ex.first$type = "first_exon"
    ex.alternative.first = ex.first.expr[-length(ex.first.expr), ]
    ex.alternative.first$type = "alternative_first_exon"
}
##  Internal exons.
#   Development a list of internal exons as a GRanges object.
ex.inter = do.call(what = getMethod(c, "GenomicRanges"),
                   args = lapply(X = tr.list,
                                 FUN = function(x){x[c(-1, -length(x)), ]}))
#   Development a list of clusters of internal exons as a GRanges object.
ex.inter.cluster = intersect(range(ex.inter), ex.inter)
ex.inter.cluster$cluster_id = paste("cl", c(1:length(ex.inter.cluster)), sep = "")
#   Summarization of expression of internal exons.
ex.inter.expr = as.data.frame(x = ex.inter[, names(elementMetadata(ex.inter)) == expr])
ex.inter.expr[, expr] = as.numeric(x = ex.inter.expr[, expr])
ex.inter.expr = aggregate(x = ex.inter.expr[, expr],
                          by = list(ex.inter.expr$seqnames,
                                    ex.inter.expr$start,
                                    ex.inter.expr$end,
                                    ex.inter.expr$strand),
                          FUN = sum)
colnames(ex.inter.expr) = c("seqnames", "start", "end", "strand", "expr")
ex.inter.expr = makeGRangesFromDataFrame(df = ex.inter.expr, keep.extra.columns = TRUE)
#   Development a final list of internal exons as a GRanges object.
hits = findOverlaps(query = ex.inter.expr, subject = ex.inter.cluster, type = "within")
ex.inter.expr$cluster_id = ""
ex.inter.expr[queryHits(hits), ]$cluster_id = ex.inter.cluster[subjectHits(hits), ]$cluster_id
ex.inter.expr = sort(do.call(what = getMethod(c, "GenomicRanges"),
                             args = unlist(lapply(X = split(x = ex.inter.expr,
                                                            f = ex.inter.expr$cluster_id),
                                                  FUN = function(x){x[x$expr == max(x$expr), ]}))))
ex.inter.expr = ex.inter.expr[, 0]
ex.inter = ex.inter.expr
ex.inter$type = "exon"
##  Last exons.
#   Development a list of last exons as a GRanges object.
if (unique(as.character(strand(unlist(tr.list)))) == "+"){
    ex.last = do.call(what = getMethod(c, "GenomicRanges"),
                      args = lapply(X = tr.list,
                                    FUN = function(x){x[x[length(x), ], ]}))
}else{
    ex.last = do.call(what = getMethod(c, "GenomicRanges"),
                      args = lapply(X = tr.list,
                                    FUN = function(x){x[1, ]}))
}
#   Development a list of clusters of last exons as a GRanges object.
ex.last.cluster = intersect(range(ex.last), ex.last)
ex.last.cluster$cluster_id = paste("cl", c(1:length(ex.last.cluster)), sep = "")
#   Summarization of expression of last exons.
ex.last.expr = as.data.frame(x = ex.last[, names(elementMetadata(ex.last)) == expr])
ex.last.expr[, expr] = as.numeric(x = ex.last.expr[, expr])
ex.last.expr = aggregate(x = ex.last.expr[, expr],
                         by = list(ex.last.expr$seqnames,
                                   ex.last.expr$start,
                                   ex.last.expr$end,
                                   ex.last.expr$strand),
                         FUN = sum)
colnames(ex.last.expr) = c("seqnames", "start", "end", "strand", "expr")
ex.last.expr = makeGRangesFromDataFrame(df = ex.last.expr, keep.extra.columns = TRUE)
#   Development a final list of last exons as a GRanges object.
hits = findOverlaps(query = ex.last.expr, subject = ex.last.cluster, type = "within")
ex.last.expr$cluster_id = ""
ex.last.expr[queryHits(hits), ]$cluster_id = ex.last.cluster[subjectHits(hits), ]$cluster_id
ex.last.expr = sort(do.call(what = getMethod(c, "GenomicRanges"),
                            args = unlist(lapply(X = split(x = ex.last.expr,
                                                           f = ex.last.expr$cluster_id),
                                                 FUN = function(x){x[x$expr == max(x$expr), ]}))))
ex.last.expr = ex.last.expr[, 0]
if (unique(as.character(strand(unlist(ex.last.expr)))) == "+"){
    ex.last = ex.last.expr[length(ex.last.expr), ]
    ex.last$type = "last_exon"
    ex.alternative.last = ex.last.expr[-length(ex.last.expr), ]
    ex.alternative.last$type = "alternative_last_exon"
}else{
    ex.last = ex.last.expr[1, ]
    ex.last$type = "last_exon"
    ex.alternative.last = ex.last.expr[-1, ]
    ex.alternative.last$type = "alternative_last_exon"
}
##  Introns.
#   Development a list of introns as a GRangesList object.
intr = sort(c(ex.first, ex.inter, ex.last))
intr = setdiff(range(intr), intr)
intr$type = "intron"
##  Consolidated results.
model = sort(c(ex.first, ex.alternative.first, ex.inter, intr, ex.alternative.last, ex.last))
return(model)
})
gene.names = lengths(hnapRNAs, use.names = TRUE)
gene.names = rep(attr(gene.names, "names"), gene.names)
hnapRNAs = do.call(what = getMethod(c, "GenomicRanges"), args = unlist(hnapRNAs))
hnapRNAs$gene_id = gene.names
hnapRNAs$transcript_id = paste(hnapRNAs$gene_id, "hnapRNA", sep = ".")
return(hnapRNAs)
}
