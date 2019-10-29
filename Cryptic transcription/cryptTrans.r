###################################################################################################
##  High-level R function for identification of cryptic intragenic transcription initiation      ##
##  based on RNA-Seq data.                                                                       ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work    - character string giving the path to and name of work directory.
##  f.counts  - character string giving the name of TXT file in tab-delimited format containing
#               summarized RNA-Seq reads for genomic intervals of interest.
##  ctrl      - character vector giving the list of control samples.
##  exper     - character vector giving the list of experimental samples.
##  thr.rpkm  - numerical threshold value for RPKM-based filtration of a count matrix.
#               Default value is 1.
##  multiExon - numerical threshold value for filtration of multi-exon genes. Default value is 5 
#               (all genes with at least 5 exons will be saved for subsequent processing).
cryptTrans = function(d.work, f.counts, ctrl, exper, thr.rpkm = 1, multiExon = 5){
##  Loading of required auxiliary library.
#   This code was successfully tested with library edgeR v.3.26.6.
suppressMessages(require(edgeR))
##  Loading of a primary count matrix.
counts = read.table(file = paste(d.work, f.counts, sep = "/"),
                    sep = "\t",
                    header = TRUE,
                    quote = "\"",
                    as.is = TRUE)
counts = counts[, c(colnames(counts)[1:6], ctrl, exper)]
##  Development of RPKM-based count matrix.
RPKM = DGEList(counts = counts[7:ncol(counts)],
               lib.size = calcNormFactors(object = counts[7:ncol(counts)], method = "TMM") *
                          colSums(counts[7:ncol(counts)]))
RPKM$genes$length = counts$width
counts[, 7:ncol(counts)] = data.frame(rpkm(RPKM))
##  Filtering against to low sequencing depth.
rpkm.ctrl = rownames(counts[rowSums(counts[, ctrl] >= thr.rpkm) == length(ctrl), ])
rpkm.exper = rownames(counts[rowSums(counts[, exper] >= thr.rpkm) == length(exper), ])
counts = counts[sort(x = unique(x = c(rpkm.ctrl, rpkm.exper))), ]
##  Filtering in of multi-exon genes.
counts = counts[counts$gene_id %in% attr((table(counts$gene_id) >= multiExon)
               [(table(counts$gene_id) >= multiExon) == "TRUE"], "dimnames")[[1]], ]
##  Exon numbering.
counts = split(x = counts, f = counts$gene_id)
counts = lapply(X = counts,
                FUN = function(x){if (x$strand[1] == "+"){
                                      x = x[order(x$start, x$end), ]
                                      x$exon_number = 1:nrow(x)
                                 }else{
                                      x = x[order(-x$end, -x$start), ]
                                      x$exon_number = 1:nrow(x)
                                 }
                                 x = x[, ]
                                 return(x)})
counts = do.call(what = rbind, args = counts)
rownames(counts) = NULL
counts = counts[, c(1, ncol(counts), 2:(ncol(counts) - 1))]
##  Calculation of a ratio between the first exon and all other exons of a gene.
counts = split(x = counts, f = counts$gene_id)
#   The ratio "second-to-first exon".
second = lapply(X = counts,
                FUN = function(x){x[x$exon_number == "2", 8:ncol(x)]/
                                  x[x$exon_number == "1", 8:ncol(x)]})
second = do.call(what = rbind, args = second)
second$gene_id = rownames(second)
rownames(second) = NULL
second = second[, c(ncol(second), 1:(ncol(second) - 1))]
#   The ratio "third-to-first exon".
third = lapply(X = counts,
               FUN = function(x){x[x$exon_number == "3", 8:ncol(x)]/
                                 x[x$exon_number == "1", 8:ncol(x)]})
third = do.call(what = rbind, args = third)
third$gene_id = rownames(third)
rownames(third) = NULL
third = third[, c(ncol(third), 1:(ncol(third) - 1))]
#   The ratio "intermediate-to-first exon".
inter = lapply(X = counts,
               FUN = function(x){apply(X = x[, 8:ncol(x)],
                                       MARGIN = 2,
                                       FUN = function(y){mean(y[4:(length(y) - 1)])/y[1]})})
inter = data.frame(do.call(what = rbind, args = inter), stringsAsFactors = FALSE)
inter$gene_id = rownames(inter)
rownames(inter) = NULL
inter = inter[, c(ncol(inter), 1:(ncol(inter) - 1))]
#   The ratio "last-to-first exon".
last = lapply(X = counts,
              FUN = function(x){x[nrow(x), 8:ncol(x)]/x[x$exon_number == "1", 8:ncol(x)]})
last = do.call(what = rbind, args = last)
last$gene_id = rownames(last)
rownames(last) = NULL
last = last[, c(ncol(last), 1:(ncol(last) - 1))]
##  Returning final results.
return(list(second = second, third = third, inter = inter, last = last))
}
