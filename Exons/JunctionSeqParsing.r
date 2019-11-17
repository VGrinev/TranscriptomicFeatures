###################################################################################################
##  High-level R function for parsing main result files of JunctionSeq.                          ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work    - path to and name of work directory containing JunctionSeq output files.
##  f.results - character string giving the name of one of the main result files of JunctionSeq.
##  f.anno    - character string giving the name of TXT file in tab-delimited format containing
#               gene annotations. It is assuming that this file includes gene IDs as first column
#               and strand information as fifth column.
JunctionSeqParsing = function(d.work, f.results, f.anno){
##  Loading of result file as a standard object of class data frame.
res = read.table(file = paste(d.work, f.results, sep = "/"),
                 sep = "\t",
                 header = TRUE,
                 quote = "\"",
                 as.is = TRUE)
##  Selection of features with sufficient coverage.
res = res[res$testable == "TRUE", ]
##  Columns sub-setting.
res = res[, c(1:2, 14:17, 7, 23:24, 21:22, 12:13)]
colnames(res) = c("feature_id",
                  "gene_id",
                  "seqnames", "start", "end", "strand",
                  "base_mean", colnames(res)[8:9],
                  "logFC", "logFCvst",
                  "p_value", "q_value")
##  Assignment of strand to genes.
anno = read.table(file = paste(d.work, f.anno, sep = "/"),
                  sep = "\t",
                  header = TRUE,
                  quote = "\"",
                  as.is = TRUE)
rownames(anno) = anno[, 1]
anno = as.matrix(x = anno)
res = as.matrix(x = res)
idx = res[, 2] %in% anno[, 1]
res[idx, 6] = anno[res[idx, 2], 5]
res = data.frame(res, stringsAsFactors = FALSE)
res[, c(4:5, 7:13)] = apply(X = res[, c(4:5, 7:13)], MARGIN = 2, FUN = as.numeric)
##  Feature type sub-setting.
res.exons = res[gsub(".+:E.+", "E", res$feature_id) == "E", ]
res.exons$start = res.exons$start + 1
colnames(res.exons) = c("exon_id", colnames(res.exons)[-1])
res.exons = res.exons[order(res.exons$gene_id, res.exons$start, res.exons$end), ]
res.eej = res[gsub(".+:J.+", "J", res$feature_id) == "J", ]
res.eej$end = res.eej$end + 1
colnames(res.eej) = c("eej_id", colnames(res.eej)[-1])
res.eej = res.eej[order(res.eej$gene_id, res.eej$start, res.eej$end), ]
##  Returning of final object.
return(list(EXONS = res.exons, EEJs = res.eej))
}
