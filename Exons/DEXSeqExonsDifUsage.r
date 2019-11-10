###################################################################################################
##  A complete R-wrap for identification of differentially used exonic bins                      ##
##  using functionality of the DESeq2 and DEXSeq packages.                                       ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work     - character string giving the path to and name of work directory.
##  f.counts   - character string giving the name of TXT file in tab-delimited format containing
#                summarized RNA-Seq reads for exonic bins. It is typically output file of
#                exonicBinsCounts function.
##  thr.fpm    - numerical threshold value for FPM-based filtration of a count matrix.
#                Default value is 1.
##  ctrl       - character vector giving the list of control samples.
##  exper      - character vector giving the list of experimental samples.
##  f.diffBins - character string giving the name of TXT file in tab-delimited format for storing
#                of final results. Default value is NULL (no any file should be saved).
DEXSeqExonsDifUsage = function(d.work, f.counts, thr.fpm = 1, ctrl, exper, f.diffBins = NULL){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with the libraries DESeq2 v.1.24.0,
#   rtracklayer v.1.44.3 and DEXSeq v.1.30.0.
suppressMessages(require(DESeq2))
suppressMessages(require(rtracklayer))
suppressMessages(require(DEXSeq))
##  Loading of the count matrix for exonic bins.
counts = read.table(file = paste(d.work, f.counts, sep = "/"),
                    sep = "\t",
                    header = TRUE,
                    quote = "\"",
                    as.is = TRUE)
##  Filtering of the count matrix.
#   Filtering against to low sequencing depth.
FPM = data.frame(condition = factor(x = rep(x = c("control", "knockdown"),
                                            times = c(length(x = ctrl), length(x = exper))),
                                    levels = c("control", "knockdown")),
                 row.names = colnames(counts[8:ncol(counts)]))
FPM = fpm(object = DESeqDataSetFromMatrix(countData = counts[, 8:ncol(counts)],
                                          colData = FPM,
                                          design = ~ condition),
          robust = TRUE)
fpm.ctrl = counts[rowSums(FPM[, ctrl] >= thr.fpm) == length(ctrl), ]$exon_id
fpm.exper = counts[rowSums(FPM[, exper] >= thr.fpm) == length(exper), ]$exon_id
counts = counts[counts$exon_id %in% sort(x = unique(x = c(fpm.ctrl, fpm.exper))), ]
#   Filtering against one-bin genes.
counts = counts[counts$gene_id %in%
                attr((table(counts$gene_id) == 1)
                [(table(counts$gene_id) == 1) == "FALSE"], "dimnames")[[1]], ]
rownames(counts) = NULL
##  Developing an object of class DEXSeqDataSet.
#   Building a flattened file.
flat = counts[, 1:6]
flat$exonic_part_number = gsub(pattern = ".+EX", replacement = "EX", x = counts$exon_id)
flat$type = "exonic_part"
flat = flat[, c(2, 7:8, 3:6)]
flat = flat[order(flat$gene_id), ]
aggregate_gene = cbind(unique(x = flat$gene_id),
                       "",
                       "aggregate_gene",
                       unique(x = flat[, c(1, 4)])$seqnames)
aggregate_gene = cbind(aggregate_gene,
                       aggregate(x = flat$start,
                                 by = list(flat$gene_id),
                                 FUN = min)[, 2])
aggregate_gene = cbind(aggregate_gene,
                       aggregate(x = flat$end,
                                 by = list(flat$gene_id),
                                 FUN = max)[, 2])
aggregate_gene = cbind(aggregate_gene, unique(x = flat[, c(1, 7)])$strand)
colnames(aggregate_gene) = colnames(flat)
flat = rbind(flat, aggregate_gene)
flat$start = as.integer(x = flat$start)
flat$end = as.integer(x = flat$end)
flat = flat[order(flat$gene_id, flat$type, flat$start, flat$end), ]
flat = makeGRangesFromDataFrame(df = flat, keep.extra.columns = TRUE)
names(flat) = NULL
export(object = flat, con = paste(d.work, "flattened.gtf", sep = "/"), format = "gtf")
flat = list.files(path = d.work, pattern = "gtf$", full.names = TRUE)
#   Generation a set of count files.
for (i in 8:ncol(counts)){
write.table(x = counts[, c(1, i)],
            file = paste(paste(d.work, colnames(counts)[i], sep = "/"), ".txt", sep = ""),
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
}
#   List of individual files with counts.
counts = paste(d.work, paste(c(exper, ctrl), "txt", sep = "."), sep = "/")
#   Development a data frame object with a sample information.
smpls = data.frame(row.names = c(exper, ctrl),
                   condition = rep(x = c("knockdown", "control"),
                                   times = c(length(x = exper), length(x = ctrl))))
#   The DEXSeqDataSet object.
DXD = DEXSeqDataSetFromHTSeq(countfiles = counts,
                             sampleData = smpls,
                             design = ~ sample + exon + condition:exon,
                             flattenedfile = flat)
##  Identification of the differentially used exonic bins.
#   Estimation of size factors.
DXD = estimateSizeFactors(object = DXD)
#   Calculation of dispersion estimates.
DXD = estimateDispersions(object = DXD)
#   Testing for differential usage of exonic bins.
DXD = testForDEU(object = DXD)
DXD = estimateExonFoldChanges(object = DXD)
##  Consolidation and saving of final results.
res.table = DEXSeqResults(DXD)
res.table = cbind(res.table@rownames,
                  res.table@listData$groupID,
                  res.table@listData$log2fold_knockdown_control,
                  res.table@listData$pvalue)
colnames(res.table) = c("exon_id", "gene_id", "logFC", "p_value")
res.table = data.frame(res.table, stringsAsFactors = FALSE)
res.table$q_value = p.adjust(p = res.table$p_value, method = "BH", n = length(res.table$p_value))
res.table$exon_id = gsub(pattern = "EEXON", replacement = "EXON", x = res.table$exon_id)
res.table = res.table[order(res.table$exon_id), ]
if (!is.null(f.diffBins)){
write.table(x = res.table,
            file = paste(d.work, f.diffBins, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
##  Deleting of temporary files.
unlink(x = c(counts, paste(d.work, "flattened.gtf", sep = "/")))
##  Returning of the final object of class data frame.
return(res.table)
}
