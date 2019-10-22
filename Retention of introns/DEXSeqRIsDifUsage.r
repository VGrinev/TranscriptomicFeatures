###################################################################################################
##  A complete R-wrap for identification of differentially retained introns                      ##
##  using functionality of the DESeq2 and DEXSeq packages.                                       ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work      - character string giving the path to and name of work directory.
##  f1.RIs      - character string giving the name of TXT file in tab-delimited format
#                 containing RIs identified in control cells by lmRIs function.
##  f2.RIs      - character string giving the name of TXT file in tab-delimited format
#                 containing RIs identified in experimental cells by lmRIs function.
##  f1.counts   - character string giving the name of TXT file in tab-delimited format
#                 containing summarized RNA-Seq reads for exonic bins. It is typically output
#                 file of genomicBinsCounts function.
##  f2.counts   - character string giving the name of TXT file in tab-delimited format
#                 containing summarized RNA-Seq reads for intronic bins. It is typically output
#                 file of genomicBinsCounts function.
##  ctrl        - character vector giving the list of control samples.
##  exper       - character vector giving the list of experimental samples.
##  f.exrpGenes - character string giving the name of TXT file in tab-delimited format
#                 with list of genes expressed in cells of interest.
##  thr.fpm     - numerical threshold value for FPM-based filtration of a count matrix.
#                 Default value is 1.
##  f.diffBins  - character string giving the name of TXT file in tab-delimited format for storing
#                 of final results. Default value is NULL (no any file should be saved).
DEXSeqRIsDifUsage = function(d.work, f1.RIs, f2.RIs, f1.counts, f2.counts, ctrl, exper,
                             f.exrpGenes, thr.fpm = 1, f.diffBins = NULL){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with the libraries DESeq2 v.1.24.0,
#   rtracklayer v.1.44.3 and DEXSeq v.1.30.0.
suppressMessages(require(DESeq2))
suppressMessages(require(rtracklayer))
suppressMessages(require(DEXSeq))
##  Loading of RIs from control cells.
RIs.ctrl = sort(x = unique(x = read.table(file = paste(d.work, f1.RIs, sep = "/"),
                                          sep = "\t",
                                          header = TRUE,
                                          quote = "\"",
                                          as.is = TRUE)$intron_id))
##  Loading of the RIs from experimental cells.
RIs.exper = sort(x = unique(x = read.table(file = paste(d.work, f1.RIs, sep = "/"),
                                           sep = "\t",
                                           header = TRUE,
                                           quote = "\"",
                                           as.is = TRUE)$intron_id))
##  Loading of the count matrix for exonic bins.
counts.exon = read.table(file = paste(d.work, f1.counts, sep = "/"),
                         sep = "\t",
                         header = TRUE,
                         quote = "\"",
                         as.is = TRUE)
colnames(counts.exon) = c("bin_id", colnames(counts.exon)[-1])
##  Loading of the count matrix for intronic bins.
counts.intr = read.table(file = paste(d.work, f2.counts, sep = "/"),
                         sep = "\t",
                         header = TRUE,
                         quote = "\"",
                         as.is = TRUE)
colnames(counts.intr) = c("bin_id", colnames(counts.intr)[-1])
counts.intr = counts.intr[, c(-2, -9)]
counts.intr = counts.intr[counts.intr$bin_id %in% sort(x = unique(x = c(RIs.ctrl, RIs.exper))), ]
##  Development of the final (combined) count matrix.
counts = rbind(counts.exon[, -3:-7], counts.intr[, -3:-7])
counts = counts[, c("bin_id", "gene_id", ctrl, exper)]
##  Filtering of the count matrix.
#   Filtering against non expressed genes.
exrpGenes = read.table(file = paste(d.work, f.exrpGenes, sep = "/"),
                       sep = "\t",
                       header = TRUE,
                       quote = "\"",
                       as.is = TRUE)$gene_id
counts = counts[counts$gene_id %in% exrpGenes, ]
#   Filtering against to low sequencing depth.
FPM = data.frame(condition = factor(x = rep(x = c("control", "knockdown"),
                                            times = c(length(x = ctrl), length(x = exper))),
                                    levels = c("control", "knockdown")),
                 row.names = colnames(counts[3:ncol(counts)]))
FPM = fpm(object = DESeqDataSetFromMatrix(countData = counts[, 3:ncol(counts)],
                                          colData = FPM,
                                          design = ~ condition),
          robust = TRUE)
fpm.ctrl = counts[rowSums(FPM[, ctrl] >= thr.fpm) == length(ctrl), ]$bin_id
fpm.exper = counts[rowSums(FPM[, exper] >= thr.fpm) == length(exper), ]$bin_id
counts = counts[counts$bin_id %in% sort(x = unique(x = c(fpm.ctrl, fpm.exper))), ]
#   Filtering against one-bin genes.
counts = counts[counts$gene_id %in%
                attr((table(counts$gene_id) == 1)
                [(table(counts$gene_id) == 1) == "FALSE"], "dimnames")[[1]], ]
rownames(counts) = NULL
##  Developing an object of class DEXSeqDataSet.
#   Building a flattened file.
flat = rbind(counts.exon, counts.intr)[, c(1:6)]
flat = flat[flat$bin_id %in% counts$bin_id, ]
flat$exonic_part_number = gsub(pattern = ".+EX", replacement = "EX", x = flat$bin_id)
flat$exonic_part_number = gsub(pattern = ".+IN", replacement = "IN", x = flat$exonic_part_number)
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
for (i in 3:ncol(counts)){
write.table(x = counts[, c(1, i)],
            file = paste(paste(d.work, colnames(counts)[i], sep = "/"), ".txt", sep = ""),
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
}
#   List of individual files with counts.
f.counts = paste(d.work, paste(c(exper, ctrl), "txt", sep = "."), sep = "/")
#   Development a data frame object with a sample information.
smpls = data.frame(row.names = c(exper, ctrl),
                   condition = rep(x = c("knockdown", "control"),
                                   times = c(length(x = exper), length(x = ctrl))))
#   The DEXSeqDataSet object.
DXD = DEXSeqDataSetFromHTSeq(countfiles = f.counts,
                             sampleData = smpls,
                             design = ~ sample + exon + condition:exon,
                             flattenedfile = flat)
##  Identification of the differentially used retained introns.
#   Estimation of size factors.
DXD = estimateSizeFactors(object = DXD)
#   Calculation of dispersion estimates.
DXD = estimateDispersions(object = DXD)
#   Testing for differential intron usage.
DXD = testForDEU(object = DXD)
DXD = estimateExonFoldChanges(object = DXD)
##  Consolidation and saving of final results.
res.table = DEXSeqResults(DXD)
res.table = cbind(res.table@rownames,
                  res.table@listData$groupID,
                  res.table@listData$log2fold_knockdown_control,
                  res.table@listData$pvalue)
colnames(res.table) = c("bin_id", "gene_id", "logFC", "p_value")
res.table = data.frame(res.table, stringsAsFactors = FALSE)
res.table$q_value = p.adjust(p = res.table$p_value, method = "BH", n = length(res.table$p_value))
res.table$bin_id = gsub(pattern = "EEXON", replacement = "EXON", x = res.table$bin_id)
res.table$bin_id = gsub(pattern = "EINTR", replacement = "INTR", x = res.table$bin_id)
res.table = res.table[order(res.table$bin_id), ]
if (!is.null(f.diffBins)){
write.table(x = res.table,
            file = paste(d.work, f.diffBins, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
##  Deleting of temporary files.
unlink(x = c(f.counts, paste(d.work, "flattened.gtf", sep = "/")))
##  Returning of the final object of class data frame.
return(res.table)
}
### A simple example of function use.
res = DEXSeqRIsDifUsage(d.work = "D:/Vasily Grinev",
                        f1.RIs = "The NML project, DESeq2-based RIs, NC transcriptome.txt",
#                        f2.RIs = "The NML project, DESeq2-based RIs, KAT8 transcriptome.txt",
                        f2.RIs = "The NML project, DESeq2-based RIs, KANSL2 transcriptome.txt",
                        f1.counts = "The NML project, exonic bins, primary matrix of read counts.txt",
                        f2.counts = "The NML project, intronic bins, primary matrix of read counts.txt",
                        ctrl = c("NC_S1", "NC_S2", "NC_S3"),
#                        exper = c("KAT8_S1", "KAT8_S2", "KAT8_S3"),
                        exper = c("KANSL2_S1", "KANSL2_S2"),
                        f.exrpGenes = "The NSL project, list of expressed genes.txt",
                        thr.fpm = 1,
#                        f.diffBins = "The NML project, DEXSeq-based diffBins, KAT8 versus NC.txt")
                        f.diffBins = "The NML project, DEXSeq-based diffBins, KANSL2 versus NC.txt")
