###################################################################################################
##  A complete R-wrap for identification                                                         ##
##  of retained introns in RNA-Seq data with the DESeq2 package.                                 ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work      - character string giving the path to and name of work directory.
##  f.counts    - character string giving the name of TXT file in tab-delimited format with
#                 summarized RNA-Seq reads for intronic bins. It is typically output file of
#                 genomicBinsCounts function.
##  f.exrpGenes - character string giving the name of TXT file in tab-delimited format
#                 with list of genes expressed in cells of interest.
##  rL          - integer specifying the length of RNA-Seq reads.
##  f.null      - character string giving the name of TXT file in tab-delimited format for storing
#                 null distributions. Default value is NULL (no any file should be saved).
##  groups      - a list of named groups of samples to be analysed.
##  thr.prob    - numeric value of probability for quantile-based threshold of intron abundance.
#                 Default value is 0.97.
##  thr.logFC   - numeric value for fold-change threshold. Default value is 1.
##  thr.q.value - numerical threshold value for statistical significance. Default value is 0.01.
##  f.RIs       - character string giving the name of TXT file in tab-delimited format for storing
#                 of final results. Default value is NULL (no any file should be saved).
DESeq2RIs = function(d.work, f.counts, f.exrpGenes, rL, f.null = NULL, groups,
                     thr.prob = 0.97, thr.logFC = 1, thr.q.value = 0.01, f.RIs = NULL){
##  Loading of a required auxiliary library.
#   This code was successfully tested with the library DESeq2 v.1.24.
suppressMessages(require(DESeq2))
##  Loading of a count matrix for intronic bins.
counts = read.table(file = paste(d.work, f.counts, sep = "/"),
                    sep = "\t",
                    header = TRUE,
                    quote = "\"",
                    as.is = TRUE)
##  Filtering of the count matrix.
#   Filtering against non expressed genes.
exrpGenes = read.table(file = paste(d.work, f.exrpGenes, sep = "/"),
                       sep = "\t",
                       header = TRUE,
                       quote = "\"",
                       as.is = TRUE)$gene_id
counts = counts[counts$gene_id %in% exrpGenes, ]
#   Filtering against to short effective length of introns.
counts = counts[counts$effective_length >= rL, ]
#   Filtering against one-intron genes.
counts = counts[counts$gene_id %in%
                        attr((table(counts$gene_id) == 1)[(table(counts$gene_id) == 1) == "FALSE"],
                             "dimnames")[[1]], ]
rownames(counts) = NULL
##  Extending of the count matrix with columns of null distributions.
n.cols = ncol(counts)
counts[c(paste(colnames(counts)[10:n.cols], "NULL", sep = "."))] = 0
##  Converting of the count matrix in gene-wise list of data frames.
counts = split(x = counts, f = counts$gene_id)
##  Generation of the null distributions.
counts = lapply(X = counts,
                FUN = function(z){iL = sqrt(z[9] - rL)
                                  r.sum = colSums(z[10:n.cols])
                                  null = do.call(what = cbind,
                                                 args = lapply(X = r.sum,
                                                 FUN = function(y){round(x = (y * iL/sum(iL)))}))
                                  z[, (n.cols + 1):ncol(z)] = null
                                  return(z)})
counts = do.call(what = rbind, args = counts)
rownames(counts) = NULL
if (!is.null(f.null)){
write.table(x = counts,
            file = paste(d.work, f.null, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
##  Wrapping of the count matrix and sample information in a DESeqDataSet object.
#   Development a data frame object with a sample information.
smpls = data.frame(condition = factor(x = rep(x = c(names(groups),
                                                    paste(names(groups), "NULL", sep = ".")),
                                              times = rep(lengths(x = groups), 2)),
                                      levels = c(names(groups),
                                                 paste(names(groups), "NULL", sep = "."))),
                   row.names = colnames(counts[10:ncol(counts)]))
#   The DESeqDataSet object.
DDS = DESeqDataSetFromMatrix(countData = counts[, 10:ncol(counts)],
                             colData = smpls,
                             design = ~ condition)
mcols(DDS) = DataFrame(mcols(DDS), data.frame(gene.id = counts[, 3],
                                              intron.id = counts[, 1]))
##  Differential expression analysis.
DDS = DESeq(object = DDS)
##  Development of FPKM-based expression matrix for introns.
FPKM = data.frame(condition = factor(x = rep(x = names(groups), times = lengths(x = groups)),
                                     levels = names(groups)),
                  row.names = colnames(counts[10:n.cols]))
FPKM = DESeqDataSetFromMatrix(countData = counts[, 10:n.cols],
                              colData = FPKM,
                              design = ~ condition)
mcols(FPKM) = DataFrame(mcols(FPKM), data.frame(basepairs = counts[, 9]))
FPKM = fpkm(object = FPKM, robust = TRUE)
##  Consolidation and filtering of the results for the differential coverage of introns.
res.list = list()
for (i in 1:length(groups)){
res = results(object = DDS,
              contrast = c("condition",
                           names(groups)[i],
                           paste(names(groups)[i], "NULL", sep = ".")))
res = cbind(counts[, 1:9],
            counts[, groups[[i]]],
            rowMeans(counts[, groups[[i]]]),
            FPKM[, groups[[i]]],
            rowMeans(FPKM[, groups[[i]]]),
            res@listData$log2FoldChange,
            res@listData$pvalue,
            res@listData$padj)
colnames(res) = c(colnames(counts[1:9]),
                  paste(groups[[i]], "RAW", sep = "."),
                  "Mean.RAW",
                  paste(groups[[i]], "RPKM", sep = "."),
                  "Mean.FPKM",
                  "logFC",
                  "p.value",
                  "q.value")
res = res[res$Mean.FPKM > quantile(x = res$Mean.FPKM[res$Mean.FPKM > 0], probs = thr.prob)[[1]], ]
res = res[res$logFC > thr.logFC, ]
res = res[res$q.value < thr.q.value, ]
rownames(res) = NULL
res.list[[i]] = res
}
names(res.list) = names(groups)
##  Saving of the final results in TXT files of tab-delimited format.
if (!is.null(f.RIs)){
    for (i in 1:length(res.list)){
         write.table(x = res.list[[i]],
                     file = paste(paste(paste(d.work, f.RIs, sep = "/"),
                                              names(res.list[i]),
                                              sep = ", "), " transcriptome.txt", sep = ""),
                     sep = "\t",
                     quote = FALSE,
                     col.names = TRUE,
                     row.names = FALSE)
        }
}
##  Returning of the final object of class list.
return(res.list)
}
### A simple example of function use.
res = DESeq2RIs(d.work = "D:/Vasily Grinev",
                f.counts = "The NML project, intronic bins, primary matrix of read counts.txt",
                f.exrpGenes = "The NSL project, list of expressed genes.txt",
                rL = 76,
                f.null = "The NML project, intronic bins, count matrix with null distributions.txt",
                groups = list(NC = c("NC_S1", "NC_S2", "NC_S3"),
                              KANSL2 = c("KANSL2_S1", "KANSL2_S2"),
                              KANSL3 = c("KANSL3_S1", "KANSL3_S2", "KANSL3_S3"),
                              KAT8 = c("KAT8_S1", "KAT8_S2", "KAT8_S3")),
                thr.logFC = 1,
                thr.q.value = 0.01,
                thr.prob = 0.97,
                f.RIs = "The NML project, DESeq-based RIs")
