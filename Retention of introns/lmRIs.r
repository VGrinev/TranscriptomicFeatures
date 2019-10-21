###################################################################################################
##  A complete R-wrap for identification                                                         ##
##  of retained introns in RNA-Seq data with the edgeR and limma packages.                       ##
##  (c) GNU GPL Vasily V. Grinev, 2017-2019. grinev_vv[at]bsu.by                                 ##
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
lmRIs = function(d.work, f.counts, f.exrpGenes, rL, f.null = NULL, groups,
                 thr.prob = 0.97, thr.logFC = 1, thr.q.value = 0.01, f.RIs = NULL){
##  Loading of a required auxiliary library.
#   This code was successfully tested with the library edgeR v.3.26.6 and library limma v.3.40.6.
suppressMessages(require(edgeR))
suppressMessages(require(limma))
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
counts = counts[counts$effective_length >= (rL + 1), ]
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
                                                               FUN = function(y){y * iL/sum(iL)}))
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
##  Wrapping of the count matrix in a DGEList object.
mat = DGEList(counts = counts[, 10:ncol(counts)],
              lib.size = NULL,
              norm.factors = NULL,
              group = factor(x = rep(x = c(names(groups), paste(names(groups), "NULL", sep = ".")),
                                     times = rep(lengths(x = groups), 2)),
                             levels = c(names(groups), paste(names(groups), "NULL", sep = "."))),
              genes = counts[, "intron_id"])
##  Estimation of normalization (scaling) factors
#   using the "trimmed mean of M-values" normalization method.
mat = calcNormFactors(object = mat, method = "TMM")
##  Performing of voom normalisation and transformation of the count matrix.
#   Calculation of design matrix.
design = model.matrix(~0 + factor(x = rep(x = c(names(groups),
                                                paste(names(groups), "NULL", sep = ".")),
                                          times = rep(lengths(x = groups), 2)),
                                  levels = c(names(groups),
                                             paste(names(groups), "NULL", sep = "."))))
colnames(design) = c(names(groups), paste(names(groups), "NULL", sep = "."))
#   Converting of raw counts to log-counts per million with associated precision weights.
mat = voom(counts = mat, design = design, plot = FALSE, save.plot = FALSE)
##  Fitting of linear models to the data.
fit = lmFit(object = mat, design = design, method = "ls")
##  Test for difference between empirical and null distributions.
#   Calculation of contrasts matrix.
fit = contrasts.fit(fit = fit,
                    contrasts = makeContrasts(contrasts = as.list(paste(names(groups),
                                                                        paste(names(groups),
                                                                              "NULL",
                                                                              sep = "."),
                                                                        sep = " - ")),
                                              levels = c(names(groups),
                                                         paste(names(groups), "NULL", sep = "."))))
#   Calculation of empirical Bayes statistics for a differential coverage of introns.
fit = eBayes(fit = fit, robust = TRUE)
##  Development of RPKM-based expression matrix for introns.
RPKM = DGEList(counts = counts[10:n.cols],
               lib.size = calcNormFactors(object = counts[10:n.cols], method = "TMM") *
                          colSums(counts[10:n.cols]))
RPKM$genes$length = counts$effective_length
RPKM = data.frame(rpkm(RPKM))
##  Consolidation and filtering of the results for the differential coverage of introns.
res.list = list()
for (i in 1:length(groups)){
res = cbind(counts[, 1:9],
            counts[, groups[[i]]],
            rowMeans(counts[, groups[[i]]]),
            RPKM[, groups[[i]]],
            rowMeans(RPKM[, groups[[i]]]),
            fit$coefficients[, paste(paste(names(groups)[i],
                                           names(groups)[i],
                                           sep = " - "),
                                     "NULL",
                                     sep = ".")],
            fit$p.value[, paste(paste(names(groups)[i],
                                      names(groups)[i],
                                      sep = " - "),
                                "NULL",
                                sep = ".")])
colnames(res) = c(colnames(counts[1:9]),
                  paste(groups[[i]], "RAW", sep = "."),
                  "Mean.RAW",
                  paste(groups[[i]], "RPKM", sep = "."),
                  "Mean.RPKM",
                  "logFC",
                  "p.value")
res$q.value = p.adjust(p = res$p.value, method = "BH", n = length(res$p.value))
res = res[res$Mean.RPKM > quantile(x = res$Mean.RPKM[res$Mean.RPKM > 0], probs = thr.prob)[[1]], ]
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
