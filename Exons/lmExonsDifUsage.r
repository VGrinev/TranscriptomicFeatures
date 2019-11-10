###################################################################################################
##  A complete R-wrap for identification of differentially used exonic bins                      ##
##  using the functionality of edgeR and limma packages.                                         ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work     - character string giving the path to and name of work directory.
##  f.counts   - character string giving the name of TXT file in tab-delimited format containing
#                summarized RNA-Seq reads for exonic bins. It is typically output file of
#                exonicBinsCounts function.
##  thr.cpm    - numerical threshold value for CPM-based filtration of a count matrix.
#                Default value is 1.
##  ctrl       - character vector giving the list of control samples.
##  exper      - character vector giving the list of experimental samples.
##  f.diffBins - character string giving the name of TXT file in tab-delimited format for storing
#                of final results. Default value is NULL (no any file should be saved).
lmExonsDifUsage = function(d.work, f.counts, thr.cpm = 1, ctrl, exper, f.diffBins = NULL){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with the library edgeR v.3.26.6 and library limma v.3.40.6.
suppressMessages(require(edgeR))
suppressMessages(require(limma))
##  Loading of the count matrix for exonic bins.
counts = read.table(file = paste(d.work, f.counts, sep = "/"),
                    sep = "\t",
                    header = TRUE,
                    quote = "\"",
                    as.is = TRUE)
##  Filtering of the count matrix.
#   Filtering against to low sequencing depth.
cpm.ctrl = counts[rowSums(cpm(counts[, ctrl]) >= thr.cpm) == length(ctrl), ]$exon_id
cpm.exper = counts[rowSums(cpm(counts[, exper]) >= thr.cpm) == length(exper), ]$exon_id
counts = counts[counts$exon_id %in% sort(x = unique(x = c(cpm.ctrl, cpm.exper))), ]
#   Filtering against one-bin genes.
counts = counts[counts$gene_id %in%
                attr((table(counts$gene_id) == 1)
                [(table(counts$gene_id) == 1) == "FALSE"], "dimnames")[[1]], ]
rownames(counts) = NULL
##  Wrapping of the count matrix in a DGEList object.
#   A factors for grouping of samples.
groups = factor(x = c(rep("NC", length(ctrl)), rep("KD", length(exper))), levels = c("NC", "KD"))
#   The DGEList object.
mat = DGEList(counts = counts[, c(ctrl, exper)],
              lib.size = NULL,
              norm.factors = NULL,
              group = groups,
              genes = counts[, c("gene_id", "exon_id")],
              remove.zeros = FALSE)
##  Estimation of normalization (scaling) factors
#   using the "trimmed mean of M-values" normalization method.
mat = calcNormFactors(object = mat, method = "TMM")
##  Performing of voom normalisation and transformation of the count matrix.
#   Calculation of design matrix.
design = model.matrix(~0 + groups)
colnames(design) = levels(groups)
#   Converting of raw counts to log-counts per million with associated precision weights.
mat = voom(counts = mat, design = design, plot = FALSE, save.plot = FALSE)
##  Fitting of linear models to the data.
fit = lmFit(object = mat, design = design, method = "ls")
##  Calculation of contrasts matrix.
contrs = makeContrasts(contr = KD - NC, levels = design)
##  Re-calculation of the linear models.
fit = contrasts.fit(fit = fit, contrasts = contrs)
##  Analysis of the differential usage of exonic bins.
dBin = diffSplice(fit = fit,
                  geneid = "gene_id",
                  exonid = "exon_id",
                  robust = TRUE,
                  verbose = FALSE)
##  Consolidation and saving of final results.
res.table = cbind(dBin$genes[, c(2, 1)],
                  dBin$coefficients,
                  dBin$p.value,
                  p.adjust(p = dBin$p.value, method = "BH", n = length(dBin$p.value)))
colnames(res.table) = c("exon_id", "gene_id", "logFC", "p_value", "q_value")
res.table = res.table[order(res.table$exon_id), ]
if (!is.null(f.diffBins)){
write.table(x = res.table,
            file = paste(d.work, f.diffBins, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
##  Returning of the final object of class data frame.
return(res.table)
}
