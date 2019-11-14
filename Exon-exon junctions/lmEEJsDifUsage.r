###################################################################################################
##  A complete R-wrap for identification of differentially used exon-exon junctions              ##
##  using the functionality of edgeR and limma packages.                                         ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work     - character string giving the path to and name of work directory.
##  f.counts   - character string giving the name of TXT file in tab-delimited format containing
#                summarized RNA-Seq reads spanning exon-exon junctions. It is typically final
#                output file of functions alignSubjunc and matrixEEJs.
##  thr.cpm    - numerical threshold value for CPM-based filtration of a count matrix.
#                Default value is 1.
##  ctrl       - character vector giving the list of control samples.
##  exper      - character vector giving the list of experimental samples.
##  f.diffEEJs - character string giving the name of TXT file in tab-delimited format for storing
#                of final results. Default value is NULL (no any file should be saved).
lmEEJsDifUsage = function(d.work, f.counts, thr.cpm = 1, ctrl, exper, f.diffEEJs = NULL){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with the library edgeR v.3.26.6 and library limma v.3.40.6.
suppressMessages(require(edgeR))
suppressMessages(require(limma))
##  Loading of the count matrix for exon-exon junctions.
counts = read.table(file = paste(d.work, f.counts, sep = "/"),
                    sep = "\t",
                    header = TRUE,
                    quote = "\"",
                    as.is = TRUE)
##  Filtering of the count matrix.
#   It is assumed that the count matrix already does not contain uninformative data (for instance,
#   reads aligned with too low quality or aligned not uniquely were not included in counts matrix).
#   Filtering against to low sequencing depth.
cpm.ctrl = counts[rowSums(cpm(counts[, ctrl]) >= thr.cpm) == length(ctrl), ]$eej_id
cpm.exper = counts[rowSums(cpm(counts[, exper]) >= thr.cpm) == length(exper), ]$eej_id
counts = counts[counts$eej_id %in% sort(x = unique(x = c(cpm.ctrl, cpm.exper))), ]
#   Filtering against to short splicing distances.
counts = counts[counts$width >= 50, ]
#   Filtering against splicing one-event genes.
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
              genes = counts[, c("gene_id", "eej_id")],
              remove.zeros = FALSE)
##  Estimation of normalization (scaling) factors
#   using the "trimmed mean of M-values" normalization method.
#   This step permits to calculate an effective size for each RNA-Seq library.
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
##  Analysis of the differential usage of exon-exon junctions by diffSplice function.
dEEJs = diffSplice(fit = fit,
                   geneid = "gene_id",
                   exonid = "eej_id",
                   robust = TRUE,
                   verbose = FALSE)
##  Consolidation and saving of final results.
res.table = cbind(dEEJs$genes[, c(2, 1)],
                  dEEJs$coefficients,
                  dEEJs$p.value,
                  p.adjust(p = dEEJs$p.value, method = "BH", n = length(dEEJs$p.value)))
colnames(res.table) = c("eej_id", "gene_id", "logFC", "p_value", "q_value")
res.table = res.table[order(res.table$eej_id), ]
if (!is.null(f.diffEEJs)){
write.table(x = res.table,
            file = paste(d.work, f.diffEEJs, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
##  Returning of the final object of class data frame.
return(res.table)
}
