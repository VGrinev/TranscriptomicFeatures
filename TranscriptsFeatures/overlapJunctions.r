##############################################################################################################
## High-level R function for calculation of overlaps between                                                ##
## reference and experimentally detected exon-exon junctions.                                               ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
##  models - path to folder and name of the SQLite database with reference transcriptional models of genes.
##  experimental - path to folder and name of the txt file in tab-delimited format with experimentally
#   detected exon-exon junctions. This is input data file which should include seven mandatory fields:
#   i) Event_ID (ID of the splicing event);
#   ii) seqnames (name of chromosome or scaffold with prefix "chr");
#   iii) start (start genomic coordinate of the exon-exon junction);
#   iv) end (end genomic coordinate of the exon-exon junction);
#   v) strand (strand information about exon-exon junction).
#   vi) logFC (log2 fold change in usage of the splicing event between two conditions);
#   vii) q_value (FDR-adjusted p-value for logFC).
#   Additionally, this file should include at least one column with normalized counts that support exon-exon
#   junctions detected in control as well as in experimental cells.
##  log.FC - an integer argument. It is a threshold for log2 fold change. Default value is 1.
##  q.value - an integer argument. It is a threshold for FDR-adjusted p-value. Default value is 0.05.
##  ctrl.samples - a character vector with names of control samples to be analysed. These are names of columns
#   in input data file with normalized counts that support exon-exon junctions in control cells.
##  expr.samples - a character vector with names of experimental samples to be analysed. These are names of
#   columns in input data file with normalized counts that support exon-exon junctions in experimental cells.
##  thr - a numeric vector with thresholds for filtration of experimental data against normalized counts
#   supporting exon-exon junctions. End user should indicate the starting value of the sequence (default value
#   is 0.03125 CPM), step factor (default value is *2), and number of filtration steps (default value is 18).
overlapJunctions = function(models, experimental, log.FC = 1, q.value = 0.05, ctrl.samples, expr.samples, thr = c(0.03125, 2, 18)){
## Loading of required auxiliary libraries.
#  This code was tested only with v.0.51.0 of matrixStats library and v.1.26.3 of GenomicFeatures library.
cat("Loading of required auxiliary libraries.\n")
flush.console()
if (!suppressMessages(require(GenomicFeatures))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("GenomicFeatures")
}
if (!suppressMessages(require(matrixStats))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("matrixStats")
}
## Loading of reference transcriptional models of genes.
#  This code was tested only with Ensembl models of human genes (Ensembl release 85, GRCh38.p7 human
#  reference genome, July 2016).
cat("Loading of reference transcriptional models of genes.\n")
flush.console()
tx = loadDb(models)
## Calculation of reference exon-exon junctions.
#  This step produces an GRanges object which contain only exon-exon junctions (without
#  any gene and/or transcript annotations) for genes belong to standard human chromosomes.
cat("Calculation of reference exon-exon junctions.\n")
flush.console()
tx = exonsBy(tx, "tx", use.names = FALSE)
modelEEJ = setdiff(range(tx), tx)
modelEEJ = setNames(unlist(modelEEJ), NULL) + 1
modelEEJ = as.data.frame(modelEEJ)
modelEEJ$seqnames = paste("chr", modelEEJ$seqnames, sep = "")
modelEEJ = modelEEJ[modelEEJ$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
modelEEJ = modelEEJ[!duplicated(modelEEJ), ]
rownames(modelEEJ) = NULL
modelEEJ = makeGRangesFromDataFrame(modelEEJ, keep.extra.columns = TRUE)
modelEEJ = sort(modelEEJ)
## Loading and sub-setting of experimentally detected exon-exon junctions.
#  This step produces three GRanges objects:
#  i) differential exon-exon junctions with prevalence in control cells;
#  ii) differential exon-exon junctions with prevalence in experimental cells;
#  iii) non-differential exon-exon junctions.
#  Each GRanges object includes all annotations from input data file.
cat("Loading and sub-setting of experimentally detected exon-exon junctions.\n")
flush.console()
expEEJ = read.table(file = experimental, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
expEEJ = makeGRangesFromDataFrame(expEEJ, keep.extra.columns = TRUE)
expEEJ.dif = expEEJ[abs(expEEJ$logFC) >= log.FC & expEEJ$q_value < q.value, ]
expEEJ.dif.minus = expEEJ.dif[expEEJ.dif$logFC < -log.FC, ]
expEEJ.dif.plus = expEEJ.dif[expEEJ.dif$logFC > log.FC, ]
expEEJ.non.dif = expEEJ[!expEEJ$Event_ID %in% expEEJ.dif$Event_ID, ]
subset.logFC = list(expEEJ.dif.plus, expEEJ.dif.minus, expEEJ.non.dif)
## Detection of identical exon-exon junctions in experimental versus reference data.
cat("Detection of identical exon-exon junctions in experimental versus reference data.\n")
flush.console()
model.overlaps = matrix(ncol = 3, nrow = 3)
model.overlaps[, 1] = c("differential splicing, logFC > 1", "differential splicing, logFC < -1", "non-differential splicing")
for (i in 1:length(subset.logFC)){
    hits = findOverlaps(subset.logFC[[i]], modelEEJ, type = c("equal"), ignore.strand = FALSE)
    model.overlaps[i, 2:3] = round(c(length(hits)/hits@nLnode, 1 - length(hits)/hits@nLnode), digits = 3)
}
model.overlaps = data.frame(model.overlaps, stringsAsFactors = FALSE)
colnames(model.overlaps) = c("sub-set of junctions", "annotated", "non-annotated")
## Detection of identical exon-exon junctions in experimental versus control cells. 
#  This code compares differential and non-differential exon-exon junctions against full list of exon-exon
#  junctions detected in control or experimental cells.
cat("Detection of identical exon-exon junctions in experimental versus control cells.\n")
flush.console()
degOverlap.list = list()
samples = c(ctrl.samples, expr.samples)
thr = thr[1] * thr[2]^seq(0, thr[3] - 1)
for (j in 1:length(samples)){
    degOverlap.matrix = matrix(nrow = length(thr), ncol = length(subset.logFC) + 1)
    degOverlap.matrix[, 1] = thr
    for (k in 1:length(thr)){
        subset.expEEJ = expEEJ[expEEJ[, samples[j]]@elementMetadata@listData[[1]] >= thr[k], ]
        for (l in 1:length(subset.logFC)){
            hits = findOverlaps(subset.logFC[[l]], subset.expEEJ, type = c("equal"), ignore.strand = FALSE)
            degOverlap.matrix[k, l + 1] = length(hits)/hits@nLnode
        }
    }
    degOverlap.list[[j]] = degOverlap.matrix
}
names(degOverlap.list) = samples
options(scipen = 999)
expr.overlaps.control = matrix(ncol = length(subset.logFC) * 2 + 1, nrow = length(thr))
expr.overlaps.control[, 1] = thr
x = degOverlap.list[ctrl.samples]
m = 1
for (o in 1:length(subset.logFC)){
    expr.overlaps.control[, m + 1] = round(rowMeans(do.call(cbind, lapply(x, function(x) {x[, o + 1]}))), digits = 4)
    expr.overlaps.control[, m + 2] = round(rowSds(do.call(cbind, lapply(x, function(x) {x[, o + 1]}))), digits = 5)
    m = m + 2
}
expr.overlaps.control = data.frame(expr.overlaps.control, stringsAsFactors = FALSE)
colnames(expr.overlaps.control) = c("CPM_threshold",
                                    "difspl_logFC.more.1_mean",
                                    "difspl_logFC.more.1_sd",
                                    "difspl_logFC.less.minus.1_mean",
                                    "difspl_logFC.less.minus.1_sd",
                                    "non.difspl_mean",
                                    "non.difspl_sd")
expr.overlaps.experimental = matrix(ncol = length(subset.logFC) * 2 + 1, nrow = length(thr))
expr.overlaps.experimental[, 1] = thr
x = degOverlap.list[expr.samples]
n = 1
for (p in 1:length(subset.logFC)){
    expr.overlaps.experimental[, n + 1] = round(rowMeans(do.call(cbind, lapply(x, function(x) {x[, p + 1]}))), digits = 4)
    expr.overlaps.experimental[, n + 2] = round(rowSds(do.call(cbind, lapply(x, function(x) {x[, p + 1]}))), digits = 5)
    n = n + 2
}
expr.overlaps.experimental = data.frame(expr.overlaps.experimental, stringsAsFactors = FALSE)
colnames(expr.overlaps.experimental) = c("CPM_threshold",
                                         "difspl_logFC.more.1_mean",
                                         "difspl_logFC.more.1_sd",
                                         "difspl_logFC.less.minus.1_mean",
                                         "difspl_logFC.less.minus.1_sd",
                                         "non.difspl_mean",
                                         "non.difspl_sd")
## Consolidation results.
#  Function returns an object of class list which include all important results. 
cat("Calculations have been finished successfully.\n")
flush.console()
return(list(EEJ.models = modelEEJ,
            Models.overlaps = model.overlaps,
            Expr.overlaps.control = expr.overlaps.control,
            Expr.overlaps.experimental = expr.overlaps.experimental))
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/overlapjunctions.r")
#   modelData = "GVV/Ensembl_GRCh38.p7_release.85.sqlite"
#   exprData = "GVV/Experimental_junctions.txt"
#   ctrl = c("ctrl_1", "ctrl_2", "ctrl_3")
#   expr = c("expr_1", "expr_2", "expr_3")
#   res = overlapJunctions(models = modelData, experimental = exprData, ctrl.samples = ctrl, expr.samples = expr)
