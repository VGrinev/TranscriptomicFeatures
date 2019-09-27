###################################################################################################
##  High-level R function for identification of differentially used TSSs                         ##
##  in experimental transcriptomes.                                                              ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  in.gtf - path to folder and name of the file in GTF/GFF format with assembled transcriptome.
#   It is typically output file of Cufflinks.
##  tss.gr - path to folder and name of the file in read_group_tracking format containing summed
#   expression and counts of transcripts sharing each tss_id in each replicate. It is typically
#   output file of Cuffdiff.
##  f.tss - the name of output TXT file in tab-delimited format for storing a FPKM matrix of TSSs.
#   Default value is NULL.
##  f.difUsage - the name of output TXT file in tab-delimited format for storing the final results
#   of differential test. Default value is NULL.
diffTSSsUsage = function(in.gtf, tss.gr, f.tss = NULL, f.difUsage = NULL){
##  Loading of required auxiliary library.
#   This code was successfully tested with libraries rtracklayer v.1.44.3,
#   edgeR v.3.26.6 and limma v.3.40.6.
suppressMessages(require(rtracklayer))
suppressMessages(require(limma))
suppressMessages(require(edgeR))
##  Loading and sub-setting of input data.
in.gtf = import(con = in.gtf, format = "gtf")
in.gtf = in.gtf[, c("tss_id", "transcript_id", "gene_name")]
tr.list = makeGRangesListFromDataFrame(df = as.data.frame(x = in.gtf),
                                       split.field = "transcript_id",
                                       keep.extra.columns = TRUE)
##  Retrieving a list of TSSs and associated annotations.
TSSs = unlist(x = as(lapply(X = tr.list,
                            FUN = function(x){if (as.character(strand(x)@values) == "+"){
                                                  x = x[1, ]
                                                  end(x) = start(x)
                                             }else{
                                                  x = x[length(x), ]
                                                  start(x) = end(x)
                                             }
                                             return(x)}), "GRangesList"))
TSSs = as.data.frame(x = TSSs)
TSSs = TSSs[complete.cases(TSSs), ]
TSSs.str = aggregate(x = unique(TSSs[, c(1, 5:6)])[, 1:2],
                     by = list(unique(TSSs[, c(1, 5:6)])[, 3]),
                     FUN = function(x){paste(x, collapse = ", ")})
TSSs = cbind(TSSs.str[, 1],
             aggregate(x = unique(TSSs[, 6:7])[, 2],
                       by = list(unique(TSSs[, 6:7])[, 1]),
                       FUN = function(x){paste(x, collapse = ", ")})[, 2],
             TSSs.str[, 2],
             aggregate(x = TSSs[, 2], by = list(TSSs[, 6]), FUN = min)[, 2],
             aggregate(x = TSSs[, 3], by = list(TSSs[, 6]), FUN = max)[, 2],
             TSSs.str[, 3])
colnames(TSSs) = c("tss_id", "gene_name", "seqnames", "start", "end", "strand")
##  Development a FPKM matrix of TSSs.
TSSs.fpkm = read.table(file = tss.gr,
                       sep = "\t",
                       header = TRUE,
                       quote = "\"",
                       as.is = TRUE)[, c(1:3, 7)]
#   A factors for grouping of samples.
factors = unique(TSSs.fpkm[, c("condition", "replicate")])
factors = factor(factors[order(factors$condition), ]$condition)
TSSs.fpkm = lapply(X = split(x = TSSs.fpkm, f = TSSs.fpkm$condition),
                   FUN = function(y){split(x = y, f = y$replicate)})
TSSs.fpkm = data.frame(cbind(TSSs.fpkm[[1]][[1]][, 1],
                             do.call(what = cbind,
                                     args = unlist(x = lapply(X = TSSs.fpkm,
                                                              FUN = function(y){lapply(X = y,
                                                                  FUN = function(z){z = z[, 4]})}),
                                                   recursive = FALSE))),
                       stringsAsFactors = FALSE)
TSSs.fpkm = TSSs.fpkm[order(TSSs.fpkm$V1), ]
TSSs = TSSs[TSSs[, 1] %in% TSSs.fpkm$V1, ]
TSSs.fpkm = TSSs.fpkm[TSSs.fpkm$V1 %in% TSSs[, 1], ]
options(stringsAsFactors = FALSE)
TSSs = cbind(TSSs, TSSs.fpkm[, -1])
TSSs[, c(4, 5, 7:length(TSSs))] = apply(X = TSSs[, c(4, 5, 7:length(TSSs))],
                                        MARGIN = 2,
                                        FUN = as.numeric)
TSSs = TSSs[rowSums(TSSs[, 7:ncol(TSSs)]) > 0, ]
if (!is.null(f.tss)){
write.table(TSSs,
            file = f.tss,
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
##  Wrapping of the FPKM matrix in a DGEList object.
mat.FPKM = DGEList(counts = TSSs[, c(7:ncol(TSSs))],
                   lib.size = NULL,
                   norm.factors = NULL,
                   group = factors,
                   genes = TSSs[, c("gene_name", "tss_id")],
                   remove.zeros = FALSE)
##  Performing of the voom normalisation and transformation of the FPKM data.
design = model.matrix(~0 + factors)
colnames(design) = levels(factors)
mat.FPKM = voom(counts = mat.FPKM,
                design = design,
                lib.size = rep(x = 1e6, times = nrow(x = mat.FPKM$samples)))
##  Fitting of the linear model to the transformed counts.
lin_fit = lmFit(object = mat.FPKM, design = design, method = "ls")
#   Calculation of the contrast matrix.
contrast.matrix = makeContrasts(cntr = paste(attr(design, "dimnames")[[2]], collapse = " - "),
                                levels = design)
lin_fit = contrasts.fit(fit = lin_fit, contrasts = contrast.matrix)
##  Analysis of the differential usage of TSSs.
difUsage = diffSplice(fit = lin_fit,
                      geneid = "gene_name",
                      exonid = "tss_id",
                      robust = TRUE,
                      verbose = FALSE)
##  Consolidation and saving of the results of interest.
difUsage = data.frame(cbind(difUsage$genes,
                                difUsage$coefficients[, 1],
                                difUsage$p.value[, 1],
                                p.adjust(p = difUsage$p.value[, 1],
                                         method = "BH",
                                         n = length(x = difUsage$p.value[, 1])),
                           -log10(p.adjust(p = difUsage$p.value[, 1],
                                         method = "BH",
                                         n = length(x = difUsage$p.value[, 1])))),
                     stringsAsFactors = FALSE)
colnames(difUsage) = c("gene_name", "tss_id", "logFC", "p_value", "q_value", "neg_log10_q_value")
difUsage = difUsage[, c(2, 1, 3:6)]
if (!is.null(f.difUsage)){
write.table(difUsage,
            file = f.difUsage,
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
return(list(TSSs = TSSs, difUsage = difUsage))
}
