###################################################################################################
##  R-wrap for filtration and re-formatting of Cuffdiff-based TSSs differential expression data. ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work - path to and name of work directory.
##  in.file - Cuffdiff tss_group_exp.diff file containing TSSs differential expression data.
##  filtr - a logical argument for filtration against unconventional TSSs data.
#   Default value is FALSE.
##  anno - the name of TXT file in tab-delimited format with Ensembl-based annotations of genes.
##  f.diExpr - the name of output TXT file in tab-delimited format for storing the final results.
#   Default value is NULL.
diffTSSsExpr = function(d.work, in.file, filtr = FALSE, anno, f.diExpr = NULL){
##  Loading and sub-setting of input data.
TSSs = read.table(file = paste(d.work, in.file, sep = "/"),
                  sep = "\t",
                  header = TRUE,
                  quote = "\"",
                  as.is = TRUE)[, c(1, 3, 10, 12:13)]
TSSs$n.log = -log10(TSSs$q_value)
colnames(TSSs) = c("tss_id", "gene_name", "logFC", "p_value", "q_value", "neg_log10_q_value")
##  Filtering of input data.
TSSs = TSSs[!TSSs$gene_name == "-", ]
TSSs = TSSs[!is.infinite(TSSs$logFC), ]
if (filtr == "TRUE"){
fltr = TSSs[gsub(pattern = ".+,.+", replacement = ",", x = TSSs$gene_name) == ",", ]
TSSs = TSSs[!gsub(pattern = ".+,.+", replacement = ",", x = TSSs$gene_name) == ",", ]
fltr.list = list()
for (i in 1:nrow(fltr)){
     fltr.list[[i]] = suppressWarnings(cbind(fltr[i, ],
                                             strsplit(x = fltr[i, ]$gene_name, split = ",")[[1]]))
}
fltr = do.call(what = rbind, args = fltr.list)
fltr$gene_name = as.character(fltr[, ncol(fltr)])
fltr = fltr[, -ncol(fltr)]
TSSs = rbind(TSSs, fltr)
}
TSSs = TSSs[order(TSSs$gene_name), ]
if (!is.null(f.diExpr)){
write.table(TSSs,
            file = paste(d.work, f.diExpr, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
return(TSSs)
}
