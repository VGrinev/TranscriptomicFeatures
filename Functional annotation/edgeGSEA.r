###################################################################################################
##  High-level R function for collection of the leading edge genes from GSEA output.             ##
##  (c) GNU GPL Vasily V. Grinev, 2017-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work - path to and name of the work directory.
##  html.file1 - the name of GSEA HTML report file with positive values of NES.
##  html.file2 - the name of GSEA HTML report file with negative values of NES.
##  gmt.file - the name of GMT file in GSEA output folder "edb".
##  rnk.file - the name of RNK file in GSEA output folder "edb".
##  FDR - an integer argument. It is a threshold for FDR-adjusted p-value. Default value is 0.05.
#   All files should be in work directory.
edgeGSEA = function(d.work, html.file1, html.file2, gmt.file, rnk.file, FDR = 0.05){
##  Loading of the required library.
#   This code was successfully tested with library XML v.3.98-1.20.
suppressMessages(require(XML))
##  Loading of the GMT file as a list of gene sets.
gene.sets = strsplit(x = readLines(con = paste(d.work, gmt.file, sep = "/")), split = "\t")
names(gene.sets) = sapply(X = gene.sets, "[", 1)
gene.sets = lapply(X = lapply(X = gene.sets, "[", -1:-2), FUN = function(y){y[which(x = y != "")]})
##  Loading of the GSEA HTML report files and collection of leading edge information.
gsea.p = readHTMLTable(doc = paste(d.work, html.file1, sep = "/"), header = FALSE)[[1]]
gsea.p = data.frame(gsea.p, stringsAsFactors = FALSE)[, c(2, 8, 11)]
gsea.p$V2 = as.character(x = gsea.p$V2)
gsea.p$V8 = as.numeric(x = as.character(x = gsea.p$V8))
colnames(gsea.p) = c("id", "fdr", "edge.size")
gsea.p$edge.size = unlist(x = lapply(X = gsea.p[, 3],
                                     FUN = function(y){as.numeric(x = sub("%.+",
                                                                          "",
                                                                          sub("tags=",
                                                                              "",
                                                                              y)))/100}))
gsea.p$effect = "positive"
gsea.n = readHTMLTable(doc = paste(d.work, html.file2, sep = "/"), header = FALSE)[[1]]
gsea.n = data.frame(gsea.n, stringsAsFactors = FALSE)[, c(2, 8, 11)]
gsea.n$V2 = as.character(x = gsea.n$V2)
gsea.n$V8 = as.numeric(x = as.character(x = gsea.n$V8))
colnames(gsea.n) = c("id", "fdr", "edge.size")
gsea.n$edge.size = unlist(x = lapply(X = gsea.n[, 3],
                                     FUN = function(y){as.numeric(x = sub("%.+",
                                                                          "",
                                                                          sub("tags=",
                                                                              "",
                                                                              y)))/100}))
gsea.n$effect = "negative"
edge.sizes = rbind(gsea.n, gsea.p)
edge.sizes = edge.sizes[order(edge.sizes$id), ]
edge.sizes = edge.sizes[edge.sizes$fdr < FDR, ]
rownames(edge.sizes) = NULL
##  Loading of the RNK file as a standard data frame object.
gene.ranks = read.table(file = paste(d.work, rnk.file, sep = "/"),
                        sep = "\t",
                        header = FALSE,
                        quote = "\"",
                        as.is = TRUE)
gene.ranks = gene.ranks[order(-gene.ranks$V2), ]
##  Final sub-setting of genes.
gene.subsets = matrix(nrow = nrow(edge.sizes), ncol = 2)
gene.subsets[, 1] = edge.sizes$id
for (i in 1:nrow(edge.sizes)){
gene.set = gene.ranks[gene.ranks$V1 %in% gene.sets[names(gene.sets) == edge.sizes[i, ]$id][[1]], ]
if (edge.sizes[i, ]$effect == "positive"){
    gene.set = gene.set[1:round(nrow(gene.set) * edge.sizes[i, ]$edge.size), ]$V1
    gene.subsets[i, 2] = paste(sort(gene.set), collapse = ", ")
}else{
    gene.set = rev(gene.set[, 1])[1:round(nrow(gene.set) * edge.sizes[i, ]$edge.size)]
    gene.subsets[i, 2] = paste(sort(gene.set), collapse = ", ")
 }
}
colnames(gene.subsets) = c("id", "gene.subset")
##  Returning of results.
return(gene.subsets)
}
### A simple example of function use.
res = edgeGSEA(d.work = "D:/Vasily Grinev",
               html.file1 = "gsea_report_for_na_neg_1524849911291.html",
               html.file2 = "gsea_report_for_na_pos_1524849911291.html",
               gmt.file = "gene_sets.gmt",
               rnk.file = "NSL_all_proteins.IDs.rnk",
               FDR = 0.05)
