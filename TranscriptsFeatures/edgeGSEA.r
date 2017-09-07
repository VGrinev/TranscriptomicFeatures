###################################################################################################
##  High-level R function for collection of the leading edge genes from GSEA output.             ##
##  (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  work_dir - path to and name of the work directory.
##  html_file1 - path to and name of the GSEA HTML report file
##  with positive values of normalized enrichment scores.
##  html_file2 - path to and name of the GSEA HTML report file
##  with negative values of normalized enrichment scores.
##  gmt_file - path to and name of the GMT file in GSEA output folder "edb".
##  rnk_file - path to and name of the RNK file in GSEA output folder "edb".
##  FDR - an integer argument. It is a threshold for FDR-adjusted p-value. Default value is 0.05.
edgeGSEA = function(work_dir, html_file1, html_file2, gmt_file, rnk_file, FDR = 0.05){
##  Loading of the required library.
if (!suppressMessages(require(XML))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("XML")
}
##  Loading of the GMT file as a list of gene sets.
gene.sets = strsplit(readLines(paste(work_dir, gmt_file, sep = "/")), "\t")
names(gene.sets) = sapply(gene.sets, "[", 1)
gene.sets = lapply(lapply(gene.sets, "[", -1:-2), function(x){x[which(x != "")]})
##  Loading of the GSEA HTML report files and collection of leading edge information.
gsea.p = readHTMLTable(doc = paste(work_dir, html_file1, sep = "/"), header = FALSE)[[1]]
gsea.p = data.frame(gsea.p, stringsAsFactors = FALSE)[, c(2, 8, 11)]
gsea.p$V2 = as.character(gsea.p$V2)
gsea.p$V8 = as.numeric(as.character(gsea.p$V8))
colnames(gsea.p) = c("id", "fdr", "edge_size")
gsea.p$edge_size = unlist(lapply(gsea.p[, 3],
                                 function(x){as.numeric(sub("%.+", "", sub("tags=", "", x)))/100}))
gsea.p$effect = "positive"
gsea.n = readHTMLTable(doc = paste(work_dir, html_file2, sep = "/"), header = FALSE)[[1]]
gsea.n = data.frame(gsea.n, stringsAsFactors = FALSE)[, c(2, 8, 11)]
gsea.n$V2 = as.character(gsea.n$V2)
gsea.n$V8 = as.numeric(as.character(gsea.n$V8))
colnames(gsea.n) = c("id", "fdr", "edge_size")
gsea.n$edge_size = unlist(lapply(gsea.n[, 3],
                                 function(x){as.numeric(sub("%.+", "", sub("tags=", "", x)))/100}))
gsea.n$effect = "negative"
edge.sizes = rbind(gsea.n, gsea.p)
edge.sizes = edge.sizes[order(edge.sizes$id), ]
edge.sizes = edge.sizes[edge.sizes$fdr < FDR, ]
rownames(edge.sizes) = NULL
##  Loading of the RNK file as a standard data frame object.
gene.ranks = read.table(file = paste(work_dir, rnk_file, sep = "/"),
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
    gene.set = gene.set[1:round(nrow(gene.set) * edge.sizes[i, ]$edge_size), ]$V1
    gene.subsets[i, 2] = paste(sort(gene.set), collapse = ", ")
}else{
    gene.set = rev(gene.set[, 1])[1:round(nrow(gene.set) * edge.sizes[i, ]$edge_size)]
    gene.subsets[i, 2] = paste(sort(gene.set), collapse = ", ")
}
}
##  Returning of results.
colnames(gene.subsets) = c("id", "gene.subset")
return(gene.subsets)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/edgegsea.r")
#   res = edgeGSEA(work_dir = "/home/work",
#                  html_file1 = "gsea_report_for_na_pos_1500275393272.html",
#                  html_file2 = "gsea_report_for_na_neg_1500275393272.html",
#                  gmt_file = "gene_sets.gmt",
#                  rnk_file = "ranked_genes.rnk",
#                  FDR = 0.05)
