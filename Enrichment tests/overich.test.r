###################################################################################################
##  User-defined non-ranked gene enrichment test.                                                ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  query - path to folder and name of the TXT file in tab-delimited format with a list of query
#   gene sets to be analysed. This file should include two mandatory fields:
#   i) path (path to folder with file of query gene set);
#   ii) file (name of file with query gene set).
#   A file with query gene set should be standard TXT file in tab-delimited format with one unnamed
#   column. This column must list all members of query gene set. It is strongly recommended to name
#   the file by the name of query gene set.
##  subject - path to folder and name of the TXT file in tab-delimited format with a list of
#   subject gene sets to be analysed. This file should include two mandatory fields:
#   i) path (path to folder with file of subject gene set);
#   ii) file (name of file with subject gene set).
#   A file with subject gene set should be standard TXT file in tab-delimited format with one
#   unnamed column. This column must list all members of subject gene set. It is strongly
#   recommended to name the file by the name of subject gene set.
##  subject.size - a minimal size of gene set to be included in statistical test.
#   Default value is 10.
##  ref.size - a size of reference gene set. Default value is NULL. In that case an actual size of
#   reference gene set will be calculated automatically based on provided subject gene sets.
##  FDR - a treshold for the controled false discovery rate. Default value is 0.05.
overich.test = function(query, subject, subject.size = 10, ref.size = NULL, FDR = 0.05){
##  Preparing of query gene set by standard object of class vector.
queries = read.table(file = query, sep = "\t", quote = "\"", header = FALSE, as.is = TRUE)$V1
##  Preparing of subject gene sets by standard object of class list.
subjects = read.table(file = subject, sep = "\t", quote = "\"", header = TRUE, as.is = TRUE)
subjects = setNames(object = split(x = subjects, f = seq(nrow(x = subjects))), nm = subjects$id)
L = as.numeric(x = unlist(x = lapply(X = subjects,
                                     FUN = function(y){length(x = strsplit(x = y$gene.set,
                                                                           split = ", ")[[1]])})))
subjects = subjects[L >= subject.size]
##  Preparing of reference gene set.
ref.genes = unlist(x = lapply(X = subjects,
                              FUN = function(y){strsplit(x = y$gene.set,
                                                         split = ", ")[[1]]}))
ref.genes = sort(x = unique(x = ref.genes))
##  Defining the length of query gene set (Gq).
Gq = queries[queries %in% ref.genes]
Gq = length(x = Gq)
##  Defining the length of reference gene set (Gr).
if (is.null(ref.size)){
    Gr = length(x = ref.genes)
}else{
    Gr = ref.size
}
##  Calculation of fold enrichment and concomitant statistics.
enr = lapply(X = subjects,
             FUN = function(y){gene.set = strsplit(x = y$gene.set, split = ", ")[[1]]
#   Defining the length of subject gene set (Gs).
                               Gs = length(x = gene.set)
#   Preparing a subset of query genes that overlap with a subject gene set.
                               query.overlap = queries[queries %in% gene.set]
#   Defining the size of overlap of query in subject (Gqs).
                               Gqs = length(x = query.overlap)
#   Calculation of expected value (E-value, or E).
                               E = Gq * Gs/Gr
#   Calculation of fold enrichment (FE).
                               FE = round(x = Gqs/E, digits = 2)
#   Calculation of odds ratio and Fisher's exact test.
                               fisher = fisher.test(x = matrix(data = c(Gqs,
                                                                        Gs - Gqs,
                                                                        Gq - Gqs,
                                                                        Gr - Gs - Gq + Gqs),
                                                               nrow = 2))
                               odds = round(x = as.numeric(x = fisher$estimate), digits = 2)
                               p.value = fisher$p.value
                               res = cbind(y$id,
                                           y$name,
                                           Gqs,
                                           Gs,
                                           FE,
                                           odds,
                                           p.value,
                                           0,
                                           paste(sort(x = query.overlap), collapse = ", "))
                               return(res)
                               }
            )
enr = as.data.frame(x = do.call(what = rbind, args = enr), stringsAsFactors = FALSE)
colnames(enr) = c("gene.set.id",
                  "gene.set.name",
                  "query.size",
                  "subject.size",
                  "fold.enrichment",
                  "odds.ratio",
                  "p_value",
                  "FDR.adjusted.p_value",
                  "query.genes.found")
enr[, 3:8] = apply(X = enr[, 3:8], MARGIN = 2, FUN = as.numeric)
enr = enr[enr$query.size >= 3, ]
enr = enr[enr$fold.enrichment > 1, ]
##  Adjustment of the p-values with Benjamini and Hochberg's method.
enr$FDR.adjusted.p_value = p.adjust(p = enr$p_value, method = "BH", n = nrow(enr))
##  Final consolidation of results.
enr = enr[enr$FDR.adjusted.p_value < FDR, ]
##  Returning of results.
return(enr)
}
