##############################################################################################################
## User-defined non-ranked gene enrichment test.                                                            ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
##  query - path to folder and name of the TXT file in tab-delimited format with a list of query gene sets to
#   be analysed. This file should include two mandatory fields:
#   i) path (path to folder with file of query gene set);
#   ii) file (name of file with query gene set).
#   A file with query gene set should be standard TXT file in tab-delimited format with one unnamed column.
#   This column must list all members of query gene set. It is strongly recommended to name the file by
#   the name of query gene set.
##  subject - path to folder and name of the TXT file in tab-delimited format with a list of subject gene sets
#   to be analysed. This file should include two mandatory fields:
#   i) path (path to folder with file of subject gene set);
#   ii) file (name of file with subject gene set).
#   A file with subject gene set should be standard TXT file in tab-delimited format with one unnamed column.
#   This column must list all members of subject gene set. It is strongly recommended to name the file by
#   the name of subject gene set.
##  ref.size - a size of reference gene set. Default value is NULL. In that case an actual size of reference
#   gene set will be calculated automatically based on provided subject gene sets.
##  p.original - the type of statistical approach that will be used to assess the significance of
#   fold enrichment. The following values are allowed:
#   i) "fisher" (Fisher's exact test);
#   ii) "chisq" (Pearson's chi-squared test);
#   iii) "binom" (binomial exact test);
#   iv) "hyper" (hypergeometric test);
#   v) "random" (random sampling test).
##  ntry - number of random sampling of reference gene set. It is valid only with p.original = "random".
#   Default value is 1000.
##  p.adj - method for adjust the p-values for multiple comparisons. The following values are allowed:
#   "BH" - adjustment of the p-values with Benjamini and Hochberg's method
#   for control of the false discovery rate;
#   "holm" - adjustment of the p-values with Holm's method for control of the family-wise error rate.
enrichmentTest = function(query, subject, ref.size = NULL, p.original = "fisher", ntry = 1000, p.adj = "BH"){
## Preparing of query gene set(-s) by standard object of class list.
query.files = read.table(file = query, sep = "\t", quote = "\"", header = TRUE, as.is = TRUE)
query.names = gsub(".txt", "", query.files$file)
query.files = paste(query.files$path, query.files$file, sep = "/")
query.list = lapply(query.files, function(x){read.table(x, quote = "\"", as.is = TRUE)$V1})
## Preparing of subject gene set(-s) by standard object of class list.
subject.files = read.table(file = subject, sep = "\t", quote = "\"", header = TRUE, as.is = TRUE)
subject.names = gsub(".txt", "", subject.files$file)
subject.files = paste(subject.files$path, subject.files$file, sep = "/")
subject.list = lapply(subject.files, function(x){read.table(x, quote = "\"", as.is = TRUE)$V1})
## Preparing of reference gene set.
reference.list = unlist(subject.list)
## Defining the length of reference gene list (Gr).
if (is.null(ref.size)){
    Gr = length(unique(do.call(c, subject.list)))
}else{
    Gr = ref.size
}
## Calculation of fold enrichment and concomitant statistics.
#  Empty matrices for future data collection.
enrichment = matrix(ncol = length(query.list), nrow = length(subject.list))
odds = matrix(ncol = length(query.list), nrow = length(subject.list))
test = matrix(ncol = length(query.list), nrow = length(subject.list))
#  Query loop.
for (i in 1:length(query.list)){
     # Defining the length of query gene list (Gq).
     Gq = length(query.list[[i]])
     # Subject loop.
     for (j in 1:length(subject.list)){
          # Defining the length of subject gene list (Gs).
          Gs = length(subject.list[[j]])
          # Defining the size of overlap of query in subject (Gqs).
          Gqs = sum(query.list[[i]] %in% subject.list[[j]])
          # Calculation of E-value.
          E = Gq * Gs/Gr
          # Calculation of fold enrichment.
          enrichment[j, i] = round(Gqs/E, digits = 2)
          # Fisher's exact test and calculation of odds ratio.
          fisher.test = fisher.test(matrix(c(Gqs, Gs - Gqs, Gq - Gqs, Gr - Gs - Gq + Gqs), nrow = 2))
          odds[j, i] = round(as.numeric(fisher.test$estimate), digits = 2)
          if (p.original == "fisher"){
              test[j, i] = fisher.test$p.value
          }
          # Pearson's chi-squared test.
          if (p.original == "chisq"){
              chisq.test = suppressWarnings(chisq.test(matrix(c(Gqs,
                                                                Gs - Gqs,
                                                                Gq - Gqs,
                                                                Gr - Gs - Gq + Gqs), nrow = 2)))
              test[j, i] = chisq.test$p.value
          }
          # Binomial exact test.
          if (p.original == "binom"){
              binom.test = binom.test(x = c(Gqs, Gq - Gqs), p = Gs/Gr, alternative = "two.sided")
              test[j, i] = binom.test$p.value
          }
          # Hypergeometric test.
          if (p.original == "hyper"){
              test[j, i] = phyper(q = Gqs - 1, m = Gs, n = Gr - Gs, k = Gq, lower.tail = FALSE)
          }
          # Random test.
          if (p.original == "random"){
              random.FE = double(length(ntry))
              for (k in 1:ntry){
                   random.Gq = sample(reference.list, Gq)
                   random.Gqs = sum(random.Gq %in% subject.list[[j]])
                   random.FE[k] = round((random.Gqs * Gr)/(Gq * Gs), digits = 2)
              }
              if (enrichment[j, i] > 1){
                  p = sum(random.FE >= enrichment[j, i])/ntry
              }else{
                  p = sum(random.FE <= enrichment[j, i])/ntry
              }
          test[j, i] = p
          }
     }
}
## Adjusting the p-values for multiple comparisons.
#  Adjustment of the p-values with Benjamini and Hochberg's method.
if (p.adj == "BH"){
    p.adjusted = apply(test, 2, function(x){p.adjust(p = x, method = "BH", n = length(subject.list))})
}
#  Adjustment of the p-values with Holm's method.
if (p.adj == "holm"){
    p.adjusted = apply(test, 2, function(x){p.adjust(p = x, method = "holm", n = length(subject.list))})
}
## Consolidation of results.
test = format(test, digits = 3, scientific = TRUE)
p.adjusted = format(p.adjusted, digits = 3, scientific = TRUE)
res.df = cbind(enrichment[, 1], odds[, 1], test[, 1], p.adjusted[, 1])
for (l in 2:length(query.list)){
     res.df = cbind(res.df, cbind(enrichment[, l], odds[, l], test[, l], p.adjusted[, l]))
}
res.df = cbind(subject.names, res.df)
colnames(res.df) = c("subject",
                     as.vector(sapply(query.names,
                                      function(x){paste(x,
                                                        c("fold.enrichment",
                                                          "odds.ratio",
                                                          "p.value",
                                                          "p.value.adjusted"),
                                                        sep = "_")})))
res.enrichment = cbind(subject.names, enrichment)
colnames(res.enrichment) = c("subject", query.names)
res.odds = cbind(subject.names, odds)
colnames(res.odds) = c("subject", query.names)
res.test = cbind(subject.names, test)
colnames(res.test) = c("subject", query.names)
res.p.adjusted = cbind(subject.names, p.adjusted)
colnames(res.p.adjusted) = c("subject", query.names)
all.res = list(fold.enrichment = res.enrichment,
               odds.ratio = res.odds,
               p.value.original = res.test,
               p.value.adjusted = res.p.adjusted,
               all.results = res.df)
return(all.res)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/enrichtest.r")
#   query.list = "GVV/Query_files.txt"
#   subject.list = "GVV/Subject_files.txt"
#   res = enrichmentTest(query = query.list, subject = subject.list)
