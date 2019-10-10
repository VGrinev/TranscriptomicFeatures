###################################################################################################
##  High-level R function for ranking of TSSs according to position in a gene.                   ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work - path to and name of work directory.
##  in.file - the name of file with TSSs data. It is usually output file of diffTSSsUsage function.
##  f.rankTSSs - the name of output TXT file in tab-delimited format for storing the final results.
#   Default value is NULL.
rankTSSs = function(d.work, in.file, f.rankTSSs = NULL){
##  Loading, filtering and sub-setting of input data.
in.file = read.table(file = paste(d.work, in.file, sep = "/"),
                     sep = "\t",
                     header = TRUE,
                     quote = "\"",
                     as.is = TRUE)
d.str = cbind(in.file$strand, in.file$gene_name)
d.str = table(d.str[!duplicated(x = d.str), ][, 2])
in.file = in.file[!in.file$gene_name %in% attr(x = d.str[d.str > 1], which = "dimnames")[[1]], ]
in.file = in.file[!gsub(pattern = ".+,.+", replacement = ",", x = in.file$gene_name) == ",", ]
in.file = split(x = in.file, f = in.file$gene_name)
##  Calculation of the TSSs ranks.
ranks = do.call(what = rbind,
                args = lapply(X = in.file,
                              FUN = function(y){if (unique(x = y$strand) == "+"){
                                                    y = y[order(y$start), ]
                                                    y$tss_rank = c(1:nrow(y))
                                               }else{
                                                    y = y[order(-y$end), ]
                                                    y$tss_rank = c(1:nrow(y))
                                               }
                                               return(y)}))
ranks = ranks[order(ranks$tss_id), ]
rownames(ranks) = NULL
ranks = ranks[, c(1, ncol(ranks), 2:ncol(ranks))]
##  Saving of the final results.
if (!is.null(f.rankTSSs)){
write.table(ranks,
            file = paste(d.work, f.rankTSSs, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
}
return(ranks)
}
### A simple example of function use.
res = rankTSSs(d.work = "D:/Vasily Grinev",
               in.file = "THP-1 cells, KAT8 transcriptome, TSSs matrix.txt",
               f.rankTSSs = "THP-1 cells, KAT8 transcriptome, TSSs matrix, with ranks.txt")
