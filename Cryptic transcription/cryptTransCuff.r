###################################################################################################
##  High-level R function for identification of cryptic intragenic transcription initiation      ##
##  based on Cufflinks assembled and Cuffdiff quantified transcripts.                            ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work          - character string giving the path to and name of work directory.
##  f.isoforms      - character string giving the name of *.fpkm_tracking file with transcript
#                     FPKMs. It is Cuffdiff output file.
##  f.transcriptome - character string giving the name of file in GTF/GFF format with assembled
#                     transcriptome. It is typically output file of Cufflinks/Cuffmerge.
##  multiTrans      - numerical threshold value for filtration of multi-transcript genes. Default
#                     value is 1 (all genes with at least 1 transcript will be saved
#                     for subsequent processing).
##  ctrl            - character string giving the name of control condition.
##  exper           - character string giving the name of experimental condition.
cryptTransCuff = function(d.work, f.isoforms, f.transcriptome, multiTrans = 1, ctrl, exper){
##  Loading of required auxiliary library.
#   This code was successfully tested with library rtracklayer v.1.44.3.
suppressMessages(require(rtracklayer))
##  Loading and filtering of RNA isoforms FPKM tracking data.
isoforms = read.table(file = paste(d.work, f.isoforms, sep = "/"),
                      sep = "\t",
                      header = TRUE,
                      quote = "\"",
                      as.is = TRUE)[, c(1, 5, 10, 14)]
colnames(isoforms) = c("transcript_id", "gene_name", colnames(isoforms)[3:4])
isoforms = isoforms[rowSums(isoforms[, 3:4] >= 1) > 0, ]
isoforms = isoforms[!isoforms$gene_name == "-", ]
isoforms = isoforms[order(isoforms$transcript_id), ]
##  Loading and sub-setting of transcriptome data.
trans = import(con = paste(d.work, f.transcriptome, sep = "/"), format = "gtf")
trans = trans[, c("gene_name", "transcript_id")]
trans = trans[trans$gene_name %in% isoforms$gene_name, ]
trans = trans[trans$transcript_id %in% isoforms$transcript_id, ]
trans = trans[trans$transcript_id %in% attr((table(trans$transcript_id) >= 5)
             [(table(trans$transcript_id) >= 5) == "TRUE"], "dimnames")[[1]], ]
if (multiTrans > 0){
    multi.trans = trans[, 1:2]
    multi.trans = multi.trans[!duplicated(multi.trans), ]
    multi.trans = attr((table(multi.trans$gene_name) >= multiTrans)
                       [(table(multi.trans$gene_name) >= multiTrans) == "TRUE"], "dimnames")[[1]]
    trans = trans[trans$gene_name %in% multi.trans, ]
}
##  Annotation of transcriptome data with FPKM tracking data.
trans$exon_id = paste(paste(paste(seqnames(trans),
                                  start(trans),
                                  sep = ":"),
                            end(trans),
                            sep = "-"),
                      strand(trans),
                      sep = "_str")
trans = as.data.frame(x = trans)
trans$seqnames = as.character(trans$seqnames)
trans$strand = as.character(trans$strand)
trans = trans[, c(6:8, 1:3, 5)]
trans$ctrl = 0
trans$exper = 0
colnames(trans) = c(colnames(trans)[1:7], ctrl, exper)
rownames(isoforms) = isoforms[, 1]
isoforms = as.matrix(x = isoforms)
trans = as.matrix(x = trans)
idx = trans[, 2] %in% isoforms[, 1]
trans[idx, 8] = isoforms[trans[idx, 2], 3]
trans[idx, 9] = isoforms[trans[idx, 2], 4]
trans = data.frame(trans, stringsAsFactors = FALSE)
trans[, c(5, 6, 8, 9)] = apply(X = trans[, c(5, 6, 8, 9)], MARGIN = 2, FUN = as.numeric)
exonFPKM = aggregate(x = trans[, 8:9], by = list(exon_id = trans$exon_id), FUN = sum)
rownames(exonFPKM) = exonFPKM[, 1]
exonFPKM = as.matrix(x = exonFPKM)
trans = trans[, c(1, 3:7)]
trans = as.matrix(x = trans[!duplicated(trans), ])
trans = cbind(trans, 0, 0)
idx = trans[, 2] %in% exonFPKM[, 1]
trans[idx, 7] = exonFPKM[trans[idx, 2], 2]
trans[idx, 8] = exonFPKM[trans[idx, 2], 3]
trans = data.frame(trans, stringsAsFactors = FALSE)
colnames(trans) = c(colnames(trans)[1:6], ctrl, exper)
trans[, c(4, 5, 7, 8)] = apply(X = trans[, c(4, 5, 7, 8)], MARGIN = 2, FUN = as.numeric)
rownames(trans) = NULL
##  Exon numbering.
trans = split(x = trans, f = trans$gene_name)
trans = lapply(X = trans,
               FUN = function(x){if (x$strand[1] == "+"){
                                     x = x[order(x$start, x$end), ]
                                     x$ex.n = 1:nrow(x)
                                }else{
                                     x = x[order(-x$end, -x$start), ]
                                     x$ex.n = 1:nrow(x)
                                }
                                x = x[, c(7:9)]
                                return(x)})
##  Calculation of a ratio between the first exon and all other exons of a gene.
#   The ratio "second-to-first exon".
second = lapply(X = trans,
                FUN = function(x){ctrl = x[x$ex.n == "2", 1]/x[x$ex.n == "1", 1]
                                  exper = x[x$ex.n == "2", 2]/x[x$ex.n == "1", 2]
                                  return(c(ctrl, exper))})
second = data.frame(do.call(what = rbind, args = second), stringsAsFactors = FALSE)
second$gene_name = rownames(second)
rownames(second) = NULL
second = second[, c(3, 1:2)]
colnames(second) = c("gene_name", "ctrl", "exper")
second = second[!is.na(second[, 2]) & !is.infinite(second[, 2]) & !is.nan(second[, 2]), ]
second = second[!is.na(second[, 3]) & !is.infinite(second[, 3]) & !is.nan(second[, 3]), ]
#   The ratio "third-to-first exon".
third = lapply(X = trans,
               FUN = function(x){ctrl = x[x$ex.n == "3", 1]/x[x$ex.n == "1", 1]
                                 exper = x[x$ex.n == "3", 2]/x[x$ex.n == "1", 2]
                                 return(c(ctrl, exper))})
third = data.frame(do.call(what = rbind, args = third), stringsAsFactors = FALSE)
third$gene_name = rownames(third)
rownames(third) = NULL
third = third[, c(3, 1:2)]
colnames(third) = c("gene_name", "ctrl", "exper")
third = third[!is.na(third[, 2]) & !is.infinite(third[, 2]) & !is.nan(third[, 2]), ]
third = third[!is.na(third[, 3]) & !is.infinite(third[, 3]) & !is.nan(third[, 3]), ]
#   The ratio "intermediate-to-first exon".
inter = lapply(X = trans,
               FUN = function(x){ctrl = mean(x[4:(nrow(x) - 1), 1])/x[x$ex.n == "1", 1]
                                 exper = mean(x[4:(nrow(x) - 1), 2])/x[x$ex.n == "1", 2]
                                 return(c(ctrl, exper))})
inter = data.frame(do.call(what = rbind, args = inter), stringsAsFactors = FALSE)
inter$gene_name = rownames(inter)
rownames(inter) = NULL
inter = inter[, c(3, 1:2)]
colnames(inter) = c("gene_name", "ctrl", "exper")
inter = inter[!is.na(inter[, 2]) & !is.infinite(inter[, 2]) & !is.nan(inter[, 2]), ]
inter = inter[!is.na(inter[, 3]) & !is.infinite(inter[, 3]) & !is.nan(inter[, 3]), ]
#   The ratio "last-to-first exon".
last = lapply(X = trans,
              FUN = function(x){ctrl = x[nrow(x), 1]/x[x$ex.n == "1", 1]
                                exper = x[nrow(x), 2]/x[x$ex.n == "1", 2]
                                return(c(ctrl, exper))})
last = data.frame(do.call(what = rbind, args = last), stringsAsFactors = FALSE)
last$gene_name = rownames(last)
rownames(last) = NULL
last = last[, c(3, 1:2)]
colnames(last) = c("gene_name", "ctrl", "exper")
last = last[!is.na(last[, 2]) & !is.infinite(last[, 2]) & !is.nan(last[, 2]), ]
last = last[!is.na(last[, 3]) & !is.infinite(last[, 3]) & !is.nan(last[, 3]), ]
##  Returning final results.
return(list(second = second, third = third, inter = inter, last = last))
}
