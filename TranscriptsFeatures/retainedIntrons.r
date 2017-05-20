##############################################################################################################
## High-level R function for reference-based annotation of exon bins as retained introns.                   ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
##  models.hmRNA - path to folder and name of the GTF/GFF file with models of hypothetical "non-alternative"
#   precursor RNAs. It is output file of hnapRNAgenerator function.
##  exBins - path to folder and name of the TXT file in tab-delimited format with experimentally detected exon
#   bins. It is input file which is based on output data of JunctionSeq (or DEXseq) package. This file should
#   include six mandatory fields:
#   i) seqnames (name of chromosome or scaffold with prefix "chr");
#   ii) start (start genomic coordinate of the exon bin);
#   iii) end (end genomic coordinate of the exon bin);
#   iv) strand (strand information about exon bin);
#   v) logFC (log2 fold change in usage of the exon bin between two conditions);
#   vi) q_value (FDR-adjusted p-value for logFC).
##  log.FC - an integer argument. It is a threshold for log2 fold change. Default value is 1.
##  q.value - an integer argument. It is a threshold for FDR-adjusted p-value. Default value is 0.05.
##  overlap.type - type of overlap between exon bin and retained intron. Default value is "within".
##  stat - path to folder and name of the TXT output file in tab-delimited format with summary statistics.
retainedIntrons = function(models.hmRNA, exBins, log.FC = 1, q.value = 0.05, overlap.type = "within", stat){
## Loading of required auxiliary libraries.
#  This code was successfully tested with library GenomicFeatures v.1.26.3 (and higher)
#  and library rtracklayer v.1.34.2 (and higher).
if (!suppressMessages(require(GenomicFeatures))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("GenomicFeatures")
}
if (!suppressMessages(require(rtracklayer))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
}
## Loading of the preRNA models as a GRanges object.
models.hmRNA = import(con = models.hmRNA)
## Retrieving a list of retained introns from preRNA models.
RI = models.hmRNA[models.hmRNA$type == "retained_intron", ]
RI = RI[!duplicated(RI), ][, 0]
## Loading of experimentally detected exon bins as a GRanges object.
ex.bins = read.table(file = exBins, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
ex.bins = makeGRangesFromDataFrame(ex.bins, keep.extra.columns = TRUE)
## Sub-setting of experimentally detected exon bins.
dif.a = ex.bins[abs(ex.bins$logFC) >= log.FC & ex.bins$q_value < q.value, ]
dif.m = dif.a[dif.a$logFC <= -log.FC, ][, 0]
dif.p = dif.a[dif.a$logFC >= log.FC, ][, 0]
dif.n = ex.bins[!ex.bins$event_id %in% dif.a$event_id, ]
sub.sets = list(dif.n = dif.n, dif.a = dif.a, dif.m = dif.m, dif.p = dif.p)
## Calculation of overall statistics.
sum.stat = lapply(sub.sets, function(x){c(length(x),
                                          length(x) - length(x[x %in% subsetByOverlaps(query = x,
                                                                                       subject = RI,
                                                                                       type = overlap.type), ]),
                                          100 - round(length(x[x %in% subsetByOverlaps(query = x,
                                                                                       subject = RI,
                                                                                       type = overlap.type), ])/
                                                      length(x), digits = 3) * 100,
                                          length(x[x %in% subsetByOverlaps(query = x,
                                                                           subject = RI,
                                                                           type = overlap.type), ]),
                                          round(length(x[x %in% subsetByOverlaps(query = x,
                                                                                 subject = RI,
                                                                                 type = overlap.type), ])/
                                                      length(x), digits = 3) * 100)})
sum.stat = cbind(names(sub.sets), do.call(rbind, sum.stat))
rownames(sum.stat) = NULL
colnames(sum.stat) = c("sub.set",
                       "all_counts",
                       "non_ri.counts", "non_ri.frequency",
                       "ri.counts", "ri.frequency")
write.table(sum.stat,
            file = stat,
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
return(sum.stat)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/retainedintrons.r")
#   models = "GVV/Ensembl_hnapRNA.gtf"
#   exons = "GVV/Experimental_exon_bins.txt"
#   out.data = "GVV/RI summary statistics.txt"
#   retainedIntrons(models.hmRNA = models, exBins = exons, stat = out.data)
