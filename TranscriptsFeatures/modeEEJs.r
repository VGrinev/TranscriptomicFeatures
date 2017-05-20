##############################################################################################################
## High-level R functions for classification of experimentally detected exon-exon junctions                 ##
## according to modes (types) of alternative splicing.                                                      ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of functions:
##  models.hmRNA - path to folder and name of the GTF/GFF file with models of hypothetical "non-alternative"
#   precursor RNAs.
##  expEEJ.in - path to folder and name of the TXT file in tab-delimited format with experimentally detected
#   exon-exon junctions. This is input data file which should include four mandatory fields:
#   i) seqnames (name of chromosome or scaffold with prefix "chr");
#   ii) start (start genomic coordinate of the exon-exon junction);
#   iii) end (end genomic coordinate of the exon-exon junction);
#   iv) strand (strand information about exon-exon junction).
##  expEEJ.out - path to folder and name of the TXT output file in tab-delimited format.
modeStatistics = function(x){
CI = x[x$CI == 1, ]
CI = nrow(CI[rowSums(CI) == 1, ])
CE = x[x$CE == 1, ]
CE = nrow(CE[rowSums(CE) == 1, ])
MCE = x[x$MCE == 1, ]
MCE = nrow(MCE[rowSums(MCE) == 1, ])
MEE = x[x$MEE == 1, ]
MEE = nrow(MEE[rowSums(MEE[, -2:-3]) == 1, ])
A5SS = x[x$A5SS == 1, ]
A5SS = nrow(A5SS[rowSums(A5SS) == 1, ])
A3SS = x[x$A3SS == 1, ]
A3SS = nrow(A3SS[rowSums(A3SS) == 1, ])
II = x[x$II == 1, ]
II = nrow(II[rowSums(II[, -5:-6]) == 1, ])
AFE = x[x$AFE == 1, ]
AFE = nrow(AFE[rowSums(AFE[, -5]) == 1, ])
ALE = x[x$ALE == 1, ]
ALE = nrow(ALE[rowSums(ALE[, -6]) == 1, ])
CE.all = nrow(x[x$CE > 0 | x$MCE > 0 | x$MEE > 0, ])
A5SS.all = nrow(x[x$A5SS > 0 | x$II > 0 | x$AFE > 0, ])
A3SS.all = nrow(x[x$A3SS > 0 | x$II > 0 | x$ALE > 0, ])
sum.mode = cbind(c("All",
                   "CI",
                   "CE", "MCE", "MEE", "CE.all",
                   "A5SS", "A5SS.all", "A3SS", "A3SS.all", "II",
                   "AFE", "ALE",
                   "CoE",
                   "ND"),
                 c(nrow(x),
                   c(CI, CE, MCE, MEE, CE.all, A5SS, A5SS.all, A3SS, A3SS.all, II, AFE, ALE),
                   nrow(x) - sum(CI, CE, MCE, MEE, A5SS, A3SS, II, AFE, ALE, nrow(x[rowSums(x) == 0, ])),
                   nrow(x[rowSums(x) == 0, ])))
colnames(sum.mode) = c("mode", "count")
sum.mode = data.frame(sum.mode, stringsAsFactors = FALSE)
sum.mode$count = as.numeric(sum.mode$count)
sum.mode$frequency = sum.mode$count/sum.mode$count[1]
return(sum.mode)
}
modeEEJs = function(models.hmRNA, expEEJ.in, expEEJ.out){
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
## Loading of experimentally detected exon-exon junctions as a GRanges object.
expEEJ = read.table(file = expEEJ.in, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
expEEJ = makeGRangesFromDataFrame(expEEJ, keep.extra.columns = TRUE)
start(expEEJ) = start(expEEJ) + 1
end(expEEJ) = end(expEEJ) - 1
## Auxiliary objects.
#  A list of internal exons from preRNA models.
ex.internal = models.hmRNA[models.hmRNA$type == "exon", ]
ex.internal = ex.internal[!ex.internal %in%
                          subsetByOverlaps(query = ex.internal,
                                           subject = c(models.hmRNA[models.hmRNA$type == "first_exon", ],
                                                       models.hmRNA[models.hmRNA$type == "last_exon", ]),
                                           type = "equal"), ]
ex.internal = reduce(ex.internal)
#  A list of introns from preRNA models.
I = models.hmRNA[models.hmRNA$type == "intron", ]
I = I[!duplicated(I), ][, 0]
#  A list of retained introns from preRNA models.
IR = models.hmRNA[models.hmRNA$type == "retained_intron", ]
IR = IR[!duplicated(IR), ][, 0]
#  An GRanges object with genomic coordinates of 5' splice sites of experimentally detected EEJs.
fiveSS = expEEJ
end(fiveSS[as.character(strand(fiveSS)) == "+", ]) = start(fiveSS[as.character(strand(fiveSS)) == "+", ])
start(fiveSS[as.character(strand(fiveSS)) == "-", ]) = end(fiveSS[as.character(strand(fiveSS)) == "-", ])
#  An GRanges object with genomic coordinates of 3' splice sites of experimentally detected EEJs.
threeSS = expEEJ
start(threeSS[as.character(strand(threeSS)) == "+", ]) = end(threeSS[as.character(strand(threeSS)) == "+", ])
end(threeSS[as.character(strand(threeSS)) == "-", ]) = start(threeSS[as.character(strand(threeSS)) == "-", ])
#  An GRanges object of EEJs with cassette exons.
CE = expEEJ[subjectHits(findOverlaps(query = ex.internal, subject = expEEJ, type = "within")), ]
CE = CE[!duplicated(CE), ]
## Annotation of canonical EEJs (CI).
expEEJ$CI = countOverlaps(query = expEEJ, subject = I, type = "equal")
## Annotation of EEJs associated with intron retention (IR).
expEEJ$IR = countOverlaps(query = expEEJ, subject = IR, type = "equal")
## Annotation of EEJs with one cassette exon (CE).
hits = subjectHits(findOverlaps(query = ex.internal, subject = expEEJ, type = "within"))
expEEJ$CE = 0
expEEJ[as.numeric(attr(table(hits)[table(hits) == 1], "dimnames")[[1]]), ]$CE = 1
## Annotation of EEJs with multiple cassette exons (MCE).
expEEJ$MCE = 0
expEEJ[as.numeric(attr(table(hits)[table(hits) > 1], "dimnames")[[1]]), ]$MCE = 1
## Annotation of EEJs associated with mutually exclusive exons (MEE).
MEE = findOverlaps(query = CE, type = "any")
MEE = CE[unique(queryHits(MEE[!(MEE@from - MEE@to) == 0, ])), ]
MEE.regions = reduce(MEE)[countOverlaps(query = reduce(MEE),
                                        subject = ex.internal[countOverlaps(query = ex.internal,
                                                                            subject = MEE,
                                                                            type = "within") == 1, ],
                                        type = "any") == 2, ]
MEE.trimmed = MEE.regions
start(MEE.trimmed) = start(MEE.trimmed) + 10
end(MEE.trimmed) = end(MEE.trimmed) - 10
MEE.introns = expEEJ[expEEJ %in% subsetByOverlaps(query = expEEJ,
                                                  subject = MEE.trimmed,
                                                  type = "within"), ]
MEE.introns = MEE.introns[!MEE.introns %in% subsetByOverlaps(query = MEE.introns,
                                                             subject = MEE,
                                                             type = "equal"), ]
MEE.regions = MEE.regions[!MEE.regions %in% subsetByOverlaps(query = MEE.regions,
                                                             subject = MEE.introns,
                                                             type = "any"), ]
MEE = MEE[MEE %in% subsetByOverlaps(query = MEE,
                                     subject = MEE.regions,
                                     type = "within"), ]
expEEJ$MEE = countOverlaps(query = expEEJ,
                           subject = MEE,
                           type = "equal")
## Annotation of EEJs with alternative 5' splice sites (A5SS).
expEEJ$A5SS = 0
expEEJ[as.character(strand(expEEJ)) == "+", ]$A5SS = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "+", ],
                                                                   subject = I, type = "start")
expEEJ[as.character(strand(expEEJ)) == "-", ]$A5SS = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "-", ],
                                                                   subject = I, type = "end")
expEEJ$A5SS = c(0, 1, expEEJ$A5SS)[match(expEEJ$A5SS, c(1, 0, expEEJ$A5SS))]
expEEJ[countOverlaps(query = fiveSS, subject = models.hmRNA, type = "any") == 0, ]$A5SS = 0
## Annotation of EEJs with alternative 3' splice sites (A3SS).
expEEJ$A3SS = 0
expEEJ[as.character(strand(expEEJ)) == "+", ]$A3SS = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "+", ],
                                                                   subject = I, type = "end")
expEEJ[as.character(strand(expEEJ)) == "-", ]$A3SS = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "-", ],
                                                                   subject = I, type = "start")
expEEJ$A3SS = c(0, 1, expEEJ$A3SS)[match(expEEJ$A3SS, c(1, 0, expEEJ$A3SS))]
expEEJ[countOverlaps(query = threeSS, subject = models.hmRNA, type = "any") == 0, ]$A3SS = 0
## Annotation of EEJs with alternative both 5' and 3' splice sites (II).
expEEJ$II = 0
expEEJ[expEEJ$A5SS > 0 & expEEJ$A3SS > 0, ]$II = 1
## Annotation of EEJs as alternative first exons (AFE).
AFE = models.hmRNA[models.hmRNA$type == "alternative_first_exon", ][, 0]
end(AFE[as.character(strand(AFE)) == "+", ]) = end(AFE[as.character(strand(AFE)) == "+", ]) + 1
start(AFE[as.character(strand(AFE)) == "+", ]) = end(AFE[as.character(strand(AFE)) == "+", ])
start(AFE[as.character(strand(AFE)) == "-", ]) = start(AFE[as.character(strand(AFE)) == "-", ]) - 1
end(AFE[as.character(strand(AFE)) == "-", ]) = start(AFE[as.character(strand(AFE)) == "-", ])
expEEJ$AFE.full = 0
expEEJ[as.character(strand(expEEJ)) == "+", ]$AFE.full = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "+", ],
                                                                       subject = AFE,
                                                                       type = "start")
expEEJ[as.character(strand(expEEJ)) == "-", ]$AFE.full = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "-", ],
                                                                       subject = AFE,
                                                                       type = "end")
AFE = AFE[!AFE %in% subsetByOverlaps(query = AFE,
                                     subject = I,
                                     type = "start"), ]
AFE = AFE[!AFE %in% subsetByOverlaps(query = AFE,
                                     subject = I,
                                     type = "end"), ]
expEEJ$AFE = 0
expEEJ[as.character(strand(expEEJ)) == "+", ]$AFE = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "+", ],
                                                                  subject = AFE,
                                                                  type = "start")
expEEJ[as.character(strand(expEEJ)) == "-", ]$AFE = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "-", ],
                                                                  subject = AFE,
                                                                  type = "end")
## Annotation of EEJs as alternative last exons (ALE).
ALE = models.hmRNA[models.hmRNA$type == "alternative_last_exon", ][, 0]
start(ALE[as.character(strand(ALE)) == "+", ]) = start(ALE[as.character(strand(ALE)) == "+", ]) - 1
end(ALE[as.character(strand(ALE)) == "+", ]) = start(ALE[as.character(strand(ALE)) == "+", ])
end(ALE[as.character(strand(ALE)) == "-", ]) = end(ALE[as.character(strand(ALE)) == "-", ]) + 1
start(ALE[as.character(strand(ALE)) == "-", ]) = end(ALE[as.character(strand(ALE)) == "-", ])
expEEJ$ALE.full = 0
expEEJ[as.character(strand(expEEJ)) == "+", ]$ALE.full = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "+", ],
                                                                       subject = ALE,
                                                                       type = "end")
expEEJ[as.character(strand(expEEJ)) == "-", ]$ALE.full = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "-", ],
                                                                       subject = ALE,
                                                                       type = "start")
ALE = ALE[!ALE %in% subsetByOverlaps(query = ALE,
                                     subject = I,
                                     type = "start"), ]
ALE = ALE[!ALE %in% subsetByOverlaps(query = ALE,
                                     subject = I,
                                     type = "end"), ]
expEEJ$ALE = 0
expEEJ[as.character(strand(expEEJ)) == "+", ]$ALE = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "+", ],
                                                                  subject = ALE,
                                                                  type = "end")
expEEJ[as.character(strand(expEEJ)) == "-", ]$ALE = countOverlaps(query = expEEJ[as.character(strand(expEEJ)) == "-", ],
                                                                  subject = ALE,
                                                                  type = "start")
## Non-annotated EEJs (ND).
expEEJ$ND = 0
expEEJ[rowSums(cbind(expEEJ$CI,
                     expEEJ$CE,
                     expEEJ$MCE,
                     expEEJ$MEE,
                     expEEJ$A5SS,
                     expEEJ$A3SS,
                     expEEJ$II,
                     expEEJ$AFE,
                     expEEJ$ALE)) == 0, ]$ND = 1
## Calculation of overall statistics.
full = elementMetadata(expEEJ)
full = data.frame(full[, (length(full) - 12):length(full)])[, c(1, 3:8, 10, 12)]
full[full > 0] = 1
dif.a = expEEJ[abs(expEEJ$logFC) >= 1 & expEEJ$q_value < 0.05, ]
dif.m = elementMetadata(dif.a[dif.a$logFC <= -1, ])
dif.m = data.frame(dif.m[, (length(dif.m) - 12):length(dif.m)])[, c(1, 3:8, 10, 12)]
dif.m[dif.m > 0] = 1
dif.p = elementMetadata(dif.a[dif.a$logFC >= 1, ])
dif.p = data.frame(dif.p[, (length(dif.p) - 12):length(dif.p)])[, c(1, 3:8, 10, 12)]
dif.p[dif.p > 0] = 1
dif.n = elementMetadata(expEEJ[!expEEJ$Event_ID %in% dif.a$Event_ID, ])
dif.n = data.frame(dif.n[, (length(dif.n) - 12):length(dif.n)])[, c(1, 3:8, 10, 12)]
dif.n[dif.n > 0] = 1
dif.a = elementMetadata(dif.a)
dif.a = data.frame(dif.a[, (length(dif.a) - 12):length(dif.a)])[, c(1, 3:8, 10, 12)]
dif.a[dif.a > 0] = 1
sub.sets = list(full, dif.a, dif.m, dif.p, dif.n)
sum.stat = lapply(sub.sets, function(x){modeStatistics(x)})
sum.stat = data.frame(Reduce(cbind, sum.stat))
sum.stat = sum.stat[, c(1:3, 5:6, 8:9, 11:12, 14:15)]
colnames(sum.stat) = c("mode",
                       "all_count", "all_frequency",
                       "all.dif_count", "all.dif_frequency",
                       "minus.dif_count", "minus.dif_frequency",
                       "plus.dif_count", "plus.dif_frequency",
                       "non.dif_count", "non.dif_frequency")
write.table(sum.stat,
            file = sub(".txt", ", summary statistics.txt", out.data),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
## Saving final results.
write.table(expEEJ,
            file = expEEJ.out,
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
return(expEEJ)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/modeeejs.r")
#   models = "GVV/Ensembl_hnapRNA.gtf"
#   in.data = "GVV/Experimental_junctions.txt"
#   out.data = "GVV/Experimental_junctions, Ensembl-based modes.txt"
#   res = modeEEJs(models.hmRNA = models, expEEJ.in = in.data, expEEJ.out = out.data)
