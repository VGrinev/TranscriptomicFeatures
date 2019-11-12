###################################################################################################
## High-level R function for classification of exon-exon junctions                               ##
## according to modes (or types) of alternative splicing.                                        ##
## (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                       ##
###################################################################################################
### Arguments of functios:
##  models.hmRNA - path to folder and name of the GTF/GFF file with models of
#   hypothetical "non-alternative" precursor RNAs;
##  eej.in - path to folder and name of the TXT file in tab-delimited format with
#   exon-exon junctions. This is input data file which should include four mandatory fields:
#   i) seqnames (name of chromosome or scaffold with prefix "chr");
#   ii) start (start genomic coordinate of the exon-exon junction);
#   iii) end (end genomic coordinate of the exon-exon junction);
#   iv) strand (strand information about exon-exon junction).
##  eej.out - path to folder and name of the TXT output file in tab-delimited format.
modeEEJs.classifier = function(models.hmRNA, eej.in, eej.out){
## Loading of required auxiliary libraries.
#  This code was successfully tested with library GenomicFeatures v.1.26.3 (and higher)
#  and library rtracklayer v.1.34.2 (and higher).
if (!suppressMessages(require(rtracklayer))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
}
## Loading of the preRNA models as a GRanges object.
models.hmRNA = import(con = models.hmRNA)
## Loading of exon-exon junctions to be classified as a GRanges object.
eej = read.table(file = eej.in, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
eej = makeGRangesFromDataFrame(eej, keep.extra.columns = TRUE)
start(eej) = start(eej) + 1
end(eej) = end(eej) - 1
eej = eej[width(eej) >= 50, ]
## Auxiliary objects.
#  A list of internal exons from preRNA models.
ex.inter = models.hmRNA[models.hmRNA$type == "exon", ]
ex.inter = ex.inter[!ex.inter %in%
                    subsetByOverlaps(x = ex.inter,
                                     ranges = c(models.hmRNA[models.hmRNA$type == "first_exon", ],
                                                models.hmRNA[models.hmRNA$type == "last_exon", ]),
                                     type = "equal"), 0]
#  A list of introns from preRNA models.
intr = models.hmRNA[models.hmRNA$type == "intron", ]
#  An GRanges object with genomic coordinates of 5' splice sites of exon-exon junctions.
f.ss = eej
end(f.ss[as.character(strand(f.ss)) == "+", ]) = start(f.ss[as.character(strand(f.ss)) == "+", ])
start(f.ss[as.character(strand(f.ss)) == "-", ]) = end(f.ss[as.character(strand(f.ss)) == "-", ])
#  An GRanges object with genomic coordinates of 3' splice sites of exon-exon junctions.
t.ss = eej
start(t.ss[as.character(strand(t.ss)) == "+", ]) = end(t.ss[as.character(strand(t.ss)) == "+", ])
end(t.ss[as.character(strand(t.ss)) == "-", ]) = start(t.ss[as.character(strand(t.ss)) == "-", ])
#  An GRanges object of exon-exon junctions with cassette exons.
CE = eej[subjectHits(findOverlaps(query = ex.inter, subject = eej, type = "within")), ]
CE = CE[!duplicated(CE), ]
## Annotation of canonical exon-exon junctions (CI).
eej$CI = countOverlaps(query = eej, subject = intr, type = "equal")
## Annotation of exon-exon junctions associated with retention of introns (RI).
RI = c(models.hmRNA[models.hmRNA$type == "first_exon", ],
       models.hmRNA[models.hmRNA$type == "alternative_first_exon", ],
       models.hmRNA[models.hmRNA$type == "exon", ],
       models.hmRNA[models.hmRNA$type == "alternative_last_exon", ],
       models.hmRNA[models.hmRNA$type == "last_exon", ])[, 0]
RI = RI[!duplicated(RI), ]
eej$RI = 0
eej$RI = countOverlaps(query = eej, subject = RI, type = "within")
## Annotation of exon-exon junctions with one cassette exon (CE).
hits = subjectHits(findOverlaps(query = ex.inter, subject = eej, type = "within"))
eej$CE = 0
eej[as.numeric(attr(table(hits)[table(hits) == 1], "dimnames")[[1]]), ]$CE = 1
## Annotation of exon-exon junctions with multiple cassette exons (MCE).
eej$MCE = 0
eej[as.numeric(attr(table(hits)[table(hits) > 1], "dimnames")[[1]]), ]$MCE = 1
## Annotation of exon-exon junctions associated with mutually exclusive exons (MEE).
MEE = findOverlaps(query = CE, type = "any")
MEE = CE[unique(queryHits(MEE[!(MEE@from - MEE@to) == 0, ])), ]
MEE.regions = reduce(MEE)[countOverlaps(query = reduce(MEE),
                                        subject = ex.inter[countOverlaps(query = ex.inter,
                                                                         subject = MEE,
                                                                         type = "within") == 1, ],
                                        type = "any") == 2, ]
MEE.trimmed = MEE.regions
start(MEE.trimmed) = start(MEE.trimmed) + 10
end(MEE.trimmed) = end(MEE.trimmed) - 10
MEE.introns = eej[eej %in% subsetByOverlaps(x = eej,
                                            ranges = MEE.trimmed,
                                            type = "within"), ]
MEE.introns = MEE.introns[!MEE.introns %in% subsetByOverlaps(x = MEE.introns,
                                                             ranges = MEE,
                                                             type = "equal"), ]
MEE.regions = MEE.regions[!MEE.regions %in% subsetByOverlaps(x = MEE.regions,
                                                             ranges = MEE.introns,
                                                             type = "any"), ]
MEE = MEE[MEE %in% subsetByOverlaps(x = MEE,
                                    ranges = MEE.regions,
                                    type = "within"), ]
eej$MEE = countOverlaps(query = eej,
                           subject = MEE,
                           type = "equal")
## Annotation of exon-exon junctions with alternative 5' splice sites (A5SS).
eej$A5SS = 0
eej[as.character(strand(eej)) == "+", ]$A5SS = countOverlaps(query = eej[as.character(strand(eej)) == "+", ],
                                                             subject = intr, type = "start")
eej[as.character(strand(eej)) == "-", ]$A5SS = countOverlaps(query = eej[as.character(strand(eej)) == "-", ],
                                                             subject = intr, type = "end")
eej$A5SS = c(0, 1, eej$A5SS)[match(eej$A5SS, c(1, 0, eej$A5SS))]
if (length(table(countOverlaps(query = f.ss, subject = models.hmRNA, type = "any") == 0)) > 1){
eej[countOverlaps(query = f.ss, subject = models.hmRNA, type = "any") == 0, ]$A5SS = 0
}
## Annotation of exon-exon junctions with alternative 3' splice sites (A3SS).
eej$A3SS = 0
eej[as.character(strand(eej)) == "+", ]$A3SS = countOverlaps(query = eej[as.character(strand(eej)) == "+", ],
                                                             subject = intr, type = "end")
eej[as.character(strand(eej)) == "-", ]$A3SS = countOverlaps(query = eej[as.character(strand(eej)) == "-", ],
                                                             subject = intr, type = "start")
eej$A3SS = c(0, 1, eej$A3SS)[match(eej$A3SS, c(1, 0, eej$A3SS))]
if (length(table(countOverlaps(query = t.ss, subject = models.hmRNA, type = "any") == 0)) > 1){
eej[countOverlaps(query = t.ss, subject = models.hmRNA, type = "any") == 0, ]$A3SS = 0
}
## Annotation of exon-exon junctions with alternative both 5' and 3' splice sites (II).
eej$II = 0
eej[eej$A5SS > 0 & eej$A3SS > 0, ]$II = 1
eej[eej$RI > 0, ]$A5SS = 0
eej[eej$RI > 0, ]$A3SS = 0
eej[eej$RI > 0, ]$II = 0
## Annotation of exon-exon junctions as alternative first exons (AFE).
AFE = models.hmRNA[models.hmRNA$type == "alternative_first_exon", ][, 0]
end(AFE[as.character(strand(AFE)) == "+", ]) = end(AFE[as.character(strand(AFE)) == "+", ]) + 1
start(AFE[as.character(strand(AFE)) == "-", ]) = start(AFE[as.character(strand(AFE)) == "-", ]) - 1
eej$AFE = 0
eej$AFE = countOverlaps(query = f.ss, subject = AFE, type = "any")
## Annotation of exon-exon junctions as alternative last exons (ALE).
ALE = models.hmRNA[models.hmRNA$type == "alternative_last_exon", ][, 0]
start(ALE[as.character(strand(ALE)) == "+", ]) = start(ALE[as.character(strand(ALE)) == "+", ]) - 1
end(ALE[as.character(strand(ALE)) == "-", ]) = end(ALE[as.character(strand(ALE)) == "-", ]) + 1
eej$ALE = 0
eej$ALE = countOverlaps(query = t.ss, subject = ALE, type = "any")
## Non-annotated EEJs (ND).
eej$ND = 1
eej[rowSums(cbind(eej$CI,
                  eej$RI,
                  eej$CE,
                  eej$MCE,
                  eej$MEE,
                  eej$A5SS,
                  eej$A3SS,
                  eej$II,
                  eej$AFE,
                  eej$ALE)) > 0, ]$ND = 0
## Saving final results.
write.table(eej,
            file = eej.out,
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
return(eej)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/modeeejs.classifier.r")
models = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, expression-based hnapRNA of the RUNX1-RUNX1T1 oncogene.gtf"
in.data = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, all EEJs.txt"
out.data = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, all EEJs, expression-based modes.txt"
res = modeEEJs.classifier(models.hmRNA = models, eej.in = in.data, eej.out = out.data)
