##############################################################################################################
## High-level R function for annotation of experimentally detected exon-exon junctions                      ##
## with distances to the nearest epigenetic marks.                                                          ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
##  eej - path to folder and name of the txt file in tab-delimited format with experimentally detected
#   exon-exon junctions. This is input data file which should include eight mandatory fields:
#   i) event_id (ID of the splicing event);
#   ii) ensembl_gene_id (Ensembl ID of the specified genes);
#   iii) seqnames (name of chromosome or scaffold with prefix "chr");
#   iv) start (start genomic coordinate of the exon-exon junction);
#   v) end (end genomic coordinate of the exon-exon junction);
#   vi) strand (strand information about exon-exon junction);
#   vii) logFC (log2 fold change in usage of the splicing event between two conditions);
#   viii) q_value (FDR-adjusted p-value for logFC).
##  features - path to folder and name of the txt file in tab-delimited format with a list of epigenetic marks
#   to be analysed. This file should include two mandatory fields:
#   i) path (path to folder with file(-s) of epigenetic marks);
#   ii) file (name of file with genomic coordinates of peaks of specified epigenetic marks).
#   A file with epigenetic marks should be standard txt file in tab-delimited format with following fields:
#   i) seqnames (name of chromosome or scaffold with prefix "chr");
#   ii) start (start genomic coordinate of the epigenetic peak);
#   iii) end (end genomic coordinate of the epigenetic peak).
#   It is strongly recommended to name the file by the name of the epigenetic mark. 
##  log.FC - an integer argument. It is a threshold for log2 fold change. Default value is 1.
##  q.value - an integer argument. It is a threshold for FDR-adjusted p-value. Default value is 0.05.
annotateEpigeneticFeatures = function(eej, features, log.FC = 1, q.value = 0.05){
## Loading of required auxiliary library.
if (!suppressMessages(require(GenomicFeatures))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("GenomicFeatures")
}
## Loading of experimentally detected EEJs as a GRanges object.
expEEJ = read.table(file = eej, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
expEEJ = makeGRangesFromDataFrame(expEEJ, keep.extra.columns = TRUE)
start(expEEJ) = start(expEEJ) + 1
end(expEEJ) = end(expEEJ) - 1
#  An GRanges object with genomic coordinates of 5' splice sites of experimentally detected EEJs.
f.ss = expEEJ[, 0]
end(f.ss[as.character(strand(f.ss)) == "+", ]) = start(f.ss[as.character(strand(f.ss)) == "+", ])
start(f.ss[as.character(strand(f.ss)) == "-", ]) = end(f.ss[as.character(strand(f.ss)) == "-", ])
#  An GRanges object with genomic coordinates of 3' splice sites of experimentally detected EEJs.
t.ss = expEEJ[, 0]
start(t.ss[as.character(strand(t.ss)) == "+", ]) = end(t.ss[as.character(strand(t.ss)) == "+", ])
end(t.ss[as.character(strand(t.ss)) == "-", ]) = start(t.ss[as.character(strand(t.ss)) == "-", ])
## Loading of epigenetic features as a GRanges object.
listFeatures = read.table(file = features, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
listFeatures = apply(listFeatures,
                     1,
                     function(x){makeGRangesFromDataFrame(read.table(file = paste(x[1], x[2], sep = "/"),
                                                                     sep = "\t",
                                                                     header = TRUE,
                                                                     quote = "\"",
                                                                     as.is = TRUE),
                                                           keep.extra.columns = TRUE)})
names(listFeatures) = gsub(".txt", "", read.table(file = features, sep = "\t", header = TRUE)$file)
## Annotation of splice sites with distances to epigenetic features.
for(i in 1:length(listFeatures)){
    u = follow(x = f.ss, subject = listFeatures[[i]], select = "all", ignore.strand = TRUE)
    d = precede(x = f.ss, subject = listFeatures[[i]], select = "all", ignore.strand = TRUE)
    o = findOverlaps(query = f.ss, subject = listFeatures[[i]], type = "any", ignore.strand = TRUE)
    expEEJ$u = "NA"
    expEEJ[queryHits(u), ]$u = start(f.ss[queryHits(u), ]) - end(listFeatures[[i]][subjectHits(u), ]) - 1
    expEEJ[queryHits(o), ]$u = 0
    expEEJ$d = "NA"
    expEEJ[queryHits(d), ]$d = start(listFeatures[[i]][subjectHits(d), ]) - end(f.ss[queryHits(d), ]) - 1
    expEEJ[queryHits(o), ]$d = 0
    strand = as.character(strand(expEEJ))
    expEEJ[strand == "-", c("u", "d")] = expEEJ[strand == "-", c("d", "u")]
    names(mcols(expEEJ))[ncol(mcols(expEEJ)) - 1] = paste("fiveSS.up", names(listFeatures[i]), sep = ".")
    names(mcols(expEEJ))[ncol(mcols(expEEJ))] = paste("fiveSS.down", names(listFeatures[i]), sep = ".")
    u = follow(x = t.ss, subject = listFeatures[[i]], select = "all", ignore.strand = TRUE)
    d = precede(x = t.ss, subject = listFeatures[[i]], select = "all", ignore.strand = TRUE)
    o = findOverlaps(query = t.ss, subject = listFeatures[[i]], type = "any", ignore.strand = TRUE)
    expEEJ$u = "NA"
    expEEJ[queryHits(u), ]$u = start(t.ss[queryHits(u), ]) - end(listFeatures[[i]][subjectHits(u), ]) - 1
    expEEJ[queryHits(o), ]$u = 0
    expEEJ$d = "NA"
    expEEJ[queryHits(d), ]$d = start(listFeatures[[i]][subjectHits(d), ]) - end(t.ss[queryHits(d), ]) - 1
    expEEJ[queryHits(o), ]$d = 0
    expEEJ[strand == "-", c("u", "d")] = expEEJ[strand == "-", c("d", "u")]
    names(mcols(expEEJ))[ncol(mcols(expEEJ)) - 1] = paste("threeSS.up", names(listFeatures[i]), sep = ".")
    names(mcols(expEEJ))[ncol(mcols(expEEJ))] = paste("threeSS.down", names(listFeatures[i]), sep = ".")
}
## Consolidation and saving all results as a standard data frame object.
res.all = as.data.frame(expEEJ)
fiveSS.up = paste("fiveSS.up", names(listFeatures), sep = ".")
fiveSS.down = paste("fiveSS.down", names(listFeatures), sep = ".")
threeSS.up = paste("threeSS.up", names(listFeatures), sep = ".")
threeSS.down = paste("threeSS.down", names(listFeatures), sep = ".")
columns = c("event_id", "ensembl_gene_id", "seqnames", "start", "end", "strand",
            "logFC", "q_value",
            fiveSS.up, fiveSS.down, threeSS.up, threeSS.down)
res.all = res.all[, c(columns)]
write.table(res.all,
            file = sub(".txt", ", distances to epimarks.txt", eej),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
## Consolidation and saving PDFs as a standard data frame object.
#  logFC- and q value-based sub-setting of data.
dif.a = res.all[abs(res.all$logFC) >= log.FC & res.all$q_value < q.value, ]
dif.n = res.all[!res.all$event_id %in% dif.a$event_id, ]
dif.m = dif.a[dif.a$logFC <= -1, ]
dif.p = dif.a[dif.a$logFC >= 1, ]
sub.sets = list(dif.non = dif.n, dif.minus = dif.m, dif.plus = dif.p)
#  Calculation and saving individual PDFs for each sub-set of data.
for (j in 1:length(listFeatures)){
     dist.fiveSS = list()
     dist.threeSS = list()
     for (k in 1:3){
          up = suppressWarnings(na.omit(as.numeric(sub.sets[[k]][, fiveSS.up[j]])))
          up = up[up > 0]
          down = suppressWarnings(na.omit(as.numeric(sub.sets[[k]][, fiveSS.down[j]])))
          zero = down[down == 0]
          down = down[down > 0]
          dist.fiveSS[[k]] = cbind(density(c(-log10(up + 1), log10(zero + 1), log10(down + 1)))$x,
                                   density(c(-log10(up + 1), log10(zero + 1), log10(down + 1)))$y)
          up = suppressWarnings(na.omit(as.numeric(sub.sets[[k]][, threeSS.up[j]])))
          up = up[up > 0]
          down = suppressWarnings(na.omit(as.numeric(sub.sets[[k]][, threeSS.down[j]])))
          zero = down[down == 0]
          down = down[down > 0]
          dist.threeSS[[k]] = cbind(density(c(-log10(up + 1), log10(zero + 1), log10(down + 1)))$x,
                                    density(c(-log10(up + 1), log10(zero + 1), log10(down + 1)))$y)
     }
     dist.fiveSS = do.call(cbind, dist.fiveSS)
     colnames(dist.fiveSS) = c("density", "fiveSS.distances_dif.non",
                               "density", "fiveSS.distances_dif.minus",
                               "density", "fiveSS.distances_dif.plus")
     dist.threeSS = do.call(cbind, dist.threeSS)
     colnames(dist.threeSS) = c("density", "threeSS.distances_dif.non",
                                "density", "threeSS.distances_dif.minus",
                                "density", "threeSS.distances_dif.plus")
     file.name = paste(paste(sub(".txt", ", distances to", eej), names(listFeatures)[j]), ".txt", sep = "")
     write.table(cbind(dist.fiveSS, dist.threeSS),
                 file = file.name,
                 sep = "\t",
                 quote = FALSE,
                 col.names = TRUE,
                 row.names = FALSE)
}
## Consolidation and saving summary statistics for each sub-set of data.
options(scipen = 999)
for (l in 1:3){
     s.set = sub.sets[[l]][, -1:-8]
     sum.stat = apply(s.set,
                      2,
                      function(x){x = suppressWarnings(na.omit(as.numeric(x)))
                                  x = x[x > 0]
                                  stat = c(min(x),
                                           round(mean(x), digits = 1),
                                           round(median(x), digits = 1),
                                           max(x),
                                           round(sd(x), digits = 1),
                                           round((sd(x)/mean(x)) * 100, digits = 1),
                                           round(((nrow(s.set) - length(x))/nrow(s.set)) * 100, digits = 1))})
     sum.stat = cbind(c("min", "mean", "median", "max", "sd", "cv", "overlap.freq"), sum.stat)
     colnames(sum.stat) = c("metric", colnames(sum.stat)[-1])
     file.name = paste(paste(sub(".txt", ", distances to epimarks, ", eej),
                             names(sub.sets[l]), ", summary.txt", sep = ""))
     write.table(sum.stat,
                 file = file.name,
                 sep = "\t",
                 quote = FALSE,
                 col.names = TRUE,
                 row.names = FALSE)
}
return(expEEJ)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/annotateepigeneticfeatures.r")
#   eej.file = "GVV/Experimental_junctions.txt"
#   features.file = "GVV/List of files.txt"
#   res = annotateEpigeneticFeatures(eej = eej.file, features = features.file)
