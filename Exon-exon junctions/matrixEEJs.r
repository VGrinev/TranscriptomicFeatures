###################################################################################################
##  High-level R function for development of count matrix of exon-exon junctions                 ##
##  from subjunc based BED files.                                                                ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work - path to and name of work directory.
##  d.bed - the name of directory with BED files.
##  f.pmEEJs - the name of file to save primary count matrix. This matrix contains count data
#   on exon-exon junctions not mapped to genes. Default value is FALSE, and no saving performig.
##  f.genes - the name of TXT file in tab-delimited format with genomic coordinates of genes.
#   This file should include five mandatory fields:
#   i) gene_id - IDs of genes;
#   ii) seqnames - the name of chromosome or scaffold with prefix "chr";
#   iii) start - start genomic coordinate of the gene;
#   iv) end - end genomic coordinate of the gene;
#   v) strand  - strand information about gene location.
##  f.grEEJs - the name of file to save finnal count matrix of exon-exon junctions mapped to genes.
### Important notes. This code is adapted to process BED files produced using subjunc function
#   from R/Bioconductor package Rsubread. In finnal count matrix, the start coordinate of exon-exon
#   junction correspond to last 3'-end nucleotide of upstream exon and the end coordinate is first
#   5'-end nucleotide of downstream exon.
matrixEEJs = function(d.work, d.bed, f.pmEEJs = FALSE, f.genes, f.grEEJs){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries data.table v.1.12.2 and GenomicFeatures v.1.36.4.
suppressMessages(require(data.table))
suppressMessages(require(GenomicFeatures))
##  Listing of BED file(-s).
bed.files = list.files(path = paste(d.work, d.bed, sep = "/"), pattern = "[.]bed$")
##  Loading of BED file(-s) and extraction of exon-exon junctions.
list.EEJs = list()
for (i in 1:length(bed.files)){
#   Loading of BED file.
bed.file = read.table(file = paste(paste(d.work, d.bed, sep = "/"), bed.files[i], sep = "/"),
                      sep = "\t",
                      as.is = TRUE)
#   Only canonical chromosomes.
chr = paste("chr", c(1:22, "X", "Y"), sep = "")
bed.file = bed.file[bed.file$V1 %in% chr, ]
#   Extraction of exon-exon junctions.
EEJs = data.frame(cbind(bed.file$V1,
                        bed.file$V2 + as.numeric(x = do.call(what = rbind, 
                                                             args = strsplit(x = bed.file$V11,
                                                                             split = ","))[, 1]),
                        bed.file$V3 - as.numeric(x = do.call(what = rbind,
                                                             args = strsplit(x = bed.file$V11,
                                                                             split = ","))[, 2]) + 1,
                        bed.file$V6,
                        bed.file$V5), stringsAsFactors = FALSE)
colnames(EEJs) = c("seqnames", "start", "end", "strand", "coverage")
EEJs[, c(2, 3, 5)] = apply(X = EEJs[, c(2, 3, 5)], MARGIN = 2, FUN = as.numeric)
EEJs$eej_id = paste(paste(EEJs$seqnames,
                          paste(EEJs$start, EEJs$end, sep = "-"), sep = ":"),
                    EEJs$strand, sep = "_str")
EEJs = EEJs[order(EEJs$eej_id), ]
#   Correction of coverage for exon-exon junctions.
if (length(x = unique(x = EEJs$eej_id)) < length(x = EEJs$eej_id)){
    coverage = data.table(EEJs[, c(6, 5)])
    coverage = coverage[, list(sum(coverage)), by = c("eej_id")]
    coverage = data.frame(coverage)
    coverage = coverage[order(coverage$eej_id), ]
    EEJs = EEJs[, -5]
    EEJs = EEJs[!duplicated(EEJs), ]
    EEJs = EEJs[order(EEJs$eej_id), ]
    EEJs$coverage = coverage$V1
    EEJs = EEJs[, c(1:4, 6, 5)]
   }
EEJs = EEJs[, c(6, 1:4, 5)]
list.EEJs[[i]] = EEJs
if (i%%1 == 0){
    cat(i, "BED file of", length(x = bed.files), "has been processed...\n")
    flush.console()
   }
}
##  Development of primary count matrix of exon-exon junctions.
pmEEJs = do.call(what = rbind, args = list.EEJs)[, -6]
pmEEJs = as.matrix(x = pmEEJs[!duplicated(pmEEJs), ])
pmEEJs = cbind(pmEEJs, matrix(0L, ncol = length(x = bed.files), nrow = nrow(pmEEJs)))
for (k in 1:length(x = bed.files)){
     cover.values = list.EEJs[[k]][, c(1, 6)]
     rownames(cover.values) = cover.values[, 1]
     cover.values = as.matrix(x = cover.values)
     index = pmEEJs[, 1] %in% cover.values[, 1]
     pmEEJs[index, 5 + k] = cover.values[pmEEJs[index, 1], 2]
}
colnames(pmEEJs) = c(colnames(pmEEJs)[1:5],
                     gsub(pattern = ".bed", replacement = "", x = bed.files))
pmEEJs = data.frame(pmEEJs, stringsAsFactors = FALSE)
pmEEJs[, c(3:4, 6:(5 + length(x = bed.files)))] = apply(X = pmEEJs[, c(3:4, 6:(5 + length(x = bed.files)))],
                                                        MARGIN = 2,
                                                        FUN = as.numeric)
if (f.pmEEJs != FALSE){
write.table(x = pmEEJs,
            file = paste(d.work, f.pmEEJs, sep = "/"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
    }
##  Gene mapping of exon-exon junctions.
genes = makeGRangesFromDataFrame(df = read.table(file = paste(d.work, f.genes, sep = "/"),
                                                 sep = "\t",
                                                 header = TRUE,
                                                 quote = "\"",
                                                 as.is = TRUE),
                                keep.extra.columns = TRUE)
grEEJs = makeGRangesFromDataFrame(df = pmEEJs,
                                  keep.extra.columns = TRUE)
hits = findOverlaps(query = grEEJs, subject = genes, type = c("within"), ignore.strand = FALSE)
grEEJs = grEEJs[queryHits(hits), ]
grEEJs$gene_id = genes[subjectHits(hits), ]$gene_id
names(grEEJs) = NULL
grEEJs = as.data.frame(grEEJs)
grEEJs = grEEJs[, c(6, ncol(grEEJs), 1:5, 7:(ncol(grEEJs) - 1))]
grEEJs[, c(4:6, 8:(7 + length(x = bed.files)))] = apply(X = grEEJs[, c(4:6, 8:(7 + length(x = bed.files)))],
                                                        MARGIN = 2,
                                                        FUN = as.numeric)
write.table(x = grEEJs,
            file = paste(d.work, f.grEEJs, sep = "/"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
return(list(pmEEJs = pmEEJs, grEEJs = grEEJs))
}
