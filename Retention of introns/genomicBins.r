###################################################################################################
##  High-level R function for developing a list of non-overlapping intronic and exonic bins      ##
##  based on reference models.                                                                   ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work       - path to and name of work directory.
##  f.sqlite     - the name of SQLite database with reference annotations. It is typically SQLite
#                  version of Ensembl annotations.
##  tRNAs_genes  - the name of TXT file in tab-delimited format with genomic coordinates of the
#                  exons of human tRNAs. This file should include five mandatory fields:
#                  i) gene_id   - IDs of genes;
#                  ii) seqnames - the name of chromosome or scaffold with prefix "chr";
#                  iii) start   - start genomic coordinate of the exon;
#                  iv) end      - end genomic coordinate of the exon;
#                  v) strand    - strand information about exon location.
##  thr          - integer giving the number of threads used for running of function.
#                  Default value is 1.
##  f.mapability - the name of file with reference genome mapability data.
##  m.over       - integer giving minimal overlap within the intron for the RNA-Seq reads that
#                  spanned an exon-intron junction. Default value is 5.
genomicBins = function(d.work, f.sqlite, tRNAs_genes, thr = 1, f.mapability, m.over = 5){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries GenomicFeatures v.1.36.4, data.table v.1.12.2,
#   and pbapply v.1.4-2.
suppressMessages(require(GenomicFeatures))
suppressMessages(require(data.table))
suppressMessages(require(pbapply))
suppressMessages(require(parallel))
##  Loading of the reference annotations as a TxDb object.
TxDb = loadDb(file = paste(d.work, f.sqlite, sep = "/"))
##  Format of the seqnames in a TxDb object.
format.seqnames = unique(x = gsub("chr.+", "chr", seqnames(seqinfo(TxDb))))
##  Development a matrix of the gene-wise names of transcripts.
tr.names = unlist(x = transcriptsBy(x = TxDb, by = "gene"))
tr.names$gene_id = names(x = tr.names)
tr.names = as.matrix(x = unique(x = data.frame(tr.names, stringsAsFactors = FALSE)[, 7:8]))
rownames(tr.names) = tr.names[, 1]
##  Development a list of the reference one-exon and multi-exon transcripts.
tr.list = exonsBy(x = TxDb, by = "tx", use.names = TRUE)
tr.one.exon = tr.list[lengths(tr.list) == 1, ]
tr.one.exon = data.frame(sort(x = unique(x = unlist(x = tr.one.exon)[, 0])),
                         stringsAsFactors = FALSE)
if (!format.seqnames[1] == "chr"){
tr.one.exon$seqnames = paste("chr", tr.one.exon$seqnames, sep = "")
}
tr.one.exon = tr.one.exon[tr.one.exon$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
tr.one.exon = makeGRangesFromDataFrame(df = tr.one.exon, keep.extra.columns = TRUE)
tr.list = tr.list[lengths(tr.list) > 1, ]
tr.list = unlist(x = tr.list)[, 3]
tr.list$transcript_id = names(tr.list)
tr.list = data.frame(tr.list, stringsAsFactors = FALSE)[, c(7, 1:3, 5, 4, 6)]
if (!format.seqnames[1] == "chr"){
tr.list$seqnames = paste("chr", tr.list$seqnames, sep = "")
}
tr.list = tr.list[tr.list$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
tr.list = makeGRangesFromDataFrame(df = tr.list, keep.extra.columns = TRUE)
##  Retrieving the genomic coordinates of exons of the reference transcripts.
#   Development a list of genomic coordinates of the exons of tRNAs.
#   This list is based on the genomic coordinates of the human tRNA genes
#   deposited in the GtRNAdb database at http://gtrnadb.ucsc.edu/GtRNAdb2/index.html
exons.tRNAs = read.table(file = paste(d.work, tRNAs_genes, sep = "/"),
                         sep = "\t",
                         header = TRUE,
                         quote = "\"",
                         as.is = TRUE)[, -1]
#   Retrieving of the genomic coordinates of all exons.
exons = unique(x = data.frame(exons(x = TxDb))[, c(1:3, 5, 4)])
if (!format.seqnames[1] == "chr"){
exons$seqnames = paste("chr", exons$seqnames, sep = "")
}
exons = exons[exons$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
exons = rbind(exons, exons.tRNAs)
exons = sort(x = makeGRangesFromDataFrame(df = exons, keep.extra.columns = TRUE))
names(exons) = NULL
#   Retrieving of the genomic coordinates of first exons.
exons.first = sort(x = unique(x = tr.list[tr.list$exon_rank == 1, 0]))
#   Retrieving of the genomic coordinates of last exons.
exons.last = data.table(data.frame(tr.list))
exons.last.number = exons.last[, list(max(exon_rank), seqnames), by = transcript_id]
exons.last = unique(x = exons.last[exons.last$transcript_id == exons.last.number$transcript_id &
                                   exons.last$exon_rank == exons.last.number$V1, ][, 1:5])
exons.last = sort(x = makeGRangesFromDataFrame(df = exons.last, keep.extra.columns = TRUE))
#   Retrieving of the genomic coordinates of internal exons.
hits  = findOverlaps(query = tr.list,
                     subject = c(exons.first, exons.last),
                     type = "equal",
                     ignore.strand = FALSE)
exons.internal = sort(x = unique(x = tr.list[-queryHits(hits), 0]))
##  Retrieving of the genomic coordinates of introns of the reference transcripts.
#   Retrieving of the genomic coordinates of all introns.
introns = unlist(x = intronsByTranscript(x = TxDb, use.names = TRUE))
introns$gene_id = names(introns)
introns = as.matrix(data.frame(introns, stringsAsFactors = FALSE))
introns[introns[, 6] %in% tr.names[, 1], 6] = tr.names[introns[introns[, 6] %in%
                                                       tr.names[, 1], 6], 2]
introns = unique(x = data.frame(introns, stringsAsFactors = FALSE))[, c(6, 1:3, 5, 4)]
if (!format.seqnames[1] == "chr"){
introns$seqnames = paste("chr", introns$seqnames, sep = "")
}
introns = introns[introns$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
rownames(introns) = NULL
introns = sort(x = makeGRangesFromDataFrame(df = introns, keep.extra.columns = TRUE))
#   Development of a gene-specific list of non overlapping intronic bins based on reference models.
#   Genomic regions of small RNAs (rRNA, tRNA, miRNA, miscRNA, ribozyme, vaultRNA, sRNA, snRNA,
#   scaRNA, scRNA and snoRNA) as well as one exon non-coding and protein-coding genes is dropped
#   out from intronic intervals.
bins.introns = c(exons, introns[, 0])
bins.introns = disjoin(x = bins.introns, ignore.strand = TRUE)
hits  = findOverlaps(query = bins.introns,
                     subject = exons,
                     type = "within",
                     ignore.strand = TRUE)
bins.introns = bins.introns[-queryHits(hits), ]
hits  = findOverlaps(query = bins.introns,
                     subject = introns,
                     type = "within",
                     ignore.strand = TRUE)
bins.introns = bins.introns[queryHits(hits), ]
strand(bins.introns) = strand(introns[subjectHits(hits), ])
bins.introns$gene_id = introns[subjectHits(hits), ]$gene_id
bins.introns$intron_type = "CI"
bins.introns = bins.introns[width(bins.introns) >= 50, ]
bins.introns = unique(x = data.frame(bins.introns))
bins.introns = sort(x = makeGRangesFromDataFrame(df = bins.introns, keep.extra.columns = TRUE))
names(bins.introns) = NULL
#   Retrieving of the genomic coordinates of retained introns.
hits  = findOverlaps(query = introns, subject = exons, type = "within", ignore.strand = FALSE)
r.introns = introns[queryHits(hits), ]
r.introns = makeGRangesFromDataFrame(df = unique(x = data.frame(r.introns)),
                                     keep.extra.columns = TRUE)
names(r.introns) = NULL
cl = makeCluster(thr)
clusterExport(cl, c("findOverlaps", "queryHits", "subjectHits", "disjoin", "exons"))
pboptions(type = "timer", style = 1, char = ".")
r.introns = pblapply(X = r.introns,
                     FUN = function(x){hits  = findOverlaps(query = x,
                                                            subject = exons,
                                                            type = "any",
                                                            ignore.strand = TRUE)
                                       r.exons = exons[subjectHits(hits), ]
                                       hits  = findOverlaps(query = x,
                                                            subject = r.exons,
                                                            type = "within",
                                                            ignore.strand = FALSE)
                                       r.exons = r.exons[-subjectHits(hits), ]
                                       if (length(x = r.exons) > 0){
                                           hits  = findOverlaps(query = x,
                                                                subject = r.exons,
                                                                type = "within",
                                                                ignore.strand = TRUE)
                                           if (length(x = hits) > 0){
                                               x = x[-queryHits(hits), ]
                                           }else{dj = disjoin(x = c(x[, 0], r.exons),
                                                              ignore.strand = TRUE)
                                                 hits  = findOverlaps(query = dj,
                                                                      subject = r.exons,
                                                                      type = "within",
                                                                      ignore.strand = TRUE)
                                                 dj = dj[-queryHits(hits), ]
                                                 if (length(x = dj) > 0){
                                                     dj@strand@values = x@strand@values
                                                     dj$gene_id = x$gene_id
                                                 }
                                                 x = dj
                                           }
                                       }
                                       return(x)
                                       },
                     cl = cl
                    )
stopCluster(cl)
#   Retrieving of the genomic coordinates of bins of retained introns.
bins.rintrons = unique(x = data.frame(do.call(what = getMethod(c, "GenomicRanges"),
                                              args = r.introns)))
bins.rintrons = bins.rintrons[bins.rintrons$width >= 50, ]
bins.rintrons = sort(x = makeGRangesFromDataFrame(df = bins.rintrons, keep.extra.columns = TRUE))
names(bins.rintrons) = NULL
bins.rintrons$intron_type = "RI"
hits = findOverlaps(query = bins.rintrons,
                    subject = tr.one.exon,
                    type = "within",
                    ignore.strand = FALSE)
bins.rintrons[queryHits(hits), ]$intron_type = "RI_OneExonTranscript"
hits = findOverlaps(query = bins.rintrons,
                    subject = exons.internal,
                    type = "within",
                    ignore.strand = FALSE)
bins.rintrons[queryHits(hits), ]$intron_type = "RI_INTERNAL"
hits = findOverlaps(query = bins.rintrons,
                    subject = exons.first,
                    type = "within",
                    ignore.strand = FALSE)
bins.rintrons[queryHits(hits), ]$intron_type = "RI_FIRST"
hits = findOverlaps(query = bins.rintrons,
                    subject = exons.last,
                    type = "within",
                    ignore.strand = FALSE)
bins.rintrons[queryHits(hits), ]$intron_type = "RI_LAST"
##  Development a full list of reference intronic bins.
bins.introns = c(bins.introns, bins.rintrons)
bins.introns = data.frame(bins.introns, stringsAsFactors = FALSE)
bins.introns = bins.introns[order(bins.introns$gene_id, bins.introns$start, bins.introns$end), ]
bins.introns$intron_id = paste(bins.introns$gene_id,
                               paste("INTR",
                                     formatC(x = seq(1, nrow(bins.introns), 1),
                                             format = "d",
                                             width = 6,
                                             flag = "0"),
                                     sep = ""),
                               sep = ":")
bins.introns = bins.introns[, c(8, 7, 6, 1:3, 5, 4)]
rownames(bins.introns) = NULL
##  Correction of effective length of intronic bins.
#   Loading of reference genome mapability data.
mapability = read.table(file = paste(d.work, f.mapability, sep = "/"),
                        sep = "\t",
                        header = TRUE,
                        quote = "\"",
                        as.is = TRUE)
#   Transformation of mapability data in a GRanges object.
mapability = makeGRangesFromDataFrame(df = mapability, keep.extra.columns = TRUE)
#   Calculation of effective length of introns.
introns.1 = makeGRangesFromDataFrame(df = bins.introns, keep.extra.columns = TRUE)
start(introns.1) = start(introns.1) + m.over - 1
end(introns.1) = end(introns.1) - m.over + 1
exclGR = disjoin(x = c(mapability, introns.1[, 0]), ignore.strand = TRUE)
hits = findOverlaps(query = exclGR, subject = introns.1, type = "within", ignore.strand = TRUE)
exclGR = exclGR[queryHits(hits), ]
hits = findOverlaps(query = exclGR, subject = mapability, type = "within", ignore.strand = TRUE)
exclGR = exclGR[queryHits(hits), ]
exclGR = exclGR[!duplicated(exclGR), ]
hits = findOverlaps(query = exclGR, subject = introns.1, type = "within", ignore.strand = TRUE)
introns.2 = introns.1[subjectHits(hits), ]
exclGR = exclGR[queryHits(hits), ]
introns.2$effect_length = width(exclGR)
introns.2 = data.frame(introns.2, stringsAsFactors = FALSE)
introns.2 = aggregate(x = introns.2[, 9], by = list(introns.2[, 6]), FUN = sum)
rownames(introns.2) = introns.2[, 1]
introns.2 = as.matrix(x = introns.2)
bins.introns$effective_length = 0
bins.introns = as.matrix(x = bins.introns)
index = bins.introns[, 1] %in% introns.2[, 1]
bins.introns[index, 9] = introns.2[bins.introns[index, 1], 2]
bins.introns = data.frame(bins.introns, stringsAsFactors = FALSE)
bins.introns[, c(5:6, 8:9)] = apply(X = bins.introns[, c(5:6, 8:9)],
                                    MARGIN = 2,
                                    FUN = function(x){as.numeric(x)})
bins.introns$effective_length = width(introns.1) - bins.introns$effective_length
##  Development a full list of reference exonic bins.
bins.exons = exonicParts(txdb = TxDb, linked.to.single.gene.only = TRUE)[, 3]
bins.exons = unique(x = data.frame(bins.exons, stringsAsFactors = FALSE))
if (!format.seqnames[1] == "chr"){
bins.exons$seqnames = paste("chr", bins.exons$seqnames, sep = "")
}
bins.exons = bins.exons[bins.exons$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
bins.exons = sort(x = makeGRangesFromDataFrame(df = bins.exons, keep.extra.columns = TRUE))
hits = findOverlaps(query = bins.exons,
                    subject = bins.rintrons,
                    type = "equal",
                    ignore.strand = TRUE)
bins.exons = bins.exons[-queryHits(hits), ]
hits = findOverlaps(query = bins.exons,
                    subject = bins.rintrons,
                    type = "any",
                    ignore.strand = TRUE)
bins.exons.1 = bins.exons[queryHits(hits), ]
bins.exons = bins.exons[-queryHits(hits), ]
exclGR = disjoin(x = c(bins.exons.1[, 0], bins.rintrons[, 0]), ignore.strand = TRUE)
hits = findOverlaps(query = exclGR,
                    subject = bins.rintrons,
                    type = "within",
                    ignore.strand = TRUE)
exclGR = unique(x = exclGR[-queryHits(hits), ])
hits = findOverlaps(query = exclGR,
                    subject = bins.exons.1,
                    type = "within",
                    ignore.strand = TRUE)
exclGR = exclGR[queryHits(hits), ]
bins.exons.1 = bins.exons.1[subjectHits(hits), ]
strand(exclGR) = strand(bins.exons.1)
exclGR$gene_id = bins.exons.1$gene_id
bins.exons = unique(x = data.frame(c(bins.exons, exclGR), stringsAsFactors = FALSE))
bins.exons = bins.exons[bins.exons$width >= 10, ]
bins.exons = bins.exons[order(bins.exons$gene_id, bins.exons$start, bins.exons$end), ]
bins.exons$exon_id = paste(bins.exons$gene_id,
                           paste("EXON",
                                 formatC(x = seq(1, nrow(bins.exons), 1),
                                         format = "d",
                                         width = 6,
                                         flag = "0"),
                                 sep = ""),
                           sep = ":")
bins.exons = bins.exons[, c(7, 6, 1:3, 5, 4)]
rownames(bins.exons) = NULL
##  Final object to be returned.
return(list(intronic.bins = bins.introns, exonic.bins = bins.exons))
}
