###################################################################################################
##  High-level R function for developing a list of non-overlapping exonic bins                   ##
##  based on reference models.                                                                   ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work   - path to and name of work directory.
##  f.sqlite - the name of SQLite database with reference annotations. It is typically SQLite
#              version of Ensembl annotations.
##  can.chr  - logical, if TRUE (default) then only exonic bins
#              belonging to the canonical chromosomes will be processed and saved.
##  exonL    - integer value giving the minimal length of exonic bins. Default value is 10.
exonicBins = function(d.work, f.sqlite, can.chr = TRUE, exonL = 10){
##  Loading of required auxiliary library.
#   This code was successfully tested with library GenomicFeatures v.1.36.4.
suppressMessages(require(GenomicFeatures))
##  Loading of the reference annotations as a TxDb object.
TxDb = loadDb(file = paste(d.work, f.sqlite, sep = "/"))
##  Development a full list of reference exonic bins.
bins.exons = exonicParts(txdb = TxDb, linked.to.single.gene.only = TRUE)[, 3]
bins.exons = unique(x = data.frame(bins.exons, stringsAsFactors = FALSE))
##  Only canonical chromosomes.
if (can.chr == "TRUE"){
format.seqnames = unique(x = gsub("chr.+", "chr", seqnames(seqinfo(TxDb))))
    if (!format.seqnames[1] == "chr"){
        bins.exons$seqnames = paste("chr", bins.exons$seqnames, sep = "")
    }
    bins.exons = bins.exons[bins.exons$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
}
bins.exons = bins.exons[bins.exons$width >= exonL, ]
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
return(bins.exons)
}
