###################################################################################################
##  High-level R function for calculation of hypothetical "non-alternative"                      ##
##  precursor RNA based on reference annotations of RNAs for gene(-s) of interest.               ##
##  (c) GNU GPL Vasily V. Grinev, 2017-2020. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work - character string giving the path to and name of work directory.
##  f.ref  - character string giving the name of file in GTF/GFF format with reference annotations
#            of RNAs for gene(-s) of interest. It is typically Ensembl annotations or similar one.
#   format - format of file with reference annotations. Default value is "gtf".
##  src    - a character vector with name of source of reference annotations.
#            Default value is "Ensembl".
hnapRNA = function(d.work, f.ref, format = "gtf", src = "Ensembl"){
### Loading of required auxiliary libraries.
#   This code was successfully tested with libraries GenomicFeatures v.1.36.4
#   and rtracklayer v.1.44.3.
suppressMessages(require(GenomicFeatures))
suppressMessages(require(rtracklayer))
### Loading of reference annotations as an object of class GRanges.
ref = import(con = paste(d.work, f.ref, sep = "/"), format = format)
### Filtration of reference annotations.
ref = sort(x = ref[ref$type == "exon", c("gene_id", "transcript_id", "exon_number")])
ref$exon_number = as.numeric(ref$exon_number)
one.exon.trans = table(ref$transcript_id)
one.exon.trans = attr(one.exon.trans[one.exon.trans == 1], "dimnames")[[1]]
ref = ref[!ref$transcript_id %in% one.exon.trans, ]
### Development a list of genes.
list.genes = makeGRangesListFromDataFrame(df = as.data.frame(ref),
                                          split.field = "gene_id",
                                          keep.extra.columns = TRUE)
### Modelling of each individual gene.
model = lapply(X = list.genes, FUN = function(x){
##  Development a list of transcripts.
list.trans = makeGRangesListFromDataFrame(df = as.data.frame(x),
                                          split.field = "transcript_id",
                                          keep.extra.columns = TRUE)
if (length(list.trans) > 1){
##  Retrieving the genomic coordinates of all unique introns.
intr.all = sort(x = unique(x = unlist(x = setdiff(x = range(list.trans), list.trans))))
intr.all = intr.all[width(intr.all) >= 50, ]
names(intr.all) = NULL
##  Retrieving the genomic coordinates of exons with retained introns.
ex.ri = x[subjectHits(findOverlaps(intr.all, x, type = "within")), 0]
##  Retrieving the genomic coordinates of internal exons.
ex.internal = GRangesList(lapply(X = list.trans,
                                 FUN = function(y){y[-c(1, length(y)), 0]}))
#   Development the final list of internal exons.
ex.internal = sort(x = unique(x = unlist(x = ex.internal)))
names(ex.internal) = NULL
if (length(ex.ri) > 0){
ex.internal = ex.internal[!ex.internal %in%
                          subsetByOverlaps(ex.internal, ex.ri, type = "equal"), ]
}
##  Retrieving the genomic coordinates of first exon and alternative first exons.
ex.first = sort(x = unique(x = x[x$exon_number == 1, 0]))
if (length(ex.ri) > 0){
ex.first = ex.first[!ex.first %in%
                    subsetByOverlaps(ex.first, ex.ri, type = "equal"), ]
}
ex.first = ex.first[!ex.first %in%
                    subsetByOverlaps(ex.first, ex.internal, type = "equal"), ]
ex.first = intersect(x = range(ex.first), ex.first)
if (as.character(strand(ex.first))[1] == "+"){
    ex.first_alter = ex.first[-1, ]
    ex.first = ex.first[1, ]
}else{
    ex.first_alter = ex.first[-length(ex.first), ]
    ex.first = ex.first[length(ex.first), ]
}
##  Retrieving the genomic coordinates of last exon and alternative last exons.
ex.last = GRangesList(lapply(X = list.trans,
                             FUN = function(y){y[y$exon_number == max(y$exon_number), 0]}))
ex.last = sort(x = unique(x = unlist(x = ex.last)))
names(ex.last) = NULL
if (length(ex.ri) > 0){
ex.last = ex.last[!ex.last %in%
                  subsetByOverlaps(ex.last, ex.ri, type = "equal"), ]
}
ex.last = ex.last[!ex.last %in%
                  subsetByOverlaps(ex.last, ex.internal, type = "equal"), ]
ex.last = intersect(x = range(ex.last), ex.last)
if (as.character(strand(ex.last))[1] == "+"){
    ex.last_alter = ex.last[-length(ex.last), ]
    ex.last = ex.last[length(ex.last), ]
}else{
    ex.last_alter = ex.last[-1, ]
    ex.last = ex.last[1, ]
}
##  Development a list of all model exons.
ex.all = c(ex.first, ex.internal, ex.last)
ex.all = intersect(x = range(ex.all), ex.all)
ex.all = ex.all[width(ex.all) >= 14, ]
##  Development a list of all model introns.
intr.all = setdiff(x = range(ex.all), ex.all)
if (min(width(intr.all)) < 50){
    ex.all = c(ex.all, intr.all[width(intr.all) < 50, ])
    ex.all = intersect(x = range(ex.all), ex.all)
    intr.all = setdiff(x = range(ex.all), ex.all)
}
##  Consolidation of results.
transcript = GRanges(seqnames = seqnames(ex.all)@values,
                     ranges = IRanges(start = start(ex.all)[1],
                                      end = end(ex.all)[length(ex.all)]),
                     strand = strand(ex.all)@values)
transcript$type = "transcript"
transcript$number = NA
ex.first$type = "first_exon"
ex.first$number = NA
if (length(ex.first_alter) > 0){
    ex.first_alter$type = "alternative_first_exon"
    if (as.character(strand(ex.first_alter))[1] == "+"){
        ex.first_alter$number = c(1:length(ex.first_alter))
   }else{
        ex.first_alter$number = c(length(ex.first_alter):1)
   }
}
ex.all$type = "exon"
if (as.character(strand(ex.all))[1] == "+"){
    ex.all$number = c(1:length(ex.all))
}else{
    ex.all$number = c(length(ex.all):1)
}
intr.all$type = "intron"
if (as.character(strand(intr.all))[1] == "+"){
    intr.all$number = c(1:length(intr.all))
}else{
    intr.all$number = c(length(intr.all):1)
}
if (length(ex.last_alter) > 0){
    ex.last_alter$type = "alternative_last_exon"
    if (as.character(strand(ex.last_alter))[1] == "+"){
        ex.last_alter$number = c(1:length(ex.last_alter))
   }else{
        ex.last_alter$number = c(length(ex.last_alter):1)
   }
}
ex.last$type = "last_exon"
ex.last$number = NA
hnap.RNA = c(transcript,
             if (as.character(strand(ex.all))[1] == "+"){
                 c(ex.first,
                   sort(x = c(ex.first_alter, ex.all, intr.all, ex.last_alter)),
                   ex.last)
            }else{
                 c(ex.last,
                   sort(x = c(ex.first_alter, ex.all, intr.all, ex.last_alter)),
                   ex.first)
            })
hnap.RNA$gene_id = ""
hnap.RNA$source = as.factor(src)
hnap.RNA$transcript_id = ""
}else{
##  Processing of one-transcript genes.
transcript = GRanges(seqnames = seqnames(x)@values,
                     ranges = IRanges(start = start(x)[1],
                                      end = end(x)[length(x)]),
                     strand = strand(x)@values)
transcript$type = "transcript"
transcript$number = NA
if (as.character(strand(x)@values) == "+"){
    ex.first = x[1, 0]
}else{
    ex.first = x[length(x), 0]
}
    ex.first$type = "first_exon"
    ex.first$number = NA
    ex.all = x[, 0]
    ex.all$type = "exon"
if (as.character(strand(x)@values) == "+"){
    ex.all$number = c(1:length(ex.all))
}else{
    ex.all$number = c(length(ex.all):1)
}
    intr.all = setdiff(x = range(ex.all), ex.all)
    intr.all$type = "intron"
if (as.character(strand(x)@values) == "+"){
    intr.all$number = c(1:length(intr.all))
}else{
    intr.all$number = c(length(intr.all):1)
}
if (as.character(strand(x)@values) == "+"){
    ex.last = x[length(x), 0]
}else{
    ex.last = x[1, 0]}
    ex.last$type = "last_exon"
    ex.last$number = NA
    hnap.RNA = c(transcript,
                 if (as.character(strand(x)@values) == "+"){
                     c(ex.first,
                       sort(x = c(ex.all, intr.all)), ex.last)
                 }else{
                     c(ex.last,
                       sort(x = c(ex.all, intr.all)), ex.first)})
    hnap.RNA$gene_id = ""
    hnap.RNA$source = as.factor(src)
    hnap.RNA$transcript_id = ""
}
return(hnap.RNA)})
### Development a final object of class GRanges
#   with hypothetical "non-alternative" precursor RNA(-s).
hnap.RNA = unlist(x = as(object = model, Class = "GRangesList"))
hnap.RNA$gene_id = names(hnap.RNA)
hnap.RNA$transcript_id = paste(hnap.RNA$gene_id, "preRNA", sep = ".")
names(hnap.RNA) = NULL
return(hnap.RNA)
}
