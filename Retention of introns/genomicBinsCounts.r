###################################################################################################
##  High-level R function for read summarization across of intronic and/or exonic bins.          ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work - path to and name of work directory.
##  f.gr   - the name of TXT file in tab-delimited format with genomic intervals of interest.
#            This file should include seven (or nine for intronic bins) mandatory fields:
#            i) intron_id, or exon_id - IDs of genomic bins (intronic or exonic, respectively);
#            ii) intron_type          - type of introns (as returned by genomicBins function).
#                                       This field is required only for a file with intronic bins;
#            iii) gene_id             - IDs of genes;
#            iv) seqnames             - the name of chromosome or scaffold with prefix "chr";
#            v) start                 - start coordinate of genomic bin;
#            vi) end                  - end coordinate of genomic bin;
#            vii) strand              - strand information about bin location;
#            viii) width              - width of genomic bin;
#            ix) effective_length     - effective length of intronic bins as calculated by
#                                       genomicBins function. This field is required only for
#                                       a file with intronic bins.
##  d.bam  - the name of BAM directory.
##  f.bam  - a character vector with list of BAM files to be treated. Default value is NULL.
#            If so, all BAM files in BAM directory will be processed.
##  m.over - an integer value giving minimal overlap within the genomic bin
#            for the RNA-Seq read to be summarized.
##  stra   - an integer value indicating if strand-specific read counting should be performed.
#            It can be following three values: 0 (unstranded), 1 (stranded) and
#            2 (reversely stranded). Default value of this parameter is 0.
##  paired - logical indicating if paired-end reads are used. If TRUE, fragments will be counted
#            instead of individual reads. FALSE by default.
##  thr    - integer giving the number of threads used for running of featureCounts function.
#            Default value is 1.
##  type   - type of genomic bins. There are two values: "intronic" or "exonic".
genomicBinsCounts = function(d.work, f.gr, d.bam, f.bam = NULL,
                             m.over = 1, stra = 0, paired = FALSE, thr = 1, type){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with
#   libraries GenomicFeatures v.1.28.5 and Rsubread v.1.26.1.
suppressMessages(require(GenomicFeatures))
suppressMessages(require(Rsubread))
##  Loading of the genomic intervals of interest as a GRanges object.
gr = read.table(file = paste(d.work, f.gr, sep = "/"),
                sep = "\t",
                header = TRUE,
                quote = "\"",
                as.is = TRUE)
if (type == "intronic"){ 
colnames(gr)[3] = "id"
}else{
colnames(gr)[2] = "id"
}
gr = makeGRangesFromDataFrame(df = gr, keep.extra.columns = TRUE)
##  Creation of annotation object.
anno = createAnnotationFile(GR = gr)
##  List of BAM files to be processed.
if (is.null(f.bam)){
f.bam = list.files(path = paste(d.work, d.bam, sep = "/"), pattern = "[.]bam$")
}
##  Read summarization.
counts = featureCounts(files = paste(paste(d.work, d.bam, sep = "/"), f.bam, sep = "/"),
                       annot.ext = anno, 
                       useMetaFeatures = FALSE,
                       allowMultiOverlap = TRUE,
                       minOverlap = m.over,
                       strandSpecific = stra,
                       isPairedEnd = paired,
                       checkFragLength = TRUE,
                       countChimericFragments = FALSE,
                       nthreads = thr)
##  Consolidation results.
if (type == "intronic"){ 
counts = suppressWarnings(cbind(elementMetadata(gr)@listData, counts[[2]][, -1], counts[[1]]))
counts = counts[, c(1:3, 5:9, 4, 10:ncol(counts))]
colnames(counts) = c("intron_id",
                     "intron_type",
                     "gene_id",
                     "seqnames",
                     "start",
                     "end",
                     "strand",
                     "width",
                     "effective_length",
                     sub(".bam", "", f.bam))
}else{
counts = suppressWarnings(cbind(elementMetadata(gr)@listData, counts[[2]][, -1], counts[[1]]))
colnames(counts) = c("exon_id",
                     "gene_id",
                     "seqnames",
                     "start",
                     "end",
                     "strand",
                     "width",
                     sub(".bam", "", f.bam))
}
rownames(counts) = NULL
return(counts)
}
### A simple example of function use.
res = genomicBinsCounts(d.work = "/home/hmglab/GVV",
                        f.gr = "Ensembl based exonic bins, GRCh38.p7, release 85.txt",
                        d.bam = "Files_BAM",
                        f.bam = NULL,
                        m.over = 1,
                        stra = 0,
                        paired = FALSE,
                        thr = 8,
                        type = "exonic")
d.work = "/home/hmglab/GVV"
f.res = "The NML project, exonic bins, primary matrix of read counts.txt"
write.table(res,
            file = paste(d.work, f.res, sep = "/"),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
