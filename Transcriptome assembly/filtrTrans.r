###################################################################################################
##  High-level R function                                                                        ##
##  for filtration of the Cufflinks assembled transcripts.                                       ##
##  (c) GNU GPL Vasily V. Grinev, 2016-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work - path to and name of work directory.
##  in.gtf - name of the input file in GTF/GFF format. It is typically output file of Cufflinks
#   (or similar ones).
##  unstr - a logical argument indicating that the unstranded transcripts should be removed.
#   Default value is TRUE.
##  difStr - a logical argument specifying that the transcripts that match two different strands
#   should be removed. Default value is TRUE.
##  canChr - a logical argument specifying that the records from non-canonical chromosomes should
#   be removed. Default value is TRUE.
##  oneEx - a logical argument indicating that the one-exon transcripts should be removed.
#   Default value is FALSE.
##  exL - an integer argument. It is a threshold for minimal length of exon(-s) in transcripts.
#   Default Ensembl-based value is 25.
##  fpkm - an integer argument. It is a threshold for minimal transcripts abundance (in FPKM).
#   Default value is 1.
##  isoforms.fpkm - name of the file with FPKM data on RNA isoforms expression.
#   It is typically *.fpkm_tracking output file of Cufflinks (or similar ones).
##  trL - an integer argument. It is a threshold for minimal length of transcripts.
#   Default value is 300.
##  intrL - an integer argument. It is a threshold for minimal length of intron(-s) in transcripts.
#   Default Ensembl-based value is 50.
##  out.gtf - name of the output file in GTF/GFF format.
##  format - format of the output file (typically "gtf" or "gff").
filtrTrans = function(d.work, in.gtf, unstr = TRUE, difStr = TRUE, canChr = TRUE, oneEx = FALSE,
                      exL = 25, fpkm = 1, isoforms.fpkm, trL = 300, intrL = 50, out.gtf, format){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with library rtracklayer v.1.44.3
#   and library Rsubread v.1.34.6
suppressMessages(library(rtracklayer))
cat("Required auxiliary library has been loaded...\n")
flush.console()
##  Loading of the primary data as a GRanges object.
tr = import(con = paste(d.work, in.gtf, sep = "/"), format = "gtf")
cat("Primary data have been loaded...\n")
flush.console()
##  Removing of unstranded transcripts.
if (unstr == TRUE){
    tr = tr[!strand(tr) == "*", ]
    cat("Unstranded transcripts have been removed...\n")
    flush.console()
}
##  Removing of transcripts that match two different strands.
if (difStr == TRUE){
    d.str = cbind(as.character(strand(tr)), tr$transcript_id)
    d.str = table(d.str[!duplicated(d.str), ][, 2])
    tr = tr[!tr$transcript_id %in% attr(x = d.str[d.str > 1], which = "dimnames")[[1]], ]
    cat("Transcripts that match two different strands have been removed...\n")
    flush.console()
}
##  Only canonical chromosomes.
if (canChr == TRUE){
    tr = as.data.frame(tr)
    chr = c(paste("chr", c(1:22, "X", "Y"), sep = ""))
    tr = tr[tr$seqnames %in% chr, ][, -8:-9]
    rownames(tr) = NULL
    tr = makeGRangesFromDataFrame(tr, keep.extra.columns = TRUE)
    cat("Records from non-canonical chromosomes have been removed...\n")
    flush.console()
}
##  Removing of one-exon transcripts.
if (oneEx == TRUE){
    tr = tr[tr$transcript_id %in%
            attr(x = table(tr$transcript_id)[table(tr$transcript_id) >
                     length(x = unique(x = tr$type))],
                 which = "dimnames")[[1]], ]
    cat("One-exon transcripts have been removed...\n")
    flush.console()
}
##  Removing of transcripts with too short exon(-s).
if (exL > 0){
    tr = tr[!tr$transcript_id %in% unique(x = tr[width(tr) < exL]$transcript_id), ]
    cat("Transcripts with too short exon(-s) have been removed...\n")
    flush.console()
}
##  Removing of transcripts with too low abundance.
if (fpkm > 0){
    fpkm.data = read.table(file = paste(d.work, isoforms.fpkm, sep = "/"),
                           sep = "\t",
                           header = TRUE,
                           quote = "\"",
                           as.is = TRUE)
    fpkm.data = fpkm.data[apply(X = fpkm.data[, grep(pattern = "FPKM",
                                                     x = colnames(x = fpkm.data),
                                                     ignore.case = TRUE)],
                                MARGIN = 1,
                                FUN = max) >= fpkm, ]$tracking_id
    tr = tr[tr$transcript_id %in% fpkm.data, ]
    cat("Transcripts with too low abundance have been removed...\n")
    flush.console()
}
##  Removing too short transcripts.
if (trL > 0){
    tr.length = aggregate(x = width(tr), by = list(tr$transcript_id), FUN = sum)
    tr = tr[tr$transcript_id %in% tr.length[tr.length$x >= trL, ]$Group.1, ]
    cat("Too short transcripts have been removed...\n")
    flush.console()
}
##  Removing of transcripts with too short intron(-s).
if (intrL > 0){
    introns = split(x = tr, f = tr$transcript_id)
    introns = unlist(x = setdiff(x = range(introns), y = introns))
    tr = tr[tr$transcript_id %in% unique(x = names(x = introns[width(introns) >= intrL, ])), ]
    cat("Transcripts with too short intron(-s) have been removed...\n")
    flush.console()
}
##  Saving of the final GRanges object as GTF/GFF file.
    cat("Saving of the final GRanges object as GTF/GFF file...\n")
    flush.console()
    export(object = tr, con = paste(d.work, out.gtf, sep = "/"), format = format)
}
