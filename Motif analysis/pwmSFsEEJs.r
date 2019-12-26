###################################################################################################
##  High-level R function for positional weight matix based identification                       ##
##  of splicing factors binding sites surrounding exon-exon junctions.                           ##
##  (c) GNU GPL Vasily V. Grinev, 2017-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work  - character string giving the path to and name of work directory.
##  f.eej   - character string giving the name of TXT file in tab-delimited format containing
#             the exon-exon junctions. This file should include four mandatory fields:
#             i)   eej_id (ID of the exon-exon junction);
#             ii)  seqnames (name of chromosome or scaffold with prefix "chr");
#             iii) start (start genomic coordinate of the exon-exon junction);
#             iv)  end (end genomic coordinate of the exon-exon junction);
#             v)   strand (strand information about exon-exon junction).
##  f.motif - character string giving the name of TXT file in tab-delimited format containing
#             positional weight matix.
pwmSFsEEJs = function(d.work, f.eej, f.motif){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries
#   GenomicFeatures v.1.36.4 and BSgenome.Hsapiens.UCSC.hg38 v.1.4.1.
suppressMessages(require(GenomicFeatures))
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))
##  Loading of exon-exon junctions as a standard object of class "data.frame".
eej = read.table(file = paste(d.work, f.eej, sep = "/"),
                 sep = "\t",
                 header = TRUE,
                 quote = "\"",
                 as.is = TRUE)
##  Sub-setting of 5' splice sites of the exon-exon junctions.
five.ss = eej
five.ss[five.ss$strand == "+", ]$start = five.ss[five.ss$strand == "+", ]$start - 99
five.ss[five.ss$strand == "+", ]$end = five.ss[five.ss$strand == "+", ]$start + 199
five.ss[five.ss$strand == "-", ]$start = five.ss[five.ss$strand == "-", ]$end - 100
five.ss[five.ss$strand == "-", ]$end = five.ss[five.ss$strand == "-", ]$start + 199
five.ss = makeGRangesFromDataFrame(df = five.ss, keep.extra.columns = TRUE)
five.ss_seq = getSeq(x = Hsapiens, five.ss)
##  Sub-setting of 5' splice sites of the exon-exon junctions.
three.ss = eej
three.ss[three.ss$strand == "+", ]$start = three.ss[three.ss$strand == "+", ]$end - 100
three.ss[three.ss$strand == "+", ]$end = three.ss[three.ss$strand == "+", ]$start + 199
three.ss[three.ss$strand == "-", ]$start = three.ss[three.ss$strand == "-", ]$start - 99
three.ss[three.ss$strand == "-", ]$end = three.ss[three.ss$strand == "-", ]$start + 199
three.ss = makeGRangesFromDataFrame(df = three.ss, keep.extra.columns = TRUE)
three.ss_seq = getSeq(x = Hsapiens, three.ss)
##  Loading of motif to be analysed as a standard object of class "matrix".
pwm = data.matrix(frame = data.frame(read.table(file = paste(d.work, f.motif, sep = "/"),
                                                sep = "\t",
                                                header = FALSE,
                                                quote = "\"",
                                                row.names = 1,
                                                as.is = TRUE)),
                  rownames.force = TRUE)
##  Motif analysis for the 5' splice sites of the exon-exon junctions.
five.ss_hits = double()
for (i in 1:length(five.ss_seq)){
     hits = matchPWM(pwm = pwm,
                     subject = five.ss_seq[[i]],
                     min.score = "0%",
                     with.score = TRUE)
     five.ss_hits[[i]] = paste(mcols(hits)$score, collapse = ", ")
}
eej$five.ss_weights = five.ss_hits
##  Motif analysis for the 3' splice sites of the exon-exon junctions.
three.ss_hits = double()
for (j in 1:length(three.ss_seq)){
     hits = matchPWM(pwm = pwm,
                     subject = three.ss_seq[[j]],
                     min.score = "0%",
                     with.score = TRUE)
     three.ss_hits[[j]] = paste(mcols(hits)$score, collapse = ", ")
}
eej$three.ss_weights = three.ss_hits
##  Returning the final object.
return(eej)
}
