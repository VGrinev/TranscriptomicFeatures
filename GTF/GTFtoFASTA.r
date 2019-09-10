###################################################################################################
##  High-level R function for converting of the GTF/GFF annotations                              ##
##  into FASTA sequences of genomic features.                                                    ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work - path to and name of the work folder;
##  in.gtf - name of the GTF/GFF file with annotations;
##  out.fasta - name of the output file in FSTA format.
GTFtoFASTA = function(d.work, in.gtf, out.fasta){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries rtracklayer v.1.44.4, GenomicRanges v.1.36.1,
#   BSgenome v.1.52.0, BSgenome.Hsapiens.UCSC.hg38 v.1.4.1 and Biostrings v.2.52.0.
suppressMessages(library(rtracklayer))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
##  Loading of the the GTF/GFF annotations as a GRanges object. 
gtf = sort(x = import(con = paste(d.work, in.gtf, sep = "/"), format = "gtf"))
##  Converting a GRanges object into transcript-splitted GRangesList object.
list.gtf = split(x = gtf, f = gtf$transcript_id)
##  Retrieving of transcript sequences as a list of DNAString instances. This step has not been
#   optimized and may take a long time. Be patient!
list.seq = lapply(X = list.gtf, FUN = function(x){if (as.character(strand(x)@values) == "-"){
                                                      unlist(x = rev(x = getSeq(Hsapiens, x)))
                                                 }else{unlist(x = getSeq(Hsapiens, x))}})
##  Consolidation of the DNAString objects into a DNAStringSet instance.
list.seq = DNAStringSet(x = list.seq)
##  Write an XStringSet object to a multiFASTA file.
writeXStringSet(x = list.seq,
                filepath = paste(d.work, out.fasta, sep = "/"),
                format = "fasta",
                width = 100)
}
GTFtoFASTA(d.work = "D:/Vasily Grinev",
           in.gtf = "Ensembl release 85, GRCh38.p7, mRNAs.gtf",
           out.fasta = "Ensembl release 85, GRCh38.p7, mRNAs1.fasta")
