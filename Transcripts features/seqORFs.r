###################################################################################################
##  High-level R function for extraction the sequences of identified ORFs.                       ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments for the function:
##  coord    - path to folder and name of the input TAB-delimited TXT file with coordinates of
#              identified ORFs. It is typically output file of ORFindeR function with three
#              mandatory fields:
#                    i)   transcript_id - IDs of transcripts;
#                    ii)  ORFStart      - start coordinates of identified ORFs;
#                    iii) ORFStart      - stop coordinates of identified ORFs.
##  in.file  - path to folder and name of the input FASTA file with sequence(-s) of interest.
##  out.file - path to folder and name of the output FASTA file to write the sequence(-s) of ORFs.
seqORF = function(coord, in.file, out.file){
##  Loading of required auxiliary library.
#   This code was successfully tested with Biostrings v.2.44.1.
if (!suppressMessages(require(Biostrings))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
}
##  Loading the coordinates of identified ORFs as a standard data.frame object.
coord.ORFs = read.table(file = coord,
                        sep = "\t",
                        header = TRUE,
                        quote = "\"",
                        as.is = TRUE)
coord.ORFs = coord.ORFs[coord.ORFs$ORFStart > 0, ]
coord.ORFs = coord.ORFs[order(coord.ORFs$transcript_id), ]
##  Loading the sequence(-s) of target RNA molecules as a DNAStringSet instance.
seq.set = readDNAStringSet(filepath = in.file, format = "fasta")
seq.set = seq.set[names(seq.set) %in% coord.ORFs$transcript_id, ]
seq.set = seq.set[order(names(seq.set)), ]
coord.ORFs = coord.ORFs[coord.ORFs$transcript_id %in% names(seq.set), ]
##  Extraction of sub-sequences.
sub.seq = subseq(x = rep(x = seq.set[names(seq.set) %in% coord.ORFs$transcript_id, ],
                         times = as.numeric(x = table(coord.ORFs$transcript_id))),
                 start = coord.ORFs$ORFStart,
                 end = coord.ORFs$ORFStop)
##  Write an XStringSet object to a file.
writeXStringSet(x = sub.seq, filepath = out.file)
##  Return the final list of sub-sequences.
return(sub.seq)
}
