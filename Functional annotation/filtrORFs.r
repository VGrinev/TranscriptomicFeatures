###################################################################################################
##  High-level R function for filtration of human Ensembl transcripts with inconsistent ORFs.    ##
##  (c) GNU GPL Vasily V. Grinev, 2020. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work       - character string giving the path to and name of work directory.
##  f.fasta      - character string giving the name of multi-FASTA file.
##  f.mtRNAs     - character string giving the name of TXT file in tab-delimited format
#                  containing the list of mitochondrial transcripts. This file should include
#                  one mandatory field transcript_id.
##  f.orf_coords - character string giving the name of TXT file in tab-delimited format
#                  containing coordinates of ORFs for transcripts to be analyzed. This file should
#                  include three mandatory fields:
#                  i) transcript_id - transcript ID;
#                  ii) start        - start coordinate of ORF;
#                  iii) end         - end coordinate of ORF.
filtrORFs = function(d.work, f.fasta, f.mtRNAs, f.orf_coords){
##  Loading of required auxiliary library.
#   This code was successfully tested with library Biostrings v.2.52.0.
suppressMessages(library(Biostrings))
##  Loading of FASTA file in an XStringSet object.
tr = readDNAStringSet(filepath = paste(d.work, f.fasta, sep = "/"),
                      format = "fasta",
                      use.names = TRUE)
##  Removing of mitochondrial transcripts.
m.tr = read.table(file = paste(d.work, f.mtRNAs, sep = "/"),
                  sep = "\t",
                  header = TRUE,
                  quote = "\"",
                  as.is = TRUE)$transcript_id
tr = tr[!names(tr) %in% m.tr, ]
##  Removing transcripts lacking a canonical start and/or stop codons inside the sequence.
orf_coords = read.table(file = paste(d.work, f.orf_coords, sep = "/"),
                        sep = "\t",
                        header = TRUE,
                        quote = "\"",
                        as.is = TRUE)
orf_coords = orf_coords[!orf_coords$transcript_id %in% m.tr, ]
orf_coords = orf_coords[order(orf_coords$transcript_id), ]
tr = tr[names(tr) %in% orf_coords$transcript_id, ]
tr = tr[order(names(tr)), ]
codon.start = substr(x = tr, start = orf_coords$start, stop = orf_coords$start + 2)
codon.start = attr(x = codon.start[codon.start == "ATG"], which = "names")
codon.stop = substr(x = tr, start = orf_coords$end - 2, stop = orf_coords$end)
codon.stop = attr(x = codon.stop[codon.stop %in% c("TAA", "TAG", "TGA")], which = "names")
tr = tr[names(tr) %in% codon.start[codon.start %in% codon.stop], ]
##  Returning the final object.
return(tr)
}
