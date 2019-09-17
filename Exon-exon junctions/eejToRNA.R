###################################################################################################
##  High-level R function for mapping of exon-exon junctions against sequence of mature RNA.     ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  eej.in - path to folder and name of the TXT file in tab-delimited format with
#   exon-exon junctions. This is input data file which should include four mandatory fields:
#   i) seqnames (name of chromosome or scaffold with prefix "chr");
#   ii) start (start genomic coordinate of the exon-exon junction);
#   iii) end (end genomic coordinate of the exon-exon junction);
#   iv) strand (strand information about exon-exon junction).
##  bs.genome - logical argument. The default value is TRUE, and BSgenome object will be used in
#   that case. Only GRCh38/hg38 assembly of H. sapiens genome is used in current implementation of
#   function.
##  file.fa - path to folder and name of the reference sequence in FASTA/FA format. This argument
#   is valid only with bs.genome = FALSE.
##  RNA - path to folder and name the FASTA/FA file with sequence of mature RNA.
eejToRNA = function(eej.in, bs.genome = TRUE, file.fa, RNA){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries
#   Rsamtools v.1.30.0 and BSgenome.Hsapiens.UCSC.hg38 v.1.4.1.
suppressMessages(require(Rsamtools))
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))
##  Loading of exon-exon junctions.
eej = read.table(file = eej.in, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
eej = eej[, which(colnames(eej)!= "width")]
##  Calculation genomic coordinates
#   of upstream and downstream sequences surrounding exon-exon junctions.
up = t(apply(eej, 1, function(x){if (x[5] == "+"){
                                     x[3] = as.numeric(x[3]) - 29
                                     x[4] = as.numeric(x[3]) + 29
                                }else{
                                     x[3] = as.numeric(x[4])
                                     x[4] = as.numeric(x[3]) + 29
                                }
                                 return(x)
                                }
            )
       )
up = makeGRangesFromDataFrame(data.frame(up), keep.extra.columns = TRUE)
down = t(apply(eej, 1, function(x){if (x[5] == "+"){
                                       x[3] = as.numeric(x[4])
                                       x[4] = as.numeric(x[3]) + 29
                                  }else{
                                       x[3] = as.numeric(x[3]) - 29
                                       x[4] = as.numeric(x[3]) + 29
                                  }
                                   return(x)
                                  }
              )
         )
down = makeGRangesFromDataFrame(data.frame(down), keep.extra.columns = TRUE)
##  Retrieving upstream and downstream sequences surrounding exon-exon junctions.
if (bs.genome == TRUE){
    seq.ref = BSgenome.Hsapiens.UCSC.hg38
}else{
    seq.ref = FaFile(file.fa)
}
up.seq = getSeq(x = seq.ref, up)
names(up.seq) = up$eej_id
up.seq = PDict(up.seq)
down.seq = getSeq(x = seq.ref, down)
names(down.seq) = down$eej_id
down.seq = PDict(down.seq)
##  Loading the sequence of RNA of interest.
seq.target = readDNAStringSet(RNA)[[1]]
##  Matching upstream and downstream sequences against the sequence of RNA of interest.
up.coord = data.frame(unlist(matchPDict(pdict = up.seq, subject = seq.target)))[, c(4, 2)]
down.coord = data.frame(unlist(matchPDict(pdict = down.seq, subject = seq.target)))[, c(4, 1)]
##  Consolidation of results.
up.coord = up.coord[up.coord$names %in% down.coord$names, ]
up.coord = up.coord[order(up.coord$names), ]
down.coord = down.coord[down.coord$names %in% up.coord$names, ]
down.coord = down.coord[order(down.coord$names), ]
coord = cbind(up.coord, down.coord[, 2])
names(coord) = c("eej_id", "start", "end")
return(coord)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/eejToRNA.r")
eej.in = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, all EEJs, filtered.txt"
file.fa = "D:/Vasily Grinev/hg38_der8.fa"
RNA.seq = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, exon 4a started RUNX1-RUNX1T1 mRNA.fasta"
res = eejToRNA(eej.in = eej.in,
               bs.genome = FALSE,
               file.fa = file.fa,
               RNA = RNA.seq)
file.out = "RUNX1-RUNX1T1 project, exon 4a started RUNX1-RUNX1T1 mRNA, coordinates of EEJs.txt"
write.table(res,
            file = file.out,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
