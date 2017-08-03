##############################################################################################################
## A set of high-level R functions for detection of significant open reading frames in nucleotide sequences ##
## and identification pre-mature translation termination codons.                                            ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
## In part, this code is based on original ideas of Avril Coghlan (Wellcome Trust Sanger Institute,         ##
## Cambridge, U.K.) published 10 June 2011 in "A little book of R for bioinformatics. Release 0.1"          ##
## (https://media.readthedocs.org/pdf/a-little-book-of-r-for-bioinformatics/latest/a-little-book-of-r-for-  ##
## bioinformatics.pdf)                                                                                      ##
##############################################################################################################
### Arguments for the main annotatePTC function:
##  in.gtf - path to folder and name of the input file in GTF/GFF format with transcripts data. It can be
##  Cufflinks output file, Ensembl annotations or similar one.
##  write.transcripts - a logical argument indicating whether to write the retrieved sequences of transcript
##  on the local computer. Default value is FALSE. If TRUE, the file in FASTA format will be saved. The name
##  of this file is based on the name of the input file.
##  threads - an integer argument means number of cores available for parallel processing. Default value is 4.
##  typeORF - the type of the start codon. Default value is "ATG". Noncanonical start codons "GTG", "TTG" and
##  "CTG" are also accepted.
##  n - the number of random sequences to be generated for statistical modeling. Default value is 100.
##  prob - a percentile threshold for identification of the true open reading frame(-s) in the empirical
##  transcript. Default value is 0.99.
##  write.ORFs - a logical argument indicating whether to write the detected significant open reading frames.
##  Default value is FALSE. If TRUE, the TXT file in tab-delimited format will be saved. The name
##  of this file is based on the name of the input file and it contains four fields: 
#   i) transcript (name of transcript);
#   ii) ORFStart (the first coordinate of the identified open reading frame);
#   iii) ORFStop (the last coordinate of the identified open reading frame);
#   iv) ORFLength (the length of the identified open reading frame).
##  write.PTCs - a logical argument indicating whether to write the annotation of transcripts with pre-mature
##  translation termination codons. Default value is FALSE. If TRUE, the TXT file in tab-delimited format will
##  be saved. The name of this file is based on the name of the input file and it contains six fields: 
#   i) name (name of transcript);
#   ii) length (length of transcript);
#   iii) n_exons (number of exons in the transcript);
#   iv) l_exon_width (the width of the last exon in the transcript);
#   v) l_EEJ_pos (the position of the last exon-exon junction in the transcript);
#   vi) status (status of transcript: with mature stop codon, with pre-mature stop codon or non-coding RNA).
### The main annotatePTC function returns an object of class list that includes the input data and all outcome
##  results. Thus, it is not necessary to record all the intermediate results as separate files.
##############################################################################################################
### Auxiliary function for retrieving of the transcript sequences.
#   This function returns a list of DNAString instances.
transcriptSeq = function(y){
    if (as.character(strand(y)@values) == "+"){
        tr.seq = unlist(getSeq(x = BSgenome.Hsapiens.UCSC.hg38, y))
    }else{
        tr.seq = unlist(rev(getSeq(x = BSgenome.Hsapiens.UCSC.hg38, y)))
    }
return(tr.seq)
}
### Auxiliary function to identify all potential start and stop codons in a nucleotide sequence.
#   The function scans a nucleotide sequence in the search of canonical start codon ATG or noncanonical start
#   codons GTG, TTG and CTG as well as canonical stop codons TAA, TAG and TGA. This function returns a list of
#   coordinates in a nucleotide sequence all codons of interest.
startAndStopCodons = function(x){
codons = c("ATG", "GTG", "TTG", "CTG", "TAA", "TAG", "TGA")
for (i in 1:7){
     codon = codons[i]
     occurrences = matchPattern(codon, x)
     codonpositions = occurrences@ranges@start
     numoccurrences = length(codonpositions)
     if (i == 1){
         positions = codonpositions
         types = rep(codon, numoccurrences)
     }
     else{
         positions = append(positions, codonpositions, after = length(positions))
         types = append(types, rep(codon, numoccurrences), after = length(types))
     }
}
indices = order(positions)
positions = positions[indices]
types = types[indices]
mylist = list(positions, types)
return(mylist)
}
### Auxiliary function to identify all possible ATG-start open reading frame(-s) in a nucleotide sequence.
#   This function returns a list of coordinates and lengths of all possible ATG-start open reading frames
#   in a nucleotide sequence of interest.
ORFFindingATG = function(x){
mylist = startAndStopCodons(x)
positions = mylist[[1]]
types = mylist[[2]]
ORFstarts = numeric()
ORFstops = numeric()
ORFlengths = numeric()
numpositions = length(positions)
if (numpositions >= 2){
    for (i in 1:(numpositions - 1)){
         posi = positions[i]
         typei = types[i]
         found = 0
         while (found == 0){
                for (j in (i + 1):numpositions){
                     posj = positions[j]
                     typej = types[j]
                     posdiff = posj - posi
                     posdiffmod3 = posdiff %% 3
                     ORFlength = posj - posi + 3
                     if ((typei == "ATG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
                         numorfs = length(ORFstops)
                         usedstop = -1
                         if (numorfs > 0){
                             for (k in 1:numorfs){
                                  ORFstopk = ORFstops[k]
                                  if (ORFstopk == (posj + 2)){usedstop <- 1}
                             }
                         }
                         if (usedstop == -1){
                             ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
                             ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
                         }
                         found = 1
                         break
                     }
                     if (j == numpositions){found = 1}
                }
         }
    }
}
indices = order(ORFstarts)
ORFstarts = ORFstarts[indices]
ORFstops = ORFstops[indices]
ORFlengths = numeric()
numorfs = length(ORFstarts)
for (i in 1:numorfs){
     ORFstart = ORFstarts[i]
     ORFstop = ORFstops[i]
     ORFlength = ORFstop - ORFstart + 1
     ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
}
mylist = list(ORFstarts, ORFstops, ORFlengths)
return(mylist)
}
### Auxiliary function to identify all possible GTG-start open reading frame(-s) in a nucleotide sequence.
#   This function returns a list of coordinates and lengths of all possible GTG-start open reading frames
#   in a nucleotide sequence of interest.
ORFFindingGTG = function(x){
mylist = startAndStopCodons(x)
positions = mylist[[1]]
types = mylist[[2]]
ORFstarts = numeric()
ORFstops = numeric()
ORFlengths = numeric()
numpositions = length(positions)
if (numpositions >= 2){
    for (i in 1:(numpositions - 1)){
         posi = positions[i]
         typei = types[i]
         found = 0
         while (found == 0){
                for (j in (i + 1):numpositions){
                     posj = positions[j]
                     typej = types[j]
                     posdiff = posj - posi
                     posdiffmod3 = posdiff %% 3
                     ORFlength = posj - posi + 3
                     if ((typei == "GTG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
                         numorfs = length(ORFstops)
                         usedstop = -1
                         if (numorfs > 0){
                             for (k in 1:numorfs){
                                  ORFstopk = ORFstops[k]
                                  if (ORFstopk == (posj + 2)){usedstop <- 1}
                             }
                         }
                         if (usedstop == -1){
                             ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
                             ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
                         }
                         found = 1
                         break
                     }
                     if (j == numpositions){found = 1}
                }
         }
    }
}
indices = order(ORFstarts)
ORFstarts = ORFstarts[indices]
ORFstops = ORFstops[indices]
ORFlengths = numeric()
numorfs = length(ORFstarts)
for (i in 1:numorfs){
     ORFstart = ORFstarts[i]
     ORFstop = ORFstops[i]
     ORFlength = ORFstop - ORFstart + 1
     ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
}
mylist = list(ORFstarts, ORFstops, ORFlengths)
return(mylist)
}
### Auxiliary function to identify all possible TTG-start open reading frame(-s) in a nucleotide sequence.
#   This function returns a list of coordinates and lengths of all possible TTG-start open reading frames
#   in a nucleotide sequence of interest.
ORFFindingTTG = function(x){
mylist = startAndStopCodons(x)
positions = mylist[[1]]
types = mylist[[2]]
ORFstarts = numeric()
ORFstops = numeric()
ORFlengths = numeric()
numpositions = length(positions)
if (numpositions >= 2){
    for (i in 1:(numpositions - 1)){
         posi = positions[i]
         typei = types[i]
         found = 0
         while (found == 0){
                for (j in (i + 1):numpositions){
                     posj = positions[j]
                     typej = types[j]
                     posdiff = posj - posi
                     posdiffmod3 = posdiff %% 3
                     ORFlength = posj - posi + 3
                     if ((typei == "TTG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
                         numorfs = length(ORFstops)
                         usedstop = -1
                         if (numorfs > 0){
                             for (k in 1:numorfs){
                                  ORFstopk = ORFstops[k]
                                  if (ORFstopk == (posj + 2)){usedstop <- 1}
                             }
                         }
                         if (usedstop == -1){
                             ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
                             ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
                         }
                         found = 1
                         break
                     }
                     if (j == numpositions){found = 1}
                }
         }
    }
}
indices = order(ORFstarts)
ORFstarts = ORFstarts[indices]
ORFstops = ORFstops[indices]
ORFlengths = numeric()
numorfs = length(ORFstarts)
for (i in 1:numorfs){
     ORFstart = ORFstarts[i]
     ORFstop = ORFstops[i]
     ORFlength = ORFstop - ORFstart + 1
     ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
}
mylist = list(ORFstarts, ORFstops, ORFlengths)
return(mylist)
}
### Auxiliary function to identify all possible CTG-start open reading frame(-s) in a nucleotide sequence.
#   This function returns a list of coordinates and lengths of all possible CTG-start open reading frames
#   in a nucleotide sequence of interest.
ORFFindingCTG = function(x){
mylist = startAndStopCodons(x)
positions = mylist[[1]]
types = mylist[[2]]
ORFstarts = numeric()
ORFstops = numeric()
ORFlengths = numeric()
numpositions = length(positions)
if (numpositions >= 2){
    for (i in 1:(numpositions - 1)){
         posi = positions[i]
         typei = types[i]
         found = 0
         while (found == 0){
                for (j in (i + 1):numpositions){
                     posj = positions[j]
                     typej = types[j]
                     posdiff = posj - posi
                     posdiffmod3 = posdiff %% 3
                     ORFlength = posj - posi + 3
                     if ((typei == "CTG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
                         numorfs = length(ORFstops)
                         usedstop = -1
                         if (numorfs > 0){
                             for (k in 1:numorfs){
                                  ORFstopk = ORFstops[k]
                                  if (ORFstopk == (posj + 2)){usedstop <- 1}
                             }
                         }
                         if (usedstop == -1){
                             ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
                             ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
                         }
                         found = 1
                         break
                     }
                     if (j == numpositions){found = 1}
                }
         }
    }
}
indices = order(ORFstarts)
ORFstarts = ORFstarts[indices]
ORFstops = ORFstops[indices]
ORFlengths = numeric()
numorfs = length(ORFstarts)
for (i in 1:numorfs){
     ORFstart = ORFstarts[i]
     ORFstop = ORFstops[i]
     ORFlength = ORFstop - ORFstart + 1
     ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
}
mylist = list(ORFstarts, ORFstops, ORFlengths)
return(mylist)
}
### Auxiliary function for generation of a set of random nucleotide sequences with the same length as original
#   nucleotide sequence using a multinomial model. This function returns a vector of new sequences.
generateSeqsWithMultinomialModel = function(x, n){
inputsequencevector = s2c(x)
mylength = length(inputsequencevector)
mytable = table(inputsequencevector)
letters = rownames(mytable)
numletters = length(letters)
probabilities = numeric()
for (i in 1:numletters){
     letter = letters[i]
     count = mytable[[i]]
     probabilities[i] = count/mylength
}
seqs = numeric(n)
for (j in 1:n){
     seq = sample(letters, mylength, rep = TRUE, prob = probabilities)
     seq = c2s(seq)
     seqs[j] = seq
}
return(seqs)
}
### An auxiliary function to identify a significant open reading frame(-s) in a sequence of interest.
#   This function returns a data frame with coorditanes of identified open reading frame(-s).
ORFindeR = function(x, typeORF, n, prob){
if (typeORF == "ATG"){
    ORF = ORFFindingATG(as.character(x))
}
if (typeORF == "CTG"){
    ORF = ORFFindingCTG(as.character(x))
}
if (typeORF == "GTG"){
    ORF = ORFFindingGTG(as.character(x))
}
if (typeORF == "TTG"){
    ORF = ORFFindingTTG(as.character(x))
}
if (length(ORF[[1]]) == 0){
    ORF = data.frame(cbind("transcript", 0, 0, 0), stringsAsFactors = FALSE)
    colnames(ORF) = c("transcript", "ORFStart", "ORFStop", "ORFLength")
}
else if (length(ORF[[1]]) >= 1){
    ORF = data.frame(cbind(ORF[[1]], ORF[[2]], ORF[[3]]))
    RANDSeq = generateSeqsWithMultinomialModel(as.character(x), n)
    RANDSeqORFlengths = numeric()
    for (k in 1:length(RANDSeq)){
         RANDSeq1 = RANDSeq[k]
         RANDSeq1 = ORFFindingATG(RANDSeq1)
         lengths = RANDSeq1[[3]]
         RANDSeqORFlengths = append(RANDSeqORFlengths, lengths, after = length(RANDSeqORFlengths))
    }
    RAND = subset(RANDSeqORFlengths, RANDSeqORFlengths > 0)
    RAND = quantile(RAND, probs = prob)
    ORF = ORF[!ORF[, 3] < RAND, ]
    if (length(ORF$X1) == 0){
        ORF[1, ] = 0
        ORF = cbind("transcript", ORF)
        colnames(ORF) = c("transcript", "ORFStart", "ORFStop", "ORFLength")
    }
    if (length(ORF$X1) >= 1){
        ORF = cbind("transcript", ORF)
        colnames(ORF) = c("transcript", "ORFStart", "ORFStop", "ORFLength")
    }
}
return(ORF)
}
### Main function to identify transcripts with PTC.
annotatePTC = function(in.gtf, write.transcripts = FALSE, threads = 4, typeORF = "ATG", n = 100, prob = 0.99, write.ORFs = FALSE, write.PTCs = FALSE){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with following libraries: Biostrings v.2.44.1, BSgenome v.1.44,
#   BSgenome.Hsapiens.UCSC.hg38 v.1.4.1, rtracklayer v.1.36.3 and seqinr v.3.3-6.
cat("Loading of required auxiliary libraries.\n")
flush.console()
if (!suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg38")
}
if (!suppressMessages(require(seqinr))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("seqinr")
}
suppressMessages(require(tools))
suppressMessages(require(pbapply))
##  Loading of GTF/GFF input data as a GRanges object.
cat("Loading of GTF/GFF input data as a GRangesList object.\n")
flush.console()
input.gtf = import(con = in.gtf)
##  Transforming a GRanges object into GRangesList object.
tr.list = sort(makeGRangesListFromDataFrame(as.data.frame(input.gtf), split.field = "transcript_id"))
##  Retrieving of the transcript sequences as a list of DNAString instances.
cat("Retrieving of the transcript sequences as a list of DNAString instances.\n")
flush.console()
tr.seq.list = lapply(tr.list, function(y){transcriptSeq(y)})
##  Setting of multithreaded computing.
cluster = makeCluster(threads)
clusterExport(cluster, c("ORFindeR",
                         "ORFFindingATG",
                         "startAndStopCodons",
                         "matchPattern",
                         "generateSeqsWithMultinomialModel",
                         "s2c", "c2s"))
##  Identification of significant open reading frames in the nucleotide sequences of interest.
cat("Identification of significant open reading frames in the nucleotide sequences of interest.\n",
"This step may be time-consuming depending from number of processed transcripts, be patient.\n")
flush.console()
pboptions(type = "timer", style = 1, char = ".")
ORFs = pblapply(tr.seq.list, function(x){ORFindeR(x, typeORF = "ATG", n = 100, prob = 0.98)}, cl = cluster)
#ORFs = pblapply(tr.seq.list, function(x){ORFindeR(x, typeORF = typeORF, n = n, prob = prob)}, cl = cluster)
stopCluster(cluster)
ORFs = data.frame(do.call(rbind, ORFs), stringsAsFactors = FALSE)
ORFs$transcript = rownames(ORFs)
rownames(ORFs) = NULL
ORFs$transcript = gsub("\\..+", "", ORFs$transcript)
ORFs[, c(2:4)] = apply(ORFs[, c(2:4)], 2, function(x){as.numeric(x)})
cat("The identification of significant open read frames has been successfully completed.\n")
flush.console()
##  Annotation of transcripts with premature translation-termination codons.
cat("Annotation of transcripts with PTCs.\n")
flush.console()
tr.annotated = data.frame(cbind(names(tr.list),
                                sum(width(tr.list)),
                                lengths(width(tr.list))),
                          stringsAsFactors = FALSE)
last.exon = unlist(lapply(tr.list,
                          function(x){if (as.character(strand(x)@values) == "+"){
                                          last.exon = width(rev(x)[1, ])
                                      }else{
                                          last.exon = width(x[1, ])
                                      }
                                     }
                         )
                  )
tr.annotated$X4 = last.exon
tr.annotated[, c(2:4)] = apply(tr.annotated[, c(2:4)], 2, function(x){as.numeric(as.character(x))})
tr.annotated$X5 = tr.annotated$X2 - tr.annotated$X4 + 1
colnames(tr.annotated) = c("name", "length", "n_exons", "l_exon_width", "l_EEJ_pos")
rownames(tr.annotated) = NULL
stop_pos = merge(aggregate(ORFLength ~ transcript, ORFs, max), ORFs, by = c("transcript", "ORFLength"))
if (max(table(stop_pos$transcript)) > 1){
    stop_pos = merge(aggregate(ORFStart ~ transcript, ORFs, min), ORFs, by = c("transcript", "ORFStart"))
}
tr.annotated$stop_pos = stop_pos$ORFStop
tr.annotated$status = ""
if (length(tr.annotated[(tr.annotated$stop_pos - tr.annotated$l_EEJ_pos - 1) > -50, ]$status) > 0){
    tr.annotated[(tr.annotated$stop_pos - tr.annotated$l_EEJ_pos - 1) > -50, ]$status = "mature_stop"
}
if (length(tr.annotated[(tr.annotated$stop_pos - tr.annotated$l_EEJ_pos - 1) < -49, ]$status) > 0){
    tr.annotated[(tr.annotated$stop_pos - tr.annotated$l_EEJ_pos - 1) < -49, ]$status = "premature_stop"
}
if (length(tr.annotated[tr.annotated$stop_pos == 0, ]$status) > 0){
    tr.annotated[tr.annotated$stop_pos == 0, ]$status = "non_coding_RNA"
}
cat("Annotation of transcripts with PTCs has been successfully completed.\n")
flush.console()
##  Consolidation and saving results.
cat("Consolidation and saving results.\n")
flush.console()
if (write.transcripts == "TRUE"){
    writeXStringSet(DNAStringSet(tr.seq.list),
                    file = sub(file_ext(in.gtf), "sequences.fasta", in.gtf),
                    format = "fasta",
                    width = 100)
}
if (write.ORFs == "TRUE"){
    write.table(ORFs,
                file = sub(file_ext(in.gtf), "ORFs.txt", in.gtf),
                sep = "\t",
                quote = FALSE,
                col.names = TRUE,
                row.names = FALSE)
}
if (write.PTCs == "TRUE"){
    write.table(tr.annotated,
                file = sub(file_ext(in.gtf), "PTCs.txt", in.gtf),
                sep = "\t",
                quote = FALSE,
                col.names = TRUE,
                row.names = FALSE)
}
annotations = list(transcripts.structure = input.gtf,
                   transcripts.sequences = tr.seq.list,
                   transcripts.ORFs = ORFs,
                   transcripts.PTCs = tr.annotated)
cat("Annotation has been successfully finished.\n")
flush.console()
return(annotations)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/annotateptc.r")
#   in.file = "GVV/Cufflinks_transcripts.gtf"
#   res = annotatePTC(in.gtf = in.file,
#                     write.transcripts = TRUE,
#                     threads = 28,
#                     typeORF = "ATG",
#                     n = 100,
#                     prob = 0.99,
#                     write.ORFs = TRUE,
#                     write.PTCs = TRUE)
