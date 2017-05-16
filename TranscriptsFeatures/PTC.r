##############################################################
## Identification                                           ##
## of premature termination codons in transcripts           ##
## (c) GNU GPL Vasily V. Grinev, 2016. grinev_vv[at]bsu.by  ##
##############################################################
rm(list = ls())
## Required libraries.
if (!require(GenomicFeatures)){
	source("http://bioconductor.org/biocLite.R")
	biocLite("GenomicFeatures")
}
if (!require(SGSeq)){
	source("http://bioconductor.org/biocLite.R")
	biocLite("SGSeq")
}
## Directory(-ies) and file(-s).
workDir = "/media/gvv/NGS_Illumina_HiSeq" #Path to work directory.
gtfDir = "/media/gvv/NGS_Illumina_HiSeq/Files_GTF" #Path to directory with GTF/GFF file(-s).
gtf = "Set of transcripts. Annotations.gtf" #GTF/GFF file(-s) with annotations.
orf = "Set of transcripts. Coordinates of ORFs.txt" #File with coordinates of open reading frames produced by function ORFindeR.
## (If required) Converting of the GTF/GFF file into local SQLite database.
txDb = makeTxDbFromGFF(file = paste(gtfDir, gtf, sep = "/"),
					   format = "auto", dataSource = "Kasumi-1 cells", organism = "Homo sapiens",
					   circ_seqs = DEFAULT_CIRC_SEQS, chrominfo = NULL, miRBaseBuild = NA)
saveDb(x = txDb, file = paste(gtfDir, sub("gtf", "sqlite", gtf, ignore.case = TRUE), sep = "/")) #Be careful and indicate correct GTF/GFF file extension in argument "pattern" of sub function.
## (If above step was ignored) Loading of the transcriptional models of genes as a TranscriptDb object. 
txDb = loadDb(paste(gtfDir, sub("gtf", "sqlite", gtf, ignore.case = TRUE), sep = "/"))
## Creation of the ordered set of annotated transcripts.
# In this set, each transcript will be annotated with name, length, number of exons, width of last exon and position of last exon-exon junction.
tr = data.frame(cbind(transcriptLengths(txDb)$tx_name, transcriptLengths(txDb)$tx_len, transcriptLengths(txDb)$nexon), stringsAsFactors = FALSE)
tr = tr[order(tr$X1), ]
txFe = convertToTxFeatures(txDb)
txFe = txFe[type(txFe)== "L", ] #Retrieving of the last exons of transcripts as a TxFeatures object.
lastExon = data.frame(do.call(rbind, lapply(txFe, function(x){return(cbind(txName(x)[[1]], width(x)))})), stringsAsFactors = FALSE)
lastExon = lastExon[order(lastExon$X1), ] #Retrieving of the width of last exons as a data frame object.
tr$X4 = lastExon$X2
tr[, c(2:4)] = apply(tr[, c(2:4)], 2, function(x) as.numeric(as.character(x)))
tr$X5 = tr$X2 - tr$X4 + 1
colnames(tr) = c("name", "length", "n_exons", "l_exon_width", "l_EEJ_pos")
rownames(tr) = NULL
## Loading of the coordinates of open reading frames as a data frame object.
ORFs = read.table(paste(workDir, orf, sep = "/"), sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
# Some transcripts may contain 2 or more open reading frames. The code below takes into account only open reading frame with the earliest termination.
ORFs = aggregate(ORFs$ORFEnd, by = list(ORFs$Transcript), FUN = min)
## Annotation of transcripts with premature termination codons.
tr$stop_pos = ORFs$x
tr$status = ""
tr[(tr$stop_pos - tr$l_EEJ_pos - 1) < -50, ]$status = "premature_stop"
tr[(tr$stop_pos - tr$l_EEJ_pos - 1) > -50, ]$status = "mature_stop"
tr[tr$stop_pos == 0, ]$status = "non_coding_RNA"
## Saving results in tab-delimited TXT file.
write.table(tr, file = paste(workDir, sub("Coordinates of ORFs", "PTC status", orf, ignore.case = TRUE), sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
