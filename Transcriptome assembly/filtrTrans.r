##############################################################
## High-level R function                                    ##
## for filtration of the Cufflinks assembled transcripts.   ##
## (c) GNU GPL Vasily V. Grinev, 2016. grinev_vv[at]bsu.by  ##
##############################################################
## Arguments of function:
# in.gtf - path to folder and name of the input file in GTF/GFF format. It is typically output file of Cufflinks (or similar ones).
# unstr - a logical argument indicating that the unstranded transcripts should be removed. Default value is TRUE.
# difStr - a logical argument specifying that the transcripts that match two different strands should be removed. Default value is TRUE.
# oneEx - a logical argument indicating that the one-exon transcripts should be removed. Default value is TRUE.
# canChr - a logical argument specifying that the records from non-canonical chromosomes should be removed. Default value is TRUE.
# exL - an integer argument. It is a threshold for minimal length of exon(-s) in transcripts. Default Ensembl-based value is 25.
# intrL - an integer argument. It is a threshold for minimal length of intron(-s) in transcripts. Default Ensembl-based value is 50.
# trL - an integer argument. It is a threshold for minimal length of transcripts. Default value is 300.
# fpkm - an integer argument. It is a threshold for minimal transcripts abundance (in FPKM). Default value is 1.
# out.gtf - path to folder and name of the output file.
# format - format of the output file (typically "gtf" or "gff").
# sqlite - a logical argument indicating that the output results should be additionally saved as a local SQLite database. Default value is TRUE.
filtrTrans = function(in.gtf, unstr = TRUE, difStr = TRUE, oneEx = TRUE, canChr = TRUE, exL = 25, intrL = 50, trL = 300, fpkm = 1, out.gtf, format, sqlite = TRUE){
## Required libraries.
	if (!suppressMessages(require(GenomicFeatures))){
		source("http://bioconductor.org/biocLite.R")
		biocLite("GenomicFeatures")
	}
	if (!suppressMessages(require(rtracklayer))){# rtracklayer version 1.32.1 or higher
		source("http://bioconductor.org/biocLite.R")
		biocLite("rtracklayer")
	}
	cat("Required auxiliary libraries have been loaded...\n")
	flush.console()
## Loading of the primary data as a GRanges object.
	tr = import(con = in.gtf)
	cat("Primary data have been loaded...\n")
	flush.console()
## Removing of unstranded transcripts.
	if (unstr == TRUE){
		tr = tr[!strand(tr) == "*", ]
		cat("Unstranded transcripts have been removed...\n")
		flush.console()
	}
## Removing of transcripts that match two different strands.
	if (difStr == TRUE){
		d.str = cbind(as.character(strand(tr)), tr$transcript_id)
		d.str = table(d.str[!duplicated(d.str), ][, 2])
		tr = tr[!tr$transcript_id %in% attr(d.str[d.str > 1], "dimnames")[[1]], ]
		cat("Transcripts that match two different strands have been removed...\n")
		flush.console()
	}
## Removing of one-exon transcripts.
	if (oneEx == TRUE){
		tr = tr[tr$transcript_id %in% attr(table(tr$transcript_id)[table(tr$transcript_id) > length(unique(tr$type))], "dimnames")[[1]], ]
		cat("One-exon transcripts have been removed...\n")
		flush.console()
	}
## Removing of transcripts with too short exon(-s).
	if (exL > 0){
		tr = tr[!tr$transcript_id %in% unique(tr[width(tr) < exL]$transcript_id), ]
		cat("Transcripts with too short exon(-s) have been removed...\n")
		flush.console()
	}
## Only canonical chromosomes.
	if (canChr == TRUE){
		tr = as.data.frame(tr)
		chr = c(paste("chr", c(1:22, "X", "Y"), sep = ""))
		tr = tr[tr$seqnames %in% chr, ][, -8:-9]
		rownames(tr) = NULL
		tr = tr[!duplicated(tr), ]
		tr = makeGRangesFromDataFrame(tr, keep.extra.columns = TRUE)
		cat("Records from non-canonical chromosomes have been removed...\n")
		flush.console()
	}
## Removing of transcripts with too low abundance.
	if (fpkm > 0){
		tr$FPKM = as.numeric(tr$FPKM)
		tr = tr[tr$FPKM >= fpkm, ]
		tr$FPKM = as.character(tr$FPKM)
		cat("Transcripts with too low abundance have been removed...\n")
		flush.console()
	}
## Saving of the GRanges object as GTF/GFF file.
	export(object = tr, con = out.gtf, format = format)
## Loading of the primary data as a TxDb object.
	cat("Creation of the TranscriptDb object...\n")
	flush.console()
	txDb = suppressMessages(makeTxDbFromGFF(file = out.gtf, format = format, circ_seqs = DEFAULT_CIRC_SEQS, chrominfo = NULL, miRBaseBuild = NA))
## Removing too short transcripts.
	if (trL > 0){
		tr = tr[!tr$transcript_id %in% transcriptLengths(txDb)[transcriptLengths(txDb)$tx_len < trL, ]$tx_name, ]
		cat("Too short transcripts have been removed...\n")
		flush.console()
	}
## Removing of transcripts with too short intron(-s).
	if (intrL > 0){
		cat("Removing of transcripts with too short intron(-s) has been started\n",
		"This step may be time-consuming depending from number of processed transcripts, be patient...\n")
		flush.console()
		exons = exonsBy(x = txDb, by = c("tx"), use.names = TRUE)
		intrLength = function(x){
			if (seqnames(x)@lengths > 1){
				if (as.character(strand(x))[1] == "+"){
					intronLength = ranges(x)@start[-1] - (ranges(x)@start + ranges(x)@width)[-length(ranges(x)@start)]
				}
				if (as.character(strand(x))[1] == "-"){
					intronLength = ranges(x)@start[-length(ranges(x)@start)] - (ranges(x)@start + ranges(x)@width)[-1]
				}
				return(min(intronLength))
			}
		}
		shortIntrons = lapply(exons, FUN = intrLength) #This step may be time-consuming depending from number of processed transcripts, be patient.
		shortIntrons = shortIntrons[lengths(shortIntrons) > 0]
		shortIntrons = names(shortIntrons[shortIntrons < intrL])
		tr = tr[!tr$transcript_id %in% shortIntrons, ]
		cat("Transcripts with too short intron(-s) have been removed...\n")
		flush.console()
	}
## Saving of the finnal GRanges object as GTF/GFF file.
	cat("Saving of the finnal GRanges object as GTF/GFF file...\n")
	flush.console()
	export(object = tr, con = out.gtf, format = format)
## Saving of the finnal GRanges object as a local SQLite database.
	if (sqlite == TRUE){
		cat("Saving of the finnal GRanges object as a local SQLite database...\n")
		flush.console()
		txDb = suppressMessages(makeTxDbFromGFF(file = out.gtf, format = format, circ_seqs = DEFAULT_CIRC_SEQS, chrominfo = NULL, miRBaseBuild = NA))
		saveDb(x = txDb, file = sub(format, "sqlite", out.gtf))
	}
}
## A simple example of function use.
# source("/media/gvv/NGS_Illumina_HiSeq/Software/filtrTrans.r")
# input = "/media/gvv/NGS_Illumina_HiSeq/Files_GTF/Set of transcripts.gtf"
# output = "/media/gvv/NGS_Illumina_HiSeq/Set of transcripts, filtered.gtf"
# format = "gtf"
# filtrTrans(in.gtf = input, out.gtf = output, format = format)
