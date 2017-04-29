##############################################################################################################
## High-level R function for reference-based functional annotation of experimentally detected exons.        ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
## models - path to folder and name of the SQLite database with reference transcriptional models of genes.
## ctd.chrs - a logical argument specifying should only standard chromosomes (no alternative assemblies)
#  be used in subsequent analysis. Default value is TRUE.
## experimental - path to folder and name of the TXT file in tab-delimited format with experimentally
#  detected exons. This is input data file which should include four mandatory fields:
#  i) seqnames (name of chromosome or scaffold with prefix "chr");
#  ii) start (start genomic coordinate of the exon);
#  iii) end (end genomic coordinate of the exon);
#  iv) strand (strand information about exon).
## overlap.type - type of overlap between experimental and reference exons. Default value is "within".
#  Alternativelly, value "any", "start", "end" or "equal" can be also used.
exonTypes = function(models, ctd.chrs = "TRUE", experimental, overlap.type = "within"){
	## Loading of required auxiliary library.
	#  This code was tested only with v.1.26.3 of GenomicFeatures library.
	cat("Loading of required auxiliary library.\n")
	flush.console()
	if (!suppressMessages(require(GenomicFeatures))){
	    source("http://bioconductor.org/biocLite.R")
	    biocLite("GenomicFeatures")
	}
	## Loading of reference transcriptional models of genes.
	#  This code was tested only with Ensembl models of human genes (Ensembl release 85, GRCh38.p7 human
	#  reference genome, July 2016).
	cat("Loading of reference transcriptional models of genes.\n")
	flush.console()
	tx = loadDb(models)
	## Development a list of reference exons as an object of class GRangesList.
	cat("Development a list of reference exons.\n")
	flush.console()
	#  Retrieving a full set of exons from reference models.
	ex = exons(tx, columns = NULL)
	ex = ex[!duplicated(ex), ]
	if ((gsub("chr.+", "chr", as.character(seqnames(ex)@values))[1] == "chr") == "FALSE"){
	    ex = as.data.frame(ex)
	    ex$seqnames = paste("chr", as.character(ex$seqnames), sep = "")
	    ex = makeGRangesFromDataFrame(ex, keep.extra.columns = TRUE)
	}
	if (ctd.chrs == "TRUE"){
	    ex = as.data.frame(ex)
	    ex = ex[ex$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
	    ex = makeGRangesFromDataFrame(ex, keep.extra.columns = TRUE)
	}
	#  Retrieving a set of 5'UTR exons from reference models.
	f.ex = unlist(fiveUTRsByTranscript(tx))[, -1:-3]
	names(f.ex) = NULL
	f.ex = f.ex[!duplicated(f.ex), ]
	if ((gsub("chr.+", "chr", as.character(seqnames(f.ex)@values))[1] == "chr") == "FALSE"){
	    f.ex = as.data.frame(f.ex)
	    f.ex$seqnames = paste("chr", as.character(f.ex$seqnames), sep = "")
	    f.ex = makeGRangesFromDataFrame(f.ex, keep.extra.columns = TRUE)
	}
	if (ctd.chrs == "TRUE"){
	    f.ex = as.data.frame(f.ex)
	    f.ex = f.ex[f.ex$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
	    f.ex = makeGRangesFromDataFrame(f.ex, keep.extra.columns = TRUE)
	}
	f.ex = ex[ex %in% subsetByOverlaps(query = ex, subject = f.ex, type = "equal"), ]
	#  Retrieving a set of 3'UTR exons from reference models.
	cds.ex = unlist(cds(tx))[, -1:-3]
	names(cds.ex) = NULL
	cds.ex = cds.ex[!duplicated(cds.ex), ]
	if ((gsub("chr.+", "chr", as.character(seqnames(cds.ex)@values))[1] == "chr") == "FALSE"){
	    cds.ex = as.data.frame(cds.ex)
	    cds.ex$seqnames = paste("chr", as.character(cds.ex$seqnames), sep = "")
	    cds.ex = makeGRangesFromDataFrame(cds.ex, keep.extra.columns = TRUE)
	}
	if (ctd.chrs == "TRUE"){
	    cds.ex = as.data.frame(cds.ex)
	    cds.ex = cds.ex[cds.ex$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
	    cds.ex = makeGRangesFromDataFrame(cds.ex, keep.extra.columns = TRUE)
	}
	cds.ex = ex[ex %in% subsetByOverlaps(query = ex, subject = cds.ex, type = "equal"), ]
	#  Retrieving a set of 3'UTR exons from reference models.
	t.ex = unlist(threeUTRsByTranscript(tx))[, -1:-3]
	names(t.ex) = NULL
	t.ex = t.ex[!duplicated(t.ex), ]
	if ((gsub("chr.+", "chr", as.character(seqnames(t.ex)@values))[1] == "chr") == "FALSE"){
	    t.ex = as.data.frame(t.ex)
	    t.ex$seqnames = paste("chr", as.character(t.ex$seqnames), sep = "")
	    t.ex = makeGRangesFromDataFrame(t.ex, keep.extra.columns = TRUE)
	}
	if (ctd.chrs == "TRUE"){
	    t.ex = as.data.frame(t.ex)
	    t.ex = t.ex[t.ex$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
	    t.ex = makeGRangesFromDataFrame(t.ex, keep.extra.columns = TRUE)
	}
	t.ex = ex[ex %in% subsetByOverlaps(query = ex, subject = t.ex, type = "equal"), ]
	#  Retrieving a set of reference exons from non-coding RNAs.
	nc.tr = transcriptLengths(tx, with.cds_len = TRUE)
	nc.tr = nc.tr[nc.tr$cds_len == 0, ]$tx_name
	nc.ex = unlist(exonsBy(tx, by = "tx", use.names = TRUE))[, -1:-3]
	nc.ex = nc.ex[names(nc.ex) %in% nc.tr, ]
	names(nc.ex) = NULL
	nc.ex = nc.ex[!duplicated(nc.ex), ]
	if ((gsub("chr.+", "chr", as.character(seqnames(nc.ex)@values))[1] == "chr") == "FALSE"){
	    nc.ex = as.data.frame(nc.ex)
	    nc.ex$seqnames = paste("chr", as.character(nc.ex$seqnames), sep = "")
	    nc.ex = makeGRangesFromDataFrame(nc.ex, keep.extra.columns = TRUE)
	}
	if (ctd.chrs == "TRUE"){
	    nc.ex = as.data.frame(nc.ex)
	    nc.ex = nc.ex[nc.ex$seqnames %in% paste("chr", c(1:22, "X", "Y"), sep = ""), ]
	    nc.ex = makeGRangesFromDataFrame(nc.ex, keep.extra.columns = TRUE)
	}
	nc.ex = ex[ex %in% subsetByOverlaps(query = ex, subject = nc.ex, type = "equal"), ]
	#  Development a final list of reference exons.
	f.ex = f.ex[!f.ex %in% subsetByOverlaps(query = f.ex,
	                                        subject = c(cds.ex, t.ex, nc.ex),
	                                        type = "equal"), ]
	cds.ex = cds.ex[!cds.ex %in% subsetByOverlaps(query = cds.ex,
	                                              subject = c(f.ex, t.ex, nc.ex),
	                                              type = "equal"), ]
	t.ex = t.ex[!t.ex %in% subsetByOverlaps(query = t.ex,
	                                        subject = c(f.ex, cds.ex, nc.ex),
	                                        type = "equal"), ]
	nc.ex = nc.ex[!nc.ex %in% subsetByOverlaps(query = nc.ex,
	                                           subject = c(f.ex, cds.ex, t.ex),
	                                           type = "equal"), ]
	m.ex = ex[!ex %in% subsetByOverlaps(query = ex,
	                                    subject = c(f.ex, cds.ex, t.ex, nc.ex),
	                                    type = "equal"), ]
	ref.ex = GRangesList("fiveUTRs" = f.ex,
	                     "CDS" = cds.ex,
	                     "threeUTRs" = t.ex,
	                     "non.coding" = nc.ex,
	                     "multi.type" = m.ex)
	## Loading of experimentally detected exons.
	#  This step produces a standard object of class GRanges.
	cat("Loading of experimentally detected exons.\n")
	flush.console()
	exp.ex = read.table(file = experimental, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
	exp.ex = makeGRangesFromDataFrame(exp.ex, keep.extra.columns = TRUE)
	## Detection of overlaps between experimental and reference exons.
	cat("Detection of overlaps between experimental and reference exons.\n")
	flush.console()
	exp.ex$fiveUTRs = countOverlaps(query = exp.ex, subject = ref.ex[[1]], type = overlap.type)
	exp.ex$CDS = countOverlaps(query = exp.ex, subject = ref.ex[[2]], type = overlap.type)
	exp.ex$threeUTRs = countOverlaps(query = exp.ex, subject = ref.ex[[3]], type = overlap.type)
	exp.ex$non.coding = countOverlaps(query = exp.ex, subject = ref.ex[[4]], type = overlap.type)
	exp.ex$multi.type = countOverlaps(query = exp.ex, subject = ref.ex[[5]], type = overlap.type)
	cat("Functional annotation of experimental exons has been successfully finished.\n")
	flush.console()
	#  This function returns an object of class GRanges which contain all input data and includes five new
	#  metadata columns with overlap counts.
	return(exp.ex)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/exontypes.r")
#   modelData = "GVV/Ensembl_GRCh38.p7_release.85.sqlite"
#   exprData = "GVV/Experimental_exons.txt"
#   res = exonTypes(models = modelData, experimental = exprData)
