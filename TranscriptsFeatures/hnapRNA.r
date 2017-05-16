##############################################################################################################
## High-level R function for calculation of hypothetical "non-alternative"                                  ##
## precursor RNA based on reference annotations of RNAs for gene(-s) of interest.                           ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
##  in.gtf - path to folder and name of the file in GTF/GFF format with reference annotations of RNAs
#   for gene(-s) of interest. It is typically Ensembl annotations or similar one.
##  src - a character vector with name of source of reference annotations. Default value is "Ensembl".
hnapRNA = function(in.gtf, src = "Ensembl"){
## Loading of required auxiliary libraries.
#  This code was successfully tested with library GenomicFeatures v.1.26.3 (and higher)
#  and library rtracklayer v.1.34.2 (and higher).
if (!suppressMessages(require(GenomicFeatures))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("GenomicFeatures")
}
if (!suppressMessages(require(rtracklayer))){
	source("http://bioconductor.org/biocLite.R")
	biocLite("rtracklayer")
}
## Sort of input data.
in.gtf = sort(in.gtf)
## Transforming a GRanges object into GRangesList object.
tr.list = makeGRangesListFromDataFrame(as.data.frame(in.gtf), split.field = "transcript_id")
if (length(tr.list) > 1){
	## Development a list of exons with retained introns.
	ex.intr.retained = lapply(tr.list,
	                          function(x){lapply(tr.list,
	                                      function(y){y[countOverlaps(query = y,
	                                                                  subject = x,
	                                                                  type = "any") > 1, ]})})
	ex.intr.retained = do.call(getMethod(c, "GenomicRanges"), unlist(ex.intr.retained))
	ex.intr.retained = sort(ex.intr.retained[!duplicated(ex.intr.retained), ])
	## Development a list of internal exons.
	ex.internal = do.call(getMethod(c, "GenomicRanges"), lapply(tr.list, function(x){x[c(-1, -length(x)), ]}))
	ex.internal = sort(ex.internal[!duplicated(ex.internal), ])
	## Development a list of start exons.
	if (as.character(strand(in.gtf)@values) == "+"){
	    ex.start = do.call(getMethod(c, "GenomicRanges"),
	                       lapply(tr.list, function(x) {x[1]}))
	    ex.start = ex.start[!ex.start %in% subsetByOverlaps(query = ex.start,
	                                                        subject = ex.internal,
	                                                        type = "equal"), ]
	    ex.first = ex.start[!ex.start %in% subsetByOverlaps(query = ex.start,
	                                                        subject = ex.intr.retained,
	                                                        type = "equal"), ]
	    ex.first = intersect(range(ex.first), ex.first)
	    ex.first = ex.first[1]
	    start(ex.first) = min(start(in.gtf))
	    ex.start = ex.start[!ex.start %in% subsetByOverlaps(query = ex.start,
	                                                        subject = ex.first,
	                                                        type = "any"), ]
	    ex.start = intersect(range(ex.start), ex.start)
	}else{
	    ex.start = do.call(getMethod(c, "GenomicRanges"),
	                       lapply(sort(tr.list, decreasing = TRUE), function(x) {x[1]}))
	    ex.start = ex.start[!ex.start %in% subsetByOverlaps(query = ex.start,
	                                                        subject = ex.internal,
	                                                        type = "equal"), ]
	    ex.first = ex.start[!ex.start %in% subsetByOverlaps(query = ex.start,
	                                                        subject = ex.intr.retained,
	                                                        type = "equal"), ]
	    ex.first = intersect(range(ex.first), ex.first)
	    ex.first = ex.first[length(ex.first)]
	    end(ex.first) = max(end(in.gtf))
	    ex.start = ex.start[!ex.start %in% subsetByOverlaps(query = ex.start,
	                                                        subject = ex.first,
	                                                        type = "any"), ]
	    ex.start = intersect(range(ex.start), ex.start)
	}
	## Development a list of end exons.
	if (as.character(strand(in.gtf)@values) == "+"){
	    ex.end = do.call(getMethod(c, "GenomicRanges"),
	                     lapply(sort(tr.list, decreasing = TRUE), function(x) {x[1]}))
	    ex.end = ex.end[!ex.end %in% subsetByOverlaps(query = ex.end,
	                                                  subject = ex.internal,
	                                                  type = "equal"), ]
	    ex.last = ex.end[!ex.end %in% subsetByOverlaps(query = ex.end,
	                                                   subject = ex.intr.retained,
	                                                   type = "equal"), ]
	    ex.last = intersect(range(ex.last), ex.last)
	    ex.last = ex.last[length(ex.last)]
	    end(ex.last) = max(end(in.gtf))
	    ex.end = ex.end[!ex.end %in% subsetByOverlaps(query = ex.end,
	                                                  subject = ex.last,
	                                                  type = "any"), ]
	    ex.end = intersect(range(ex.end), ex.end)
	}else{
	    ex.end = do.call(getMethod(c, "GenomicRanges"),
	                     lapply(tr.list, function(x) {x[1]}))
	    ex.end = ex.end[!ex.end %in% subsetByOverlaps(query = ex.end,
	                                                  subject = ex.internal,
	                                                  type = "equal"), ]
	    ex.last = ex.end[!ex.end %in% subsetByOverlaps(query = ex.end,
	                                                   subject = ex.intr.retained,
	                                                   type = "equal"), ]
	    ex.last = intersect(range(ex.last), ex.last)
	    ex.last = ex.last[1]
	    start(ex.last) = min(start(in.gtf))
	    ex.end = ex.end[!ex.end %in% subsetByOverlaps(query = ex.end,
	                                                  subject = ex.last,
	                                                  type = "any"), ]
	    ex.end = intersect(range(ex.end), ex.end)
	}
	## Development a list of all exons.
	ex.all = c(ex.first,
	           ex.internal[!ex.internal %in% subsetByOverlaps(query = ex.internal,
	                                                          subject = ex.intr.retained,
	                                                          type = "equal"), ],
	           ex.last)
	ex.all = intersect(range(ex.all), ex.all)
	## Development a list of all introns.
	intr.all = setdiff(range(ex.all), ex.all)
	## Development a list of retained introns.
	intr.retained = intr.all[queryHits(findOverlaps(query = intr.all,
	                                                subject = intersect(range(ex.intr.retained),
	                                                                    ex.intr.retained),
	                                                type = "within")), ]
	## Consolidation of results.
	transcript = GRanges(seqnames = seqnames(ex.all)@values,
	                     ranges = IRanges(start = start(ex.all)[1],
	                                      end = end(ex.all)[length(ex.all)]),
	                     strand = strand(ex.all)@values)
	transcript$type = "transcript"
	transcript$number = NA
	ex.first$type = "first_exon"
	ex.first$number = NA
	if(length(ex.start) > 0){
	   ex.start$type = "alternative_first_exon"
	   ex.start$number = c(1:length(ex.start))
	}
	ex.all$type = "exon"
	ex.all$number = c(1:length(ex.all))
	intr.all$type = "intron"
	intr.all$number = c(1:length(intr.all))
	if(length(ex.intr.retained) > 0){
	   ex.intr.retained$type = "exon_with_retained_intron"
	   ex.intr.retained$number = c(1:length(ex.intr.retained))
	}
	if(length(intr.retained) > 0){
	   intr.retained$type = "retained_intron"
	   intr.retained$number = c(1:length(intr.retained))
	}
	if(length(ex.end) > 0){
	   ex.end$type = "alternative_last_exon"
	   ex.end$number = c(1:length(ex.end))
	}
	ex.last$type = "last_exon"
	ex.last$number = NA
	hnap.RNA = c(transcript,
	             if (as.character(strand(in.gtf)@values) == "+"){
	                 c(ex.first,
	                   sort(c(ex.start, ex.all, intr.all, ex.intr.retained, intr.retained, ex.end)),
	                   ex.last)
	             }else{c(ex.last,
	                     sort(c(ex.start, ex.all, intr.all, ex.intr.retained, intr.retained, ex.end)),
	                     ex.first)})
	hnap.RNA$gene_id = unique(in.gtf$gene_id)
	hnap.RNA$gene_name = paste(unique(in.gtf$gene_name), collapse = ";")
	hnap.RNA$source = as.factor(src)
	hnap.RNA$transcript_id = paste(unique(in.gtf$gene_id), "preRNA", sep = ".")
}else{
	transcript = GRanges(seqnames = seqnames(in.gtf)@values,
	                     ranges = IRanges(start = start(in.gtf)[1],
	                                      end = end(in.gtf)[length(in.gtf)]),
	                     strand = strand(in.gtf)@values)
	transcript$type = "transcript"
	transcript$number = NA
	if (as.character(strand(in.gtf)@values) == "+"){
	    ex.first = in.gtf[1, 0]
	}else{ex.first = in.gtf[length(in.gtf), 0]}
	ex.first$type = "first_exon"
	ex.first$number = NA
	ex.all = in.gtf[, 0]
	ex.all$type = "exon"
	ex.all$number = c(1:length(ex.all))
	intr.all = setdiff(range(ex.all), ex.all)
	intr.all$type = "intron"
	intr.all$number = c(1:length(intr.all))
	if (as.character(strand(in.gtf)@values) == "+"){
	    ex.last = in.gtf[length(in.gtf), 0]
	}else{ex.last = in.gtf[1, 0]}
	ex.last$type = "last_exon"
	ex.last$number = NA
	hnap.RNA = c(transcript,
	             if (as.character(strand(in.gtf)@values) == "+"){
	                 c(ex.first, sort(c(ex.all, intr.all)), ex.last)
	             }else{c(ex.last, sort(c(ex.all, intr.all)), ex.first)})
	hnap.RNA$gene_id = unique(in.gtf$gene_id)
	hnap.RNA$gene_name = paste(unique(in.gtf$gene_name), collapse = ";")
	hnap.RNA$source = as.factor(src)
	hnap.RNA$transcript_id = paste(unique(in.gtf$gene_id), "preRNA", sep = ".")
}
return(hnap.RNA)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/hnaprna.r")
#   in.gtf = "GVV/Ensembl_genes.gtf"
#   src = "Ensembl"
#   res = hnapRNA(in.gtf = in.gtf, src = src)
