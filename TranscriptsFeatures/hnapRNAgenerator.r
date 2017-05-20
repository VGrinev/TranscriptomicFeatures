##############################################################################################################
## Generic high-level R function for calculation of hypothetical "non-alternative"                          ##
## precursor RNAs based on reference annotations of RNAs for gene(-s) of interest.                          ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
##  hnapRNA.f - path and name of file with main hnapRNA function.
##  gtf.file - path to folder and name of the file in GTF/GFF format with reference annotations of RNAs for
#   gene(-s) of interest. It is typically Ensembl annotations or similar one.
##  output.file - path to folder and name of the TXT output file in tab-delimited format.
##  src - a character vector with name of source of reference annotations. Default value is "Ensembl".
##  cl - an integer argument means number of cores available for parallel processing. Default value is the
#   maximum number of CPU cores on the current host minus 1.
hnapRNAgenerator = function(hnapRNA.f, gtf.file, output.file, src, cl = NULL){
## Connection to hnapRNA function.
source(file = hnapRNA.f)
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
suppressMessages(require(pbapply))
suppressMessages(require(parallel))
## Loading of GTF/GFF input data as a GRanges object.
gtf.file = import(con = gtf.file)
gtf.file = gtf.file[gtf.file$type == "exon", ]
## Filtering of input data against unspliced transcripts.
tr.filter = table(gtf.file$transcript_id) == 1
gtf.file = gtf.file[!gtf.file$transcript_id %in% attr(tr.filter[tr.filter == "TRUE"], "dimnames")[[1]], ]
## Calculation of hnapRNA molecules.
gtf.file = split(gtf.file, list(gtf.file$gene_id))
if (cl == "NULL"){
    cl = detectCores(all.tests = TRUE, logical = TRUE) - 1
}
cl = makeCluster(cl)
clusterExport(cl, c("hnapRNA"))
pboptions(type = "timer", style = 1, char = ".")
res = pblapply(X = gtf.file, FUN = function(x){hnapRNA(x, src = src)}, cl = cl)
stopCluster(cl)
res = do.call(getMethod(c, "GenomicRanges"), res)
## Saving final results.
export(object = res, con = output.file, format = "gtf")
return(res)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/hmrnagenerator.r")
#   main.f = "http://bio.bsu.by/genetics/files/hnaprna.r"
#   in.file = "GVV/Ensembl_genes.gtf"
#   out.file = "GVV/Ensembl_hnapRNA.gtf"
#   res = hmRNAgenerator(hnapRNA.f = main.f,
#                        gtf.file = in.file,
#                        output.file = out.file,
#                        src = "Ensembl",
#                        cl = 2L,
#                        pb = TRUE)
