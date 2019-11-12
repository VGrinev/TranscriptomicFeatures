###################################################################################################
##  High-level R function for filtering out of BAM reads with bad CIGAR strings.                 ##
##  (c) GNU GPL                                                                                  ##
##  Vasily V. Grinev, Petr V. Nazarov, 2017-2019. grinev_vv[at]bsu.by, petr.nazarov[at]lih.lu    ##
###################################################################################################
##  This is high-level R function for filtering out of reads with wonky CIGAR strings from BAM
##  files. Current experimental version of function allow to work with only BAM (not SAM) files
##  and only against two type of bad CIGAR: i) CIGAR op has zero length; ii) CIGAR M operator maps
##  off end of reference. In addition, users can specify their own list of unwanted reads.
### Arguments of function:
##  d.work      - character string giving the path to and name of work directory.
##  d.bam       - character string giving the name of BAM directory.
##  f.input     - character string giving the name of input BAM file with reads to be filtered.
##  f.index     - character string giving the name of BAM index file.
##  mappedL     - a logical argument indicating whether to use the mapped length
#                 to calculate the end of read. Default value is FALSE.
##  maxReadL    - integer value indicating the maximum read length in the RNA-Seq library.
#                 Default value is 100.
##  user.qnames - character string giving the name of a TXT file containing user provided list of
#                 qnames with bad CIGAR strings. NULL by default.
##  f.dest      - character string giving the name of BAM file to store filtered reads.
filtrBadCIGAR = function(d.work, d.bam, f.input, f.index,
                         mappedL = FALSE, maxReadL = 100, user.qnames = NULL, f.dest){
##  Loading of required auxiliary library.
#   This code was successfully tested with library Rsamtools v.2.0.0.
suppressMessages(require(Rsamtools))
##  Creation a reference to a BAM file.
bam = BamFile(file = paste(paste(d.work, d.bam, sep = "/"), f.input, sep = "/"),
              index = paste(paste(d.work, d.bam, sep = "/"), f.index, sep = "/"))
##  Extraction of information about reference sequences as a GRanges object.
gr = cbind(seqinfo(bam)@seqnames, seqinfo(bam)@seqlengths)
gr = GRanges(seqnames = gr[, 1], ranges = IRanges(start = 1, end = as.numeric(gr[, 2])))
##  Creation a list of qnames with bad CIGAR strings.
qnames = list()
for (i in 1:length(gr)){
name = scanBam(file = bam,
               param = ScanBamParam(what = c("qname", "pos", "cigar"),
                                    which = gr[i]))[[1]]
m = strsplit(x = sub("[+]$", "", gsub("[0-9]+[A-Z]", "", gsub("M", "+", name$cigar))),
             split = "+",
             fixed = TRUE)
#   Collection of qnames for which CIGAR op has zero length.
qnames[[i]] = c(name$qname[gsub(".+[0-9][a-zA-Z]+0M.+", "0M", name$cigar) == "0M"],
#   Collection of qnames for which CIGAR M operator maps off end of reference.
    if (mappedL == "TRUE"){
                name$qname[name$pos +
                           unlist(x = lapply(X = m,
                                             FUN = function(z){sum(as.numeric(x = z))})) -
                           1 > width(gr[i])]
    }else{
                name$qname[name$pos + maxReadL - 1 > width(gr[i])]
    })
}
qnames = unique(x = unlist(x = qnames[lengths(qnames) > 0]))
if (!is.null(user.qnames)){ 
    u.qnames = read.table(file = user.qnames,
                          sep = "\t",
                          header = TRUE,
                          quote = "\"",
                          as.is = TRUE)$qnames
    qnames = sort(x = unique(x = c(qnames, u.qnames)))
}
##  Filtering out of BAM reads with bad CIGAR strings.
if (!is.null(qnames)){
filtr = FilterRules(list(qName = function(x){!x$qname %in% qnames}))
filterBam(file = paste(paste(d.work, d.bam, sep = "/"), f.input, sep = "/"),
          index = paste(paste(d.work, d.bam, sep = "/"), f.index, sep = "/"),
          destination = paste(paste(d.work, d.bam, sep = "/"), f.dest, sep = "/"),
          indexDestination = FALSE,
          filter = filtr,
          param = ScanBamParam(what = "qname"))
}
return(qnames)
}
