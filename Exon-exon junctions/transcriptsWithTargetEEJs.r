###################################################################################################
##  A high-level R function for the sub-setting of transcripts                                   ##
##  containing exon-exon junctions of interest.                                                  ##
##  (c) GNU GPL Vasily V. Grinev, 2018. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  listEEJs - path to folder and name of the TXT file in tab-delimited format with exon-exon 
#   junctions of interest. This file should include four mandatory fields:
#   i) seqnames (name of chromosome or scaffold with prefix "chr");
#   ii) start (start genomic coordinate of the exon-exon junction);
#   iii) end (end genomic coordinate of the exon-exon junction);
#   iv) strand (strand information about exon-exon junction).
##  listTranscripts - path to folder and name of the GTF/GFF file
#   containing transcripts to be sub-setted.
transcriptsWithTargetEEJs = function(listEEJs, listTranscripts){
##  Loading of required auxiliary libraries.
if (!suppressMessages(require(rtracklayer))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
}
##  Loading of the list of transcripts as a GRangesList object.
trans.all = import(con = listTranscripts)
trans.all = split(x = trans.all, f = trans.all$transcript_id)
##  Development of a complete list of exon-exon junctions of transcripts as GRanges object.
EEJs.all = unlist(setdiff(range(trans.all), trans.all))
start(EEJs.all) = start(EEJs.all) - 1
end(EEJs.all) = end(EEJs.all) + 1
EEJs.all$transcript_id = names(EEJs.all)
names(EEJs.all) = NULL
##  Loading of the list of exon-exon junctions of interest as a GRanges object.
EEJs.target = read.table(listEEJs, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
EEJs.target = makeGRangesFromDataFrame(EEJs.target, keep.extra.columns = TRUE)
##  Sub-setting of the original list of transcripts.
hits  = findOverlaps(query = EEJs.all, subject = EEJs.target, type = "equal")
trans.target_names = sort(unique(EEJs.all[queryHits(hits), ]$transcript_id))
trans.target = unlist(trans.all)
trans.target = trans.target[trans.target$transcript_id %in% trans.target_names, ]
names(trans.target) = NULL
##  Final GRanges object to be returned.
return(trans.target)
}
