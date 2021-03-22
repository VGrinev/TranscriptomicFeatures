###################################################################################################
## High-level R function for calculation of overall statistics                                   ##
## on splicing modes of exon-exon junctions.                                                     ##
## (c) GNU GPL Vasily V. Grinev & Ilia M. Ilyushonak, 2019. grinev_vv[at]bsu.by                  ##
###################################################################################################
### Arguments of function:
##  x - path to folder and name of the TAB file with modes of exon-exon junctions. It is typically
#   modeEEJs.classifier output file. This file should include the columns of all main alternative
#   splicing modes: CI, RI, CE, MCE, MEE, A5SS, A3SS, AFE and ALE.
modeEEJs.statistics = function(x, with.CI = TRUE){
modes.stat = read.table(file = x, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
if (with.CI == TRUE){
modes.stat = modes.stat[, c("CI", "RI", "CE", "MCE", "MEE", "A5SS", "A3SS", "AFE", "ALE")]
modes.stat = data.frame(cbind(c("CI", "RI", "CE", "MCE", "MEE", "A5SS", "A3SS", "AFE", "ALE"),
                              round(colSums(modes.stat)/sum(modes.stat), digits = 3)))
}else{
modes.stat = modes.stat[, c("RI", "CE", "MCE", "MEE", "A5SS", "A3SS", "AFE", "ALE")]
modes.stat = data.frame(cbind(c("RI", "CE", "MCE", "MEE", "A5SS", "A3SS", "AFE", "ALE"),
                              round(colSums(modes.stat)/sum(modes.stat), digits = 3)))
}
colnames(modes.stat) = c("mode", "frequency")
return(modes.stat)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/modeeejs.statistics.r")
modes = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, all EEJs, expression-based modes.txt"
res = modeEEJs.statistics(x = modes, with.CI = FALSE)
##  Saving final results.
file.out = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, all EEJs, expression-based modes, summary.txt"
write.table(res,
            file = file.out,
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE)
