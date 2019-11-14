###################################################################################################
##  Consolidated high-level R function for local and global alignment of short RNA-Seq reads     ##
##  by seed-and-vote algorithm.                                                                  ##
##  (c) GNU GPL Vasily V. Grinev, 2018-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work - path to and the name of work directory.
##  ref.Genome - a character vector giving the path to and basename of index file.
##  d.fastq - the name of directory with FASTQ files.
##  file.fastq1 - a character vector including name of FASTQ file that contain reads to be aligned.
#   For paired-end sequencing, it is file including first reads in RNA-Seq library.
##  file.fastq2 - a character vector giving name of file that include second reads in paired-end
#   sequencing data. Default value is NULL
##  d.bam - the name of output directory for BAM file(s).
##  thr - integer giving the number of threads used for running of subjunc function.
#   Default value is 1.
### See the Rsubread package documentation for further details
#   on arguments of low-level function subjunc.
alignSubjunc = function(d.work, ref.Genome, d.fastq, file.fastq1, file.fastq2 = NULL, d.bam, thr){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries Rsubread v.1.34.6 and Rsamtools v.1.28.0.
suppressMessages(library(package = Rsubread))
suppressMessages(library(package = Rsamtools))
##  Assignment the basename of index file.
ref.Genome = paste(d.work, ref.Genome, sep = "/")
##  Assignment the path to and the name of FASTQ file including first reads in RNA-Seq library.
file.fastq1 = paste(paste(d.work, d.fastq, sep = "/"), file.fastq1, sep = "/")
##  Assignment the path to and the name of FASTQ file including second reads in RNA-Seq library.
if (is.null(file.fastq2) == FALSE){
file.fastq2 = paste(paste(d.work, d.fastq, sep = "/"), file.fastq2, sep = "/")
}
##  Assignment the name of output BAM file.
file.bam = sub("fast.+", "bam", sub(paste(d.work, d.fastq, sep = "/"), "", file.fastq1))
##  Alignment of short RNA-Seq reads against reference genome.
align = subjunc(index = ref.Genome,
                readfile1 = file.fastq1,
                readfile2 = file.fastq2,
                input_format = "gzFASTQ",
                output_format = "BAM",
                output_file = paste(paste(d.work, d.bam, sep = "/"), file.bam, sep = ""),
                unique = TRUE,
                nBestLocations = 1,
                nthreads = thr,
                reportAllJunctions = TRUE)
##  Sorting and indexing of generated BAM file.
s.bam = sortBam(file = paste(paste(d.work, d.bam, sep = "/"), file.bam, sep = ""),
                destination = sub(".bam",
                                  "",
                                  paste(paste(d.work, d.bam, sep = "/"), file.bam, sep = "")),
                byQname = FALSE)
suppressMessages(indexBam(s.bam))
}
