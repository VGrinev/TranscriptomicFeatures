###################################################################################################
##  A collection of various methods for quality assessment of BAM files.                         ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Assessment of the overall quality of BAM files with Picard v.2.21.2.
##  Test JAVA version.
java -version
##  Test Picard installation.
java -jar D:/Software/Picard/picard.jar -h
##  Validation of BAM file relative to the SAM format specification.
java -jar D:/Software/Picard/picard.jar ValidateSamFile I=D:/Software/QoRTs/BAM_files/KAT8_S3.bam MODE=SUMMARY
### Filtering out of BAM reads with bad CIGAR strings using filtrBadCIGAR function.
source("D:/Software/R codes/filtrBadCIGAR.r")
res = filtrBadCIGAR(d.work = "D:/Software/QoRTs",
                    d.bam = "Files_BAM",
                    f.input = "KAT8_S3.bam",
                    f.index = "KAT8_S3.bam.bai",
                    user.qnames = NULL,
                    f.dest = "KAT8_S3_filtered.bam")
