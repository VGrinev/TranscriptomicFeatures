##############################################################################################################
## High-level R function for consolidation of global statistics on splicing distances.                      ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by                                                  ##
##############################################################################################################
### Arguments of function:
##  input - path to folder and name of the input txt-file in tab-delimited format.
##  log.FC - an integer argument. It is a threshold for log2 fold change. Default value is 1.
##  q.value - an integer argument. It is a threshold for FDR-adjusted p-value. Default value is 0.05.
##  log.scale - a logical argument specifying that splicing distances should be log10 transformed.
#   Default value is TRUE.
##  plot - a logical argument specifying that distribution of the splicing distances should be plotted.
#   Default value is TRUE.
splDistance = function(input, log.FC = 1, q.value = 0.05, log.scale = "TRUE", plot = "TRUE"){
## Loading and preparation of the input data.
#  The input file should be standard txt-file in tab-delimited format.
#  This file should include seven mandatory fields:
#  i) Event_ID (ID of the splicing event);
#  ii) seqnames (name of chromosome or scaffold with prefix "chr");
#  iii) start (start genomic coordinate of the exon-exon junction);
#  iv) end (end genomic coordinate of the exon-exon junction);
#  v) strand (strand information about exon-exon junction).
#  vi) logFC (log2 fold change in usage of the splicing event between two conditions);
#  vii) q_value (FDR-adjusted p-value for logFC).
cat("Loading and preparation of the input data.\n")
flush.console()
EEJ = read.table(file = input, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
EEJ.dif = EEJ[abs(EEJ$logFC) >= log.FC & EEJ$q_value < q.value, ]
EEJ.no.dif = EEJ[!EEJ$Event_ID %in% EEJ.dif$Event_ID, ]
setClass("splicingDistances",
         representation = representation(data = "data.frame",
                                         distances = "numeric",
                                         summary = "table",
                                         density = "density"))
Arguments = data.frame(cbind(input, log.FC, q.value, log.scale), stringsAsFactors = FALSE)
colnames(Arguments) = c("Input.file", "log.FC.threshold", "q.value.threshold", "log.scale")
Arguments$log.FC.threshold = as.numeric(Arguments$log.FC.threshold)
Arguments$q.value.threshold = as.numeric(Arguments$q.value.threshold)
## Calculation of the splicing distances (with log10 transformation).
if (log.scale == "TRUE"){
    cat("Calculation of the splicing distances (with log10 transformation).\n")
    flush.console()
    splDis.dif = EEJ.dif$end - EEJ.dif$start - 1
    EEJ.dif$splicing.distances.nucleotides = splDis.dif
    splDis.dif = log(splDis.dif, 10)
    EEJ.dif$splicing.distances.log = splDis.dif
    splDis.no.dif = EEJ.no.dif$end - EEJ.no.dif$start - 1
    EEJ.no.dif$splicing.distances.nucleotides = splDis.no.dif
    splDis.no.dif = log(splDis.no.dif, 10)
    EEJ.no.dif$splicing.distances.log = splDis.no.dif
    ## Consolidation of the results in an object of class splicingDistances.
    #  This object contains input data as well as all calculation results.
    cat("Consolidation of the results in an object of class splicingDistances.\n")
    flush.console()
    splDis = list("Differential splicing events" = new("splicingDistances",
                                                       data = EEJ.dif,
                                                       distances = splDis.dif,
                                                       summary = summary(splDis.dif),
                                                       density = density(splDis.dif)),
                  "Non-differential splicing events" = new("splicingDistances",
                                                           data = EEJ.no.dif,
                                                           distances = splDis.no.dif,
                                                           summary = summary(splDis.no.dif),
                                                           density = density(splDis.no.dif)),
                  "Mann-Whitney-Wilcoxon Test" = wilcox.test(EEJ.dif$splicing.distances.nucleotides,
                                                             EEJ.no.dif$splicing.distances.nucleotides),
                  "Arguments" = Arguments)
}
## Calculation of the splicing distances (without log10 transformation).
if (log.scale == "FALSE"){
    cat("Calculation of the splicing distances (without log10 transformation).\n")
    flush.console()
    splDis.dif = EEJ.dif$end - EEJ.dif$start - 1
    EEJ.dif$splicing.distances.nucleotides = splDis.dif
    splDis.no.dif = EEJ.no.dif$end - EEJ.no.dif$start - 1
    EEJ.no.dif$splicing.distances.nucleotides = splDis.no.dif
    ## Consolidation of the results in an object of class splicingDistances.
    cat("Consolidation of the results in an object of class splicingDistances.\n")
    flush.console()
    splDis = list("Differential splicing events" = new("splicingDistances",
                                                       data = EEJ.dif,
                                                       distances = splDis.dif,
                                                       summary = summary(splDis.dif),
                                                       density = density(splDis.dif)),
                  "Non-differential splicing events" = new("splicingDistances",
                                                           data = EEJ.no.dif,
                                                           distances = splDis.no.dif,
                                                           summary = summary(splDis.no.dif),
                                                           density = density(splDis.no.dif)),
                  "Mann-Whitney-Wilcoxon Test" = wilcox.test(EEJ.dif$splicing.distances.nucleotides,
                                                             EEJ.no.dif$splicing.distances.nucleotides),
                  "Arguments" = Arguments)
}
if (plot == "TRUE"){
    ## Plotting results (with log10 transformation).
    if (log.scale == "TRUE"){
        cat("Plotting results (with log10 transformation).\n")
        flush.console()
        windowsFonts(A = windowsFont("Arial"))
        par(cex.main = 1.2,
            cex.axis = 1.2,
            cex.lab = 1.2,
            lwd = 2,
            family = "A",
            ann = FALSE)
        plot(density(splDis.dif),
             ylim = c(0, max(max(density(splDis.dif)$y), max(density(splDis.no.dif)$y))),
             col = rgb(0, 114, 178, maxColorValue = 255))
        lines(density(splDis.no.dif),
              col = rgb(213, 94, 0, maxColorValue = 255))
        title(main = "Distribution of splicing distances",
              ylab = "Density",
              xlab = expression("Splicing distances, log"[10]*" (nucleotides)"))
        legend(x = "topright",
               y = NULL,
               legend = c("differential splicing events", "non-differential splicing events"),
               col = c(rgb(0, 114, 178, maxColorValue = 255), rgb(213, 94, 0, maxColorValue = 255)),
               lty = 1,
               lwd = 2,
               cex = 1.1,
               bty = "n")
    }
    ## Plotting results (without log10 transformation).
    if (log.scale == "FALSE"){
        cat("Plotting results (without log10 transformation).\n")
        flush.console()
        options(scipen = 999)
        windowsFonts(A = windowsFont("Arial"))
        par(mfrow = c(1, 2),
            cex.main = 1.1,
            cex.axis = 1.1,
            cex.lab = 1.1,
            lwd = 2,
            family = "A")
        hist(x = splDis.dif,
             col = rgb(0, 114, 178, maxColorValue = 255),
             main = "Histogram of splicing distances, \ndifferential splicing events",
             xlab = "Splicing distances, nucleotides",
             ylab = "Number of events")
        hist(x = splDis.no.dif,
             col = rgb(213, 94, 0, maxColorValue = 255),
             main = "Histogram of splicing distances, \nnon-differential splicing events",
             xlab = "Splicing distances, nucleotides",
             ylab = "Number of events")
    }
}
return(splDis)
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/spldistance.r")
#   inputData = "GVV/Experimental_junctions.txt"
#   res = splDistance(input = inputData)
