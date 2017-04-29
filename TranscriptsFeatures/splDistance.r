#####################################################################
## High-level R function                                           ##
## for consolidation of global statistics on splicing distances.   ##
## (c) GNU GPL Vasily V. Grinev, 2017. grinev_vv[at]bsu.by         ##
#####################################################################
## Arguments of function:
# input - path to folder and name of the input txt-file in tab-delimited format.
# log.FC - an integer argument. It is a threshold for log2 fold change. Default value is 1.
# q.value - an integer argument. It is a threshold for FDR-adjusted p-value. Default value is 0.05.
# log.scale - a logical argument specifying that splicing distances should be log10 transformed. Default value is TRUE.
# plot - a logical argument specifying that distribution of the splicing distances should be plotted. Default value is TRUE.
splDistance = function(input, log.FC = 1, q.value = 0.05, log.scale = "TRUE", plot = "TRUE"){
	## Loading and preparation of the input data.
	# The input file should be standard txt-file in tab-delimited format. This file should include five mandatory fields:
	# 1) Event_ID (ID of the splicing event);
	# 2) Start (start genomic coordinate of the splicing event);
	# 3) End (end genomic coordinate of the splicing event);
	# 4) logFC (log2 fold change in usage of the splicing event between two conditions);
	# 5) q_value (FDR-adjusted p-value for logFC).
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
		splDis.dif = EEJ.dif$End - EEJ.dif$Start - 1
		EEJ.dif$splicing.distances.nucleotides = splDis.dif
		splDis.dif = log(splDis.dif, 10)
		EEJ.dif$splicing.distances.log = splDis.dif
		splDis.no.dif = EEJ.no.dif$End - EEJ.no.dif$Start - 1
		EEJ.no.dif$splicing.distances.nucleotides = splDis.no.dif
		splDis.no.dif = log(splDis.no.dif, 10)
		EEJ.no.dif$splicing.distances.log = splDis.no.dif
		## Consolidation of the results in an object of class splicingDistances.
		# This object contains input data as well as all calculation results.
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
					  "Mann-Whitney-Wilcoxon Test" = wilcox.test(EEJ.dif$splicing.distances.nucleotides, EEJ.no.dif$splicing.distances.nucleotides),
					  "Arguments" = Arguments)
	}
	## Calculation of the splicing distances (without log10 transformation).
	if (log.scale == "FALSE"){
		splDis.dif = EEJ.dif$End - EEJ.dif$Start - 1
		EEJ.dif$splicing.distances.nucleotides = splDis.dif
		splDis.no.dif = EEJ.no.dif$End - EEJ.no.dif$Start - 1
		EEJ.no.dif$splicing.distances.nucleotides = splDis.no.dif
		## Consolidation of the results in an object of class splicingDistances.
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
					  "Mann-Whitney-Wilcoxon Test" = wilcox.test(EEJ.dif$splicing.distances.nucleotides, EEJ.no.dif$splicing.distances.nucleotides),
					  "Arguments" = Arguments)
	}
	if (plot == "TRUE"){
		## Plotting results (with log10 transformation).
		if (log.scale == "TRUE"){
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
## A simple example of function use.
# source("/media/gvv/NGS_Illumina_HiSeq/Software/splDistance.r")
# inputData = "/media/gvv/NGS_Illumina_HiSeq/PHF5A dataset, limma.diffSplice, exon-exon junctions.gtf"
# res = splDistance(input = inputData)

inputData = "D:/Vasily Grinev/PHF5A dataset, JunctionSeq, exon-exon junctions.txt"
res = splDistance(input = inputData)#, log.FC = 1, q.value = 0.05, log.scale = "TRUE", plot = "TRUE")
l.dif = max(length(res[[1]]@data$splicing.distances.nucleotides), length(res[[1]]@density$x), length(res[[1]]@density$y))
res.dif = data.frame(splicing.distances.nucleotides = c(res[[1]]@data$splicing.distances.nucleotides, rep("", l.dif - length(res[[1]]@data$splicing.distances.nucleotides))),
					 splicing.distances = c(res[[1]]@density$x, rep("", l.dif - length(res[[1]]@density$x))),
					 density = c(res[[1]]@density$y, rep("", l.dif - length(res[[1]]@density$y))))
l.no.dif = max(length(res[[2]]@data$splicing.distances.nucleotides), length(res[[2]]@density$x), length(res[[2]]@density$y))
res.no.dif = data.frame(splicing.distances.nucleotides = c(res[[2]]@data$splicing.distances.nucleotides, rep("", l.no.dif - length(res[[2]]@data$splicing.distances.nucleotides))),
					 splicing.distances = c(res[[2]]@density$x, rep("", l.no.dif - length(res[[2]]@density$x))),
					 density = c(res[[2]]@density$y, rep("", l.no.dif - length(res[[2]]@density$y))))
write.table(res.dif,
			file = "D:/Vasily Grinev/Differential splicing events, splicing distances.txt",
			sep = "\t",
			quote = FALSE,
			col.names = TRUE,
			row.names = FALSE)
write.table(res.no.dif,
			file = "D:/Vasily Grinev/Non-differential splicing events, splicing distances.txt",
			sep = "\t",
			quote = FALSE,
			col.names = TRUE,
			row.names = FALSE)
