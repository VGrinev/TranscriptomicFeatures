setwd("D:/Vasily Grinev")
##  Input files.
eejs = "PRMT5 project, diffEEJs, full list, day 6.txt"
motif = "BP motif, PWM.txt"
##  Auxiliary libraries.
if (!suppressMessages(require(GenomicFeatures))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("GenomicFeatures")
}
if (!suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg38")
}
if (!suppressMessages(require(vioplot))){
    source("http://bioconductor.org/biocLite.R")
    biocLite("vioplot")
}
##  Motif to be analysed.
MOTIF = read.table(file = motif, sep = "\t", header = FALSE, quote = "\"", as.is = TRUE)
MOTIF = data.frame(MOTIF, row.names = c(1))
MOTIF = data.matrix(MOTIF, rownames.force = TRUE)
##  Loading of the experimental exon-exon junctions.
expEEJs = read.table(file = eejs, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
difEEJs.d = expEEJs[expEEJs$q_value < 0.1 & expEEJs$logFC < 0, c(-2, -7:-9)]
difEEJs.d = difEEJs.d[!duplicated(difEEJs.d), ]
difEEJs.u = expEEJs[expEEJs$q_value < 0.1 & expEEJs$logFC > 0, c(-2, -7:-9)]
difEEJs.u = difEEJs.u[!duplicated(difEEJs.u), ]
n.difEEJs = expEEJs[expEEJs$q_value >= 0.1, c(-2, -7:-9)]
n.difEEJs = n.difEEJs[!duplicated(n.difEEJs), ]
##  Sub-setting of the 3' splice sites of differential exon-exon junctions, logFC < 0.
difEEJs.d.3p = difEEJs.d[difEEJs.d$strand == "+", ]
difEEJs.d.3p$start = difEEJs.d.3p$end - 100
difEEJs.d.3p$end = difEEJs.d.3p$start + 85
difEEJs.d.3m = difEEJs.d[difEEJs.d$strand == "-", ]
difEEJs.d.3m$start = difEEJs.d.3m$start + 15
difEEJs.d.3m$end = difEEJs.d.3m$start + 85
difEEJs.d.3f = rbind(difEEJs.d.3p, difEEJs.d.3m)
difEEJs.d.3f = difEEJs.d.3f[order(difEEJs.d.3f$eej_id), ]
rownames(difEEJs.d.3f) = NULL
difEEJs.d.3f = makeGRangesFromDataFrame(df = difEEJs.d.3f, keep.extra.columns = FALSE)
difEEJs.d.3f = getSeq(x = Hsapiens, difEEJs.d.3f)
##  Motif analysis for the 3' splice sites of differential exon-exon junctions.
scores_difEEJs.d.3f = list()
for (i in 1:length(difEEJs.d.3f)){
     hits = matchPWM(pwm = MOTIF, subject = difEEJs.d.3f[[i]], min.score = "0%", with.score = TRUE)
     hits = cbind(start(hits), mcols(hits)$score)
     scores_difEEJs.d.3f[[i]] = hits
}
##  Sub-setting of the 3' splice sites of differential exon-exon junctions, logFC > 0.
difEEJs.u.3p = difEEJs.u[difEEJs.u$strand == "+", ]
difEEJs.u.3p$start = difEEJs.u.3p$end - 100
difEEJs.u.3p$end = difEEJs.u.3p$start + 85
difEEJs.u.3m = difEEJs.u[difEEJs.u$strand == "-", ]
difEEJs.u.3m$start = difEEJs.u.3m$start + 15
difEEJs.u.3m$end = difEEJs.u.3m$start + 85
difEEJs.u.3f = rbind(difEEJs.u.3p, difEEJs.u.3m)
difEEJs.u.3f = difEEJs.u.3f[order(difEEJs.u.3f$eej_id), ]
rownames(difEEJs.u.3f) = NULL
difEEJs.u.3f = makeGRangesFromDataFrame(df = difEEJs.u.3f, keep.extra.columns = FALSE)
difEEJs.u.3f = getSeq(x = Hsapiens, difEEJs.u.3f)
##  Motif analysis for the 3' splice sites of differential exon-exon junctions.
scores_difEEJs.u.3f = list()
for (i in 1:length(difEEJs.u.3f)){
     hits = matchPWM(pwm = MOTIF, subject = difEEJs.u.3f[[i]], min.score = "0%", with.score = TRUE)
     hits = cbind(start(hits), mcols(hits)$score)
     scores_difEEJs.u.3f[[i]] = hits
}
##  Sub-setting of the 3' splice sites of non-differential exon-exon junctions.
n.difEEJs.3p = n.difEEJs[n.difEEJs$strand == "+", ]
n.difEEJs.3p$start = n.difEEJs.3p$end - 100
n.difEEJs.3p$end = n.difEEJs.3p$start + 85
n.difEEJs.3m = n.difEEJs[n.difEEJs$strand == "-", ]
n.difEEJs.3m$start = n.difEEJs.3m$start + 15
n.difEEJs.3m$end = n.difEEJs.3m$start + 85
n.difEEJs.3f = rbind(n.difEEJs.3p, n.difEEJs.3m)
n.difEEJs.3f = n.difEEJs.3f[order(n.difEEJs.3f$eej_id), ]
rownames(n.difEEJs.3f) = NULL
n.difEEJs.3f = makeGRangesFromDataFrame(df = n.difEEJs.3f, keep.extra.columns = FALSE)
n.difEEJs.3f = getSeq(x = Hsapiens, n.difEEJs.3f)
##  Motif analysis for the 3' splice sites of non-differential exon-exon junctions.
scores_n.difEEJs.3f = list()
for (i in 1:length(n.difEEJs.3f)){
     hits = matchPWM(pwm = MOTIF, subject = n.difEEJs.3f[[i]], min.score = "0%", with.score = TRUE)
     hits = cbind(start(hits), mcols(hits)$score)
     scores_n.difEEJs.3f[[i]] = hits
}
##  Post-processing of the motif analysis results.
#   Distribution of scores.
score_difEEJs.d.3f = do.call(rbind, scores_difEEJs.d.3f)
score_difEEJs.u.3f = do.call(rbind, scores_difEEJs.u.3f)
score_n.difEEJs.3f = do.call(rbind, scores_n.difEEJs.3f)
plot(x = density(score_n.difEEJs.3f[, 2], bw = 0.2),
     main = "",
     xlab = "Branchpoint motif, score",
     ylab = "Density")
lines(density(score_difEEJs.d.3f[, 2]), type = "l", col = "blue")
lines(density(score_difEEJs.u.3f[, 2]), type = "l", col = "red")
legend(x = "topright",
       border = "white",
       lty = 1,
       col = c("black", "blue", "red"),
       legend = c("non-diffEEJs", "diffEEJs, logFC < 0", "diffEEJs, logFC > 0"))
wilcox.test(score_n.difEEJs.3f, score_difEEJs.d.3f, alternative = "greater")
wilcox.test(score_n.difEEJs.3f, score_difEEJs.u.3f, alternative = "greater")
#   Frequency of motif above of the given threshold, 3' splice sites.
thr = 0.8
start = 1
end = 87 - dim(MOTIF)[2]
number = 86 - dim(MOTIF)[2] + 1
freq_difEEJs.d.3f = lapply(scores_difEEJs.d.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, ]
                                                            return(nrow(x)/number)})
freq_difEEJs.d.3f = do.call(rbind, freq_difEEJs.d.3f)
freq_difEEJs.u.3f = lapply(scores_difEEJs.u.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, ]
                                                            return(nrow(x)/number)})
freq_difEEJs.u.3f = do.call(rbind, freq_difEEJs.u.3f)
freq_n.difEEJs.3f = lapply(scores_n.difEEJs.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, ]
                                                            return(nrow(x)/number)})
freq_n.difEEJs.3f = do.call(rbind, freq_n.difEEJs.3f)
plot(density(freq_n.difEEJs.3f, bw = 0.007),
     main = "",
     xlab = "Frequency",
     ylab = "Density",
     xlim = c(0, 0.25))
lines(density(freq_difEEJs.d.3f), type = "l", col = "blue")
lines(density(freq_difEEJs.u.3f), type = "l", col = "red")
legend(x = "topright",
       border = "white",
       lty = 1,
       col = c("black", "blue", "red"),
       legend = c("non-diffEEJs", "diffEEJs, logFC < 0", "diffEEJs, logFC > 0"))
summary(freq_n.difEEJs.3f)
summary(freq_difEEJs.d.3f)
summary(freq_difEEJs.u.3f)
wilcox.test(freq_n.difEEJs.3f, freq_difEEJs.d.3f, alternative = "greater")
wilcox.test(freq_n.difEEJs.3f, freq_difEEJs.u.3f, alternative = "greater")
boxplot(x = list(freq_n.difEEJs.3f, freq_difEEJs.d.3f, freq_difEEJs.u.3f))
vioplot(freq_n.difEEJs.3f, freq_difEEJs.d.3f, freq_difEEJs.u.3f,
        names = c("non-diffEEJs", "diffEEJs, logFC < 0", "diffEEJs, logFC > 0"))
#   Positional distribution of motif above of the given threshold, 3' splice sites.
thr = 0.8
start = 1
end = 87 - dim(MOTIF)[2]
freq_difEEJs.d.3f = lapply(scores_difEEJs.d.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, 1]})
freq_difEEJs.d.3f = data.frame(table(do.call(c, freq_difEEJs.d.3f))/length(scores_difEEJs.d.3f))
freq_difEEJs.d.3f$Var1 = as.numeric(freq_difEEJs.d.3f$Var1)
freq_difEEJs.u.3f = lapply(scores_difEEJs.u.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, 1]})
freq_difEEJs.u.3f = data.frame(table(do.call(c, freq_difEEJs.u.3f))/length(scores_difEEJs.u.3f))
freq_difEEJs.u.3f$Var1 = as.numeric(freq_difEEJs.u.3f$Var1)
freq_n.difEEJs.3f = lapply(scores_n.difEEJs.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, 1]})
freq_n.difEEJs.3f = data.frame(table(do.call(c, freq_n.difEEJs.3f))/length(scores_n.difEEJs.3f))
freq_n.difEEJs.3f$Var1 = as.numeric(freq_n.difEEJs.3f$Var1)
plot(freq_n.difEEJs.3f, type = "l")
lines(freq_difEEJs.d.3f, type = "l", col = "blue")
lines(freq_difEEJs.u.3f, type = "l", col = "red")
##  Analysis of the RIs.
RIs.diff = read.table(file = "RIs, diffEU.txt", sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
RIs.n.diff = read.table(file = "RIs, non-diffEU.txt", sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
##  Sub-setting of the 3' splice sites of differential exon-exon junctions, logFC < 0.
difRIs.3p = RIs.diff[RIs.diff$strand == "+", ]
difRIs.3p$start = difRIs.3p$end - 100
difRIs.3p$end = difRIs.3p$start + 85
difRIs.3m = RIs.diff[RIs.diff$strand == "-", ]
difRIs.3m$start = difRIs.3m$start + 15
difRIs.3m$end = difRIs.3m$start + 85
difRIs.3f = rbind(difRIs.3p, difRIs.3m)
difRIs.3f = difRIs.3f[order(difRIs.3f$exon_id), ]
rownames(difRIs.3f) = NULL
difRIs.3f = makeGRangesFromDataFrame(df = difRIs.3f, keep.extra.columns = FALSE)
difRIs.3f = getSeq(x = Hsapiens, difRIs.3f)
##  Motif analysis for the 3' splice sites of differential exon-exon junctions.
scores_difRIs.3f = list()
for (i in 1:length(difRIs.3f)){
     hits = matchPWM(pwm = MOTIF, subject = difRIs.3f[[i]], min.score = "0%", with.score = TRUE)
     hits = cbind(start(hits), mcols(hits)$score)
     scores_difRIs.3f[[i]] = hits
}
##  Sub-setting of the 3' splice sites of non-differential exon-exon junctions.
n.difRIs.3p = RIs.n.diff[RIs.n.diff$strand == "+", ]
n.difRIs.3p$start = n.difRIs.3p$end - 100
n.difRIs.3p$end = n.difRIs.3p$start + 85
n.difRIs.3m = RIs.n.diff[RIs.n.diff$strand == "-", ]
n.difRIs.3m$start = n.difRIs.3m$start + 15
n.difRIs.3m$end = n.difRIs.3m$start + 85
n.difRIs.3f = rbind(n.difRIs.3p, n.difRIs.3m)
n.difRIs.3f = n.difRIs.3f[order(n.difRIs.3f$exon_id), ]
rownames(n.difRIs.3f) = NULL
n.difRIs.3f = makeGRangesFromDataFrame(df = n.difRIs.3f, keep.extra.columns = FALSE)
n.difRIs.3f = getSeq(x = Hsapiens, n.difRIs.3f)
##  Motif analysis for the 3' splice sites of non-differential exon-exon junctions.
scores_n.difRIs.3f = list()
for (i in 1:length(n.difRIs.3f)){
     hits = matchPWM(pwm = MOTIF, subject = n.difRIs.3f[[i]], min.score = "0%", with.score = TRUE)
     hits = cbind(start(hits), mcols(hits)$score)
     scores_n.difRIs.3f[[i]] = hits
}
##  Post-processing of the motif analysis results.
#   Distribution of scores.
score_difRIs.3f = do.call(rbind, scores_difRIs.3f)
score_n.difRIs.3f = do.call(rbind, scores_n.difRIs.3f)
plot(density(score_n.difRIs.3f[, 2]),
     main = "",
     xlab = "Branchpoint motif, score",
     ylab = "Density")
lines(density(score_difRIs.3f[, 2]), type = "l", col = "blue")
legend(x = "topright",
       border = "white",
       lty = 1,
       col = c("black", "blue"),
       legend = c("non-diffEEJs", "diffEEJs"))
wilcox.test(score_n.difRIs.3f, score_difRIs.3f, alternative = "greater")
#   Frequency of motif above of the given threshold, 3' splice sites.
thr = 0.72
start = 1
end = 87 - dim(MOTIF)[2]
number = 86 - dim(MOTIF)[2] + 1
freq_difRIs.3f = lapply(scores_difRIs.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, ]
                                                      return(nrow(x)/number)})
freq_difRIs.3f = do.call(rbind, freq_difRIs.3f)
freq_n.difRIs.3f = lapply(scores_n.difRIs.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, ]
                                                          return(nrow(x)/number)})
freq_n.difRIs.3f = do.call(rbind, freq_n.difRIs.3f)
plot(x = density(freq_n.difRIs.3f),
     main = "",
     xlab = "Frequency",
     ylab = "Density",
     xlim = c(0, 0.25))
lines(density(freq_difRIs.3f), type = "l", col = "blue")
legend(x = "topright",
       border = "white",
       lty = 1,
       col = c("black", "blue"),
       legend = c("non-diffEEJs", "diffEEJs"))
summary(freq_n.difRIs.3f)
summary(freq_difRIs.3f)
wilcox.test(freq_n.difRIs.3f, freq_difRIs.3f, alternative = "less")
boxplot(x = list(freq_n.difRIs.3f, freq_difRIs.3f))
vioplot(freq_n.difRIs.3f, freq_difRIs.3f,
        names = c("non-diffEEJs", "diffEEJs"))
#   Positional distribution of motif above of the given threshold, 3' splice sites.
thr = 0
start = 1
end = 87 - dim(MOTIF)[2]
freq_difRIs.3f = lapply(scores_difRIs.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, 1]})
freq_difRIs.3f = data.frame(table(do.call(c, freq_difRIs.3f))/length(scores_difRIs.3f))
freq_difRIs.3f$Var1 = as.numeric(freq_difRIs.3f$Var1)
freq_n.difRIs.3f = lapply(scores_n.difRIs.3f, function(x){x = x[x[, 1] >= start & x[, 1] <= end & x[, 2] >= thr, 1]})
freq_n.difRIs.3f = data.frame(table(do.call(c, freq_n.difRIs.3f))/length(scores_n.difRIs.3f))
freq_n.difRIs.3f$Var1 = as.numeric(freq_n.difRIs.3f$Var1)
plot(freq_n.difRIs.3f, type = "l")
lines(freq_difRIs.3f, type = "l", col = "blue")
