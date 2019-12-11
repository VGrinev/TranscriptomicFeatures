### Random forest-based classification of the splicing events.
rm(list = ls())
setwd("D:/Vasily Grinev")
## Required libraries.
library(caret)
library(randomForest)
library(randomForest.ddR)
library(corrplot)
library(matrixStats)
## Loading of the primary features matrix.
SpliceEvents = read.table(file = "De novo assembled transcripts, all features.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "\"")
## Filtration of features.
workData = SpliceEvents[, -12:-13][, -1:-9]
nz = nearZeroVar(x = workData, freqCut = 95/5, uniqueCut = 10, saveMetrics = FALSE, foreach = FALSE, allowParallel = TRUE)
workData = workData[, -nz]
corVar = cor(workData, method = "spearman")
highCor = findCorrelation(x = corVar, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
workData = workData[, -highCor]
workData = cbind(SpliceEvents[, c(1:9, 12:13)], workData)
write.table(workData, file = "De novo assembled transcripts, filtered features.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
## Loading of the filtered features matrix.
SpliceEvents = read.table(file = "De novo assembled transcripts, filtered features.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "\"")
SpliceEvents$Class = as.factor(SpliceEvents$Class)
## Computation settings.
Events = table(SpliceEvents$Class == "DifSpl")[[2]]
ntree = 1000
ntry = 100
nExecutor = 4
run = 100
index = "Gini"
## Features importances (with Breiman-Cutler index (mean decrease accuracy or permutation importance) and Gini impurity index (mean decrease Gini or Gini importance)).
DifSpl = SpliceEvents[which(SpliceEvents$Class == "DifSpl"), ]
importanceMatrix = matrix(nrow = ncol(SpliceEvents) - 3, ncol = ntry * 4)
for (i in 1:ntry){
	# Formation of the balanced data set.
	NoDifSpl = SpliceEvents[sample(which(SpliceEvents$Class == "NoDifSpl"), Events), ]
	BalancedDataSet = rbind(DifSpl, NoDifSpl)[, -3][, -1]
	# Training of the random forest model.
	randomForestModel = drandomForest(Class ~ ., data = BalancedDataSet, type = classification, ntree = ntree, importance = TRUE, nExecutor = nExecutor, trace = FALSE)
	importanceOfFeatures = cbind(rownames(importance(randomForestModel)), importance(randomForestModel))
	importanceOfFeatures = importanceOfFeatures[order(importanceOfFeatures[, 1]), ]
	importanceMatrix[, i] = importanceOfFeatures[, 2]
	importanceMatrix[, i + ntry] = importanceOfFeatures[, 3]
	importanceMatrix[, i + (2 * ntry)] = importanceOfFeatures[, 4]
	importanceMatrix[, i + (3 * ntry)] = importanceOfFeatures[, 5]
	if(i%%1 == 0){
		cat(i, "run of", ntry, "has been finished...\n")
		flush.console()
	}
}
# Converting character matrix into numeric matrix.
class(importanceMatrix) = "numeric"
# Robustness of indices.
r = cor(importanceMatrix)
corrplot(r, method = "color", type = "upper")
# Summarization and saving of results.
Imp = cbind(rowMeans(importanceMatrix[, 1:ntry]), rowSds(importanceMatrix[, 1:ntry]),
			rowMeans(importanceMatrix[, (ntry + 1):(2 * ntry)]), rowSds(importanceMatrix[, (ntry + 1):(2 * ntry)]),
			rowMeans(importanceMatrix[, (2 * ntry + 1):(3 * ntry)]), rowSds(importanceMatrix[, (2 * ntry + 1):(3 * ntry)]),
			rowMeans(importanceMatrix[, (3 * ntry + 1):(4 * ntry)]), rowSds(importanceMatrix[, (3 * ntry + 1):(4 * ntry)]))
Imp = cbind(sort(colnames(SpliceEvents[-1:-3])), Imp)
colnames(Imp) = c("Feature", "Class_DifSpl_mean", "Class_DifSpl_SD", "Class_NoDifSpl_mean", "Class_NoDifSpl_SD", "Mean_decrease_in_accuracy_mean", "Mean_decrease_in_accuracy_SD",
				  "Gini_index_mean", "Gini_index_SD")
write.table(Imp, file = "Splicing events, importance of features.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
## Features selection with ranked features.
Features = read.table(file = "Splicing events, importance of features.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "\"")
Features = Features[order(-Features$Mean_decrease_in_accuracy_mean), ]
MDA = Features$Feature
Features = Features[order(-Features$Gini_index_mean), ]
Gini = Features$Feature
SubsetSize = seq(1, 100, 1)
Accuracy = matrix(ncol = run * 3, nrow = length(SubsetSize))
DifSpl = SpliceEvents[which(SpliceEvents$Class == "DifSpl"), ]
t1 = Sys.time()
for (j in 1:run){
	# Formation of the balanced data set.
	NoDifSpl = SpliceEvents[sample(which(SpliceEvents$Class == "NoDifSpl"), Events), ]
	BalancedDataSet = rbind(DifSpl, NoDifSpl)[, -3][, -1]
	Class = BalancedDataSet$Class
	for (k in 1:length(SubsetSize)){
		# Subset of data.
		if (index == "MDA"){
			SubsetFeatures = MDA[1:SubsetSize[k]]
			SubsetData = cbind(Class, BalancedDataSet[SubsetFeatures])
		}
		if (index == "Gini"){
			SubsetFeatures = Gini[1:SubsetSize[k]]
			SubsetData = cbind(Class, BalancedDataSet[SubsetFeatures])
		}
		# Data splitting.
		tr_vs_te = sample(2, nrow(SubsetData), replace = TRUE, prob = c(0.7, 0.3))
		train = SubsetData[tr_vs_te == 1, ]
		test = SubsetData[tr_vs_te == 2, ]
		# Model training.
		randomForestModel = drandomForest(Class ~ ., data = train, type = classification, ntree = ntree, importance = TRUE, nExecutor = nExecutor, trace = FALSE)
		# Model testing.
		pred = predict(object = randomForestModel, newdata = test, type = "response", norm.votes = TRUE, predict.all = FALSE, proximity = FALSE, nodes = FALSE)
		generalizedPrecision = round(length(subset(as.numeric(test[, 1]) - as.numeric(pred),
									 as.numeric(test[, 1]) - as.numeric(pred) == 0))/length(test[, 1]), digits = 3)
		Class_DifSpl_Precision = 1 - (round(length(subset(as.numeric(test[, 1]) - as.numeric(pred),
											as.numeric(test[, 1]) - as.numeric(pred) == -1))/summary(test$Class)[[1]], digits = 3))
		Class_NoDifSpl_Precision = 1 - (round(length(subset(as.numeric(test[, 1]) - as.numeric(pred),
											  as.numeric(test[, 1]) - as.numeric(pred) == 1))/summary(test$Class)[[2]], digits = 3))
		# Results collection.
		Accuracy[k, j] = generalizedPrecision
		Accuracy[k, j + run] = Class_DifSpl_Precision
		Accuracy[k, j + (2 * run)] = Class_NoDifSpl_Precision
		if(k%%1 == 0){
			cat("Repeat", j, "...,", "subset size", SubsetSize[k], "...\n")
			flush.console()
		}
	}
}
t2 = Sys.time()
print(t2 - t1)
# Converting character matrix into numeric matrix.
class(Accuracy) = "numeric"
# Summarization and saving of results.
Accur = cbind(SubsetSize,
			  round(rowMeans(Accuracy[, 1:run]), digits = 3), round(rowSds(Accuracy[, 1:run]), digits = 3),
			  round(rowMeans(Accuracy[, (run + 1):(2 * run)]), digits = 3), round(rowSds(Accuracy[, (run + 1):(2 * run)]), digits = 3),
			  round(rowMeans(Accuracy[, ((2 * run) + 1):(3 * run)]), digits = 3), round(rowSds(Accuracy[, ((2 * run) + 1):(3 * run)]), digits = 3))
colnames(Accur) = c("Number_of_features", "Gener_mean", "Gener_SD", "Class_DifSpl_mean", "Class_DifSpl_SD", "Class_NoDifSpl_mean", "Class_NoDifSpl_SD")
write.table(Accur, file = "Splicing events, accuracy of the random forest models, index Gini.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
## Features selection with recursive feature selection algorithm.
Features = read.table(file = "Splicing events, importance of features.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "\"")
Features = Features[order(-Features$Mean_decrease_in_accuracy_mean), ]
MDA = Features$Feature[1:100]
Features = Features[order(-Features$Gini_index_mean), ]
Gini = Features$Feature[1:100]
DifSpl = SpliceEvents[which(SpliceEvents$Class == "DifSpl"), ]
ReFeSe_Features = list()
for (l in 1:run){
	# Formation of the balanced data set.
	NoDifSpl = SpliceEvents[sample(which(SpliceEvents$Class == "NoDifSpl"), Events), ]
	BalancedDataSet = rbind(DifSpl, NoDifSpl)[, -3][, -1]
	Class = BalancedDataSet$Class
	# Subset of data.
	if (index == "MDA"){
		SubsetData = cbind(Class, BalancedDataSet[MDA])
	}
	if (index == "Gini"){
		SubsetData = cbind(Class, BalancedDataSet[Gini])
	}
	# Recursive feature selection.
	ReFeSe = rfe(x = SubsetData[, -1], y = SubsetData[, 1], sizes = 1:ncol(SubsetData[, -1]), metric = "Accuracy",
				 rfeControl = rfeControl(functions = rfFuncs, method = "cv", number = 5, verbose = FALSE, allowParallel = TRUE))
	ReFeSe_Features[[l]] = predictors(ReFeSe)
	if(l%%1 == 0){
		cat("Repeat", l, "has been finished...\n")
		flush.console()
	}
}
# Summarization and saving of results.
ListFeatures = sort(table(unlist(ReFeSe_Features)))
ListFeatures = cbind(rownames(ListFeatures), ListFeatures/run)
colnames(ListFeatures) = c("Feature", "Frequency")
write.table(ListFeatures, file = "Splicing events, recursively selected features, index Gini.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
## Final classification.
Features = read.table(file = "Splicing events, importance of features, short list.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "\"")
Features = Features$Feature
DifSpl = SpliceEvents[which(SpliceEvents$Class == "DifSpl"), ]
set.seed(100)
run = 100
Accuracy = matrix(ncol = 3, nrow = run)
for (m in 1:run){
	# Formation of the balanced data set.
	NoDifSpl = SpliceEvents[sample(which(SpliceEvents$Class == "NoDifSpl"), (Events * 3)), ]
	BalancedDataSet = rbind(DifSpl, NoDifSpl)[, -3][, -1]
	Class = BalancedDataSet$Class
	# Subset of data.
	SubsetData = cbind(Class, BalancedDataSet[Features])
	# Data splitting.
	tr_vs_te = sample(2, nrow(SubsetData), replace = TRUE, prob = c(0.7, 0.3))
	train = SubsetData[tr_vs_te == 1, ]
	test = SubsetData[tr_vs_te == 2, ]
	# Model training.
	randomForestModel = drandomForest(Class ~ ., data = train, type = classification, ntree = 1000, importance = TRUE, nExecutor = 4, trace = FALSE)
	# Model testing.
	pred = predict(object = randomForestModel, newdata = test, type = "response", norm.votes = TRUE, predict.all = FALSE, proximity = FALSE, nodes = FALSE)
	generalizedPrecision = round(length(subset(as.numeric(test[, 1]) - as.numeric(pred),
								 as.numeric(test[, 1]) - as.numeric(pred) == 0))/length(test[, 1]), digits = 3)
	Class_DifSpl_Precision = 1 - (round(length(subset(as.numeric(test[, 1]) - as.numeric(pred),
										as.numeric(test[, 1]) - as.numeric(pred) == -1))/summary(test$Class)[[1]], digits = 3))
	Class_NoDifSpl_Precision = 1 - (round(length(subset(as.numeric(test[, 1]) - as.numeric(pred),
										  as.numeric(test[, 1]) - as.numeric(pred) == 1))/summary(test$Class)[[2]], digits = 3))
	Accuracy[m, 1] = generalizedPrecision
	Accuracy[m, 2] = Class_DifSpl_Precision
	Accuracy[m, 3] = Class_NoDifSpl_Precision
	if(m%%1 == 0){
		cat("Run", m, "has been finished...\n")
		flush.console()
	}
}
# Converting character matrix into numeric matrix.
class(Accuracy) = "numeric"
# Summarization of results.
colMeans(Accuracy)
[1] 0.68775	0.67578	0.71059
colSds(Accuracy)
[1] 0.05055158	0.09594017	0.09999870
max(Accuracy)
[1] 0.782	0.957	0.938
### Marginal effects of variables.
rm(list = ls())
setwd("D:/Vasily Grinev")
## Required libraries.
library(randomForest)
## Auxiliary function.
equalBins = function(x, N){
	nx = length(x)
	nrepl = floor(nx/N)
	nplus = sample(1:N, nx - nrepl * N)
	nrep = rep(nrepl, N)
	nrep[nplus] = nrepl + 1
	x[order(x)] = rep(seq.int(N), nrep)
	x
}
## Loading of the matrix of filtered features.
SpliceEvents = read.table(file = "Splicing events, all filtered features.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "\"")
SpliceEvents$Class = as.factor(SpliceEvents$Class)
Class_DifSpl = SpliceEvents[which(SpliceEvents$Class == "DifSpl"), ]
## Loading of the short list of features.
Features = read.table(file = "Splicing events, importance of features, short list.txt", header = TRUE, sep = "\t", as.is = TRUE, quote = "\"")$Feature
## Computation settings.
ntree = 1000
ntry = 100
N = 50
# Calculation of marginal effects.
DifSpl_ME = list()
NoDifSpl_ME = list()
for (n in 1:ntry){
	Class_NoDifSpl = SpliceEvents[sample(which(SpliceEvents$Class == "NoDifSpl"), nrow(Class_DifSpl)), ]
	BalancedDataSet = rbind(Class_DifSpl, Class_NoDifSpl)[, -3][, -1]
	Class = BalancedDataSet$Class
	BalancedDataSet = BalancedDataSet[Features]
	randomForestModel = randomForest(Class ~ ., data = BalancedDataSet, type = classification, ntree = ntree, importance = FALSE)
	DifSpl_ME[[n]] = data.frame(partialPlot(x = randomForestModel, pred.data = BalancedDataSet, x.var = FiveSS_Distance_to_RUNX1_RUNX1T1_siMM, which.class = "DifSpl", plot = FALSE))
	NoDifSpl_ME[[n]] = data.frame(partialPlot(x = randomForestModel, pred.data = BalancedDataSet, x.var = FiveSS_Distance_to_RUNX1_RUNX1T1_siMM, which.class = "NoDifSpl", plot = FALSE))
	if(n%%1 == 0){
		cat("Run", n, "has been finished...\n")
		flush.console()
	}
}
DifSpl_ME = do.call(rbind, DifSpl_ME)
NoDifSpl_ME = do.call(rbind, NoDifSpl_ME)
DifSpl_ME = split(DifSpl_ME, equalBins(DifSpl_ME$x, N = N))
x = double()
y = double()
for (p in 1:N){
	x[[p]] = mean(DifSpl_ME[[p]]$x)
	y[[p]] = mean(DifSpl_ME[[p]]$y)
}
DifSpl_ME = cbind(x, y)
NoDifSpl_ME = split(NoDifSpl_ME, equalBins(NoDifSpl_ME$x, N = N))
x = double()
y = double()
for (p in 1:N){
	x[[p]] = mean(NoDifSpl_ME[[p]]$x)
	y[[p]] = mean(NoDifSpl_ME[[p]]$y)
}
NoDifSpl_ME = cbind(x, y)
ME = cbind(DifSpl_ME, NoDifSpl_ME)
colnames(ME) = c("DifSpl_Feature_value", "DifSpl_Marginal_effect", "NoDifSpl_Feature_value", "NoDifSpl_Marginal_effect")
write.table(ME, file = "Splicing events, marginal effects of FiveSS_Distance_to_RUNX1_RUNX1T1_siMM.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
