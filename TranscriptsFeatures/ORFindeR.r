#######################################################################################
## A set of R functions for identification of significant open reading frames        ##
## in nucleotide sequences using multinomial model.                                  ##
## (c) GNU GPL Vasily V. Grinev, 2016. grinev_vv[at]bsu.by                           ##
#######################################################################################
# This code is based on original ideas of Avril Coghlan (Wellcome Trust               #
# Sanger Institute, Cambridge, U.K.) published 10 June 2011 in "A little              #
# book of R for bioinformatics. Release 0.1" (https://media.readthedocs.              #
# org/pdf/a-little-book-of-r-for-bioinformatics/latest/a-little-book-of               #
# -r-for-bioinformatics.pdf                                                           #
#######################################################################################
## Identification of the start and stop codons in a nucleotide sequence.
# The function scans a nucleotide sequence in the search of
# canonical start codon ATG or noncanonical start codons GTG, TTG and CTG
# as well as canonical stop codons TAA, TAG and TGA.
StartAndStopCodons = function(Seq){
codons = c("ATG", "GTG", "TTG", "CTG", "TAA", "TAG", "TGA")
	for (i in 1:7){
		codon = codons[i]
		occurrences = matchPattern(codon, Seq)
		codonpositions = occurrences@ranges@start
		numoccurrences = length(codonpositions)
		if (i == 1){
			positions = codonpositions
			types = rep(codon, numoccurrences)
				}
		else{
			positions = append(positions, codonpositions, after = length(positions))
			types = append(types, rep(codon, numoccurrences), after = length(types))
		}
	}
indices = order(positions)
positions = positions[indices]
types = types[indices]
mylist = list(positions, types)
return(mylist)
}
## Identification of the ATG-start open reading frame(-s) in a nucleotide sequence.
# This function returns all possible ATG-start open reading frames from a nucleotide sequence.
ORFFindingATG = function(Seq){
	require(Biostrings)
	mylist = StartAndStopCodons(Seq)
	positions = mylist[[1]]
	types = mylist[[2]]
	ORFstarts = numeric()
	ORFstops = numeric()
	ORFlengths = numeric()
	numpositions = length(positions)
		if (numpositions >= 2){
			for (i in 1:(numpositions - 1)){
				posi = positions[i]
				typei = types[i]
				found = 0
			while (found == 0){
				for (j in (i + 1):numpositions){
					posj = positions[j]
					typej = types[j]
					posdiff = posj - posi
					posdiffmod3 = posdiff %% 3
	ORFlength = posj - posi + 3
		if ((typei == "ATG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
			numorfs = length(ORFstops)
			usedstop = -1
		if (numorfs > 0){
			for (k in 1:numorfs){
				ORFstopk = ORFstops[k]
			if (ORFstopk == (posj + 2)){usedstop <- 1}
							}
					}
			if (usedstop == -1){
				ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
				ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
							}
	found = 1
		break
						}
	if (j == numpositions){found = 1}
					}
				}
			}
		}
	indices = order(ORFstarts)
	ORFstarts = ORFstarts[indices]
	ORFstops = ORFstops[indices]
	ORFlengths = numeric()
	numorfs = length(ORFstarts)
	for (i in 1:numorfs){
		ORFstart = ORFstarts[i]
		ORFstop = ORFstops[i]
		ORFlength = ORFstop - ORFstart + 1
		ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
	}
	mylist = list(ORFstarts, ORFstops, ORFlengths)
	return(mylist)
}
## Identification of the GTG-start open reading frame(-s) in a nucleotide sequence.
# This function returns all possible GTG-start open reading frames from a nucleotide sequence.
ORFFindingGTG = function(Seq){
	require(Biostrings)
	mylist = StartAndStopCodons(Seq)
	positions = mylist[[1]]
	types = mylist[[2]]
	ORFstarts = numeric()
	ORFstops = numeric()
	ORFlengths = numeric()
	numpositions = length(positions)
		if (numpositions >= 2){
			for (i in 1:(numpositions - 1)){
				posi = positions[i]
				typei = types[i]
				found = 0
			while (found == 0){
				for (j in (i + 1):numpositions){
					posj = positions[j]
					typej = types[j]
					posdiff = posj - posi
					posdiffmod3 = posdiff %% 3
	ORFlength = posj - posi + 3
		if ((typei == "GTG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
			numorfs = length(ORFstops)
			usedstop = -1
		if (numorfs > 0){
			for (k in 1:numorfs){
				ORFstopk = ORFstops[k]
			if (ORFstopk == (posj + 2)){usedstop <- 1}
							}
					}
			if (usedstop == -1){
				ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
				ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
							}
	found = 1
		break
						}
	if (j == numpositions){found = 1}
					}
				}
			}
		}
	indices = order(ORFstarts)
	ORFstarts = ORFstarts[indices]
	ORFstops = ORFstops[indices]
	ORFlengths = numeric()
	numorfs = length(ORFstarts)
	for (i in 1:numorfs){
		ORFstart = ORFstarts[i]
		ORFstop = ORFstops[i]
		ORFlength = ORFstop - ORFstart + 1
		ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
	}
	mylist = list(ORFstarts, ORFstops, ORFlengths)
	return(mylist)
}
## Identification of the TTG-start open reading frame(-s) in a nucleotide sequence.
# This function returns all possible TTG-start open reading frames from a nucleotide sequence.
ORFFindingTTG = function(Seq){
	require(Biostrings)
	mylist = StartAndStopCodons(Seq)
	positions = mylist[[1]]
	types = mylist[[2]]
	ORFstarts = numeric()
	ORFstops = numeric()
	ORFlengths = numeric()
	numpositions = length(positions)
		if (numpositions >= 2){
			for (i in 1:(numpositions - 1)){
				posi = positions[i]
				typei = types[i]
				found = 0
			while (found == 0){
				for (j in (i + 1):numpositions){
					posj = positions[j]
					typej = types[j]
					posdiff = posj - posi
					posdiffmod3 = posdiff %% 3
	ORFlength = posj - posi + 3
		if ((typei == "TTG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
			numorfs = length(ORFstops)
			usedstop = -1
		if (numorfs > 0){
			for (k in 1:numorfs){
				ORFstopk = ORFstops[k]
			if (ORFstopk == (posj + 2)){usedstop <- 1}
							}
					}
			if (usedstop == -1){
				ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
				ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
							}
	found = 1
		break
						}
	if (j == numpositions){found = 1}
					}
				}
			}
		}
	indices = order(ORFstarts)
	ORFstarts = ORFstarts[indices]
	ORFstops = ORFstops[indices]
	ORFlengths = numeric()
	numorfs = length(ORFstarts)
	for (i in 1:numorfs){
		ORFstart = ORFstarts[i]
		ORFstop = ORFstops[i]
		ORFlength = ORFstop - ORFstart + 1
		ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
	}
	mylist = list(ORFstarts, ORFstops, ORFlengths)
	return(mylist)
}
## Identification of the CTG-start open reading frame(-s) in a nucleotide sequence.
# This function returns all possible CTG-start open reading frames from a nucleotide sequence.
ORFFindingCTG = function(Seq){
	require(Biostrings)
	mylist = StartAndStopCodons(Seq)
	positions = mylist[[1]]
	types = mylist[[2]]
	ORFstarts = numeric()
	ORFstops = numeric()
	ORFlengths = numeric()
	numpositions = length(positions)
		if (numpositions >= 2){
			for (i in 1:(numpositions - 1)){
				posi = positions[i]
				typei = types[i]
				found = 0
			while (found == 0){
				for (j in (i + 1):numpositions){
					posj = positions[j]
					typej = types[j]
					posdiff = posj - posi
					posdiffmod3 = posdiff %% 3
	ORFlength = posj - posi + 3
		if ((typei == "CTG" ) && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
			numorfs = length(ORFstops)
			usedstop = -1
		if (numorfs > 0){
			for (k in 1:numorfs){
				ORFstopk = ORFstops[k]
			if (ORFstopk == (posj + 2)){usedstop <- 1}
							}
					}
			if (usedstop == -1){
				ORFstarts = append(ORFstarts, posi, after = length(ORFstarts))
				ORFstops = append(ORFstops, posj + 2, after = length(ORFstops))
							}
	found = 1
		break
						}
	if (j == numpositions){found = 1}
					}
				}
			}
		}
	indices = order(ORFstarts)
	ORFstarts = ORFstarts[indices]
	ORFstops = ORFstops[indices]
	ORFlengths = numeric()
	numorfs = length(ORFstarts)
	for (i in 1:numorfs){
		ORFstart = ORFstarts[i]
		ORFstop = ORFstops[i]
		ORFlength = ORFstop - ORFstart + 1
		ORFlengths = append(ORFlengths, ORFlength, after = length(ORFlengths))
	}
	mylist = list(ORFstarts, ORFstops, ORFlengths)
	return(mylist)
}
## Generation of set of random nucleotide sequences with the same length as original nucleotide sequence using a multinomial model.
generateSeqsWithMultinomialModel = function(inputsequence, n){
	require("seqinr")
	inputsequencevector = s2c(inputsequence)
	mylength = length(inputsequencevector)
	mytable = table(inputsequencevector)
	letters = rownames(mytable)
	numletters = length(letters)
	probabilities = numeric()
		for (i in 1:numletters){
			letter = letters[i]
			count = mytable[[i]]
			probabilities[i] = count/mylength
		}
	seqs = numeric(n)
		for (j in 1:n){
			seq = sample(letters, mylength, rep = TRUE, prob = probabilities)
			seq = c2s(seq)
			seqs[j] = seq
		}
	return(seqs)
}
## Identification of significant open reading frame(-s) in nucleotide sequence.
# The high-level function that integrates work of all above mentioned functions.
# This function has four arguments:
# x - an object of class DNAStringSet which contains a set of nucleotide sequences of interest;
# typeORF - the type of the start codon;
# n - the number of random sequences to be generated for statistical modeling;
# prob - a percentile threshold for identification of the true open reading frame(-s) in the empirical transcript.
# This function returns an object of class data.frame with four fields:
# "Transcript" - name of transcript;
# "ORFStart" - the first coordinate of the identified open reading frame;
# "ORFStop" - the last coordinate of the identified open reading frame;
# "ORFLength" - length of the identified open reading frame.
ORFindeR = function(x, typeORF, n, prob){
	require(Biostrings)
	require("seqinr")
	ORFs = list()
	for (i in 1:length(x)){
		if (typeORF == "ATG"){
			ORF = ORFFindingATG(as.character(x[[i]]))
		}
		if (typeORF == "CTG"){
			ORF = ORFFindingCTG(as.character(x[[i]]))
		}
		if (typeORF == "GTG"){
			ORF = ORFFindingGTG(as.character(x[[i]]))
		}
		if (typeORF == "TTG"){
			ORF = ORFFindingTTG(as.character(x[[i]]))
		}
		if (length(ORF[[1]]) == 0){
			ORF = data.frame(cbind(names(x[i]), 0, 0, 0), stringsAsFactors = FALSE)
			colnames(ORF) = c("Transcript", "ORFStart", "ORFStop", "ORFLength")
			ORFs[[i]] = ORF
		}
		else if (length(ORF[[1]]) >= 1){
			ORF = data.frame(cbind(ORF[[1]], ORF[[2]], ORF[[3]]))
			RANDSeq = generateSeqsWithMultinomialModel(as.character(x[[i]]), n)
			RANDSeqORFlengths = numeric()
			for (k in 1:length(RANDSeq)){
				RANDSeq1 = RANDSeq[k]
				RANDSeq1 = ORFFindingATG(RANDSeq1)
				lengths = RANDSeq1[[3]]
				RANDSeqORFlengths = append(RANDSeqORFlengths, lengths, after = length(RANDSeqORFlengths))
			}
			RAND = subset(RANDSeqORFlengths, RANDSeqORFlengths > 0)
			RAND = quantile(RAND, probs = prob)
			ORF = ORF[!ORF[, 3] < RAND, ]
			if(length(ORF$X1) == 0){
				ORF[1,] = 0
				ORF = cbind(names(x[i]), ORF)
				colnames(ORF) = c("Transcript", "ORFStart", "ORFStop", "ORFLength")
				ORFs[[i]] = ORF
			}
			if(length(ORF$X1) >= 1){
				ORF = cbind(names(x[i]), ORF)
				colnames(ORF) = c("Transcript", "ORFStart", "ORFStop", "ORFLength")
				ORFs[[i]] = ORF
			}
			if(i%%1 == 0){
				cat("Transcript number", i, "of", length(x), "has been processed...\n")
				flush.console()
			}
		}
	}
	ORFs = do.call(rbind.data.frame, ORFs)
	return(ORFs)
}
