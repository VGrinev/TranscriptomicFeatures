###################################################################################################
##  High-level R function for analysis of microhomology of sites                                 ##
##  surrounding exon-exon junctions.                                                             ##
##  (c) GNU GPL Vasily V. Grinev & Ilia M. Ilyushonak, 2019. grinev_vv[at]bsu.by                 ##
###################################################################################################
### Arguments of function:
##  eej.in - path to folder and name of the TXT file in tab-delimited format with
#   exon-exon junctions. This is input data file which should include four mandatory fields:
#   i) seqnames (name of chromosome or scaffold with prefix "chr");
#   ii) start (start genomic coordinate of the exon-exon junction);
#   iii) end (end genomic coordinate of the exon-exon junction);
#   iv) strand (strand information about exon-exon junction).
##  mer1 - integer value, maximum length of microhomology sites.
##  mer2 - integer value, minimum length of microhomology sites.
##  bs.genome - logical argument. The default value is TRUE, and BSgenome object will be used in
#   that case. Only GRCh38/hg38 assembly of H. sapiens genome is used in current implementation of
#   function.
##  file.fa - path to folder and name of the reference sequence in FASTA/FA format. This argument
#   is valid only with bs.genome = FALSE.
##  random - logical argument. This argument defins whether random sites surrounding exon-exon
#   junctions will be generated or not. The default value is FALSE.
##  exc - integer value, how many times random sites surrounding exon-exon junctions should excess
#   the number of empirical ones.
eejMicroHomology = function(eej.in, mer1, mer2, bs.genome = TRUE, file.fa, random = FALSE, exc){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries Rsamtools v.1.30.0 and stringdist v.0.9.5.1.
suppressMessages(require(Rsamtools))
suppressMessages(require(stringdist)
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))
##  Loading of exon-exon junctions.
eej = read.table(file = eej.in, sep = "\t", header = TRUE, quote = "\"", as.is = TRUE)
eej = eej[, which(colnames(eej)!= "width")]
##  Calculation genomic coordinates of up and down sites surrounding exon-exon junctions.
up = t(apply(eej, 1, function(x){if (x[5] == "+"){
                                     x[3] = as.numeric(x[3]) + 1
                                     x[4] = as.numeric(x[3]) + mer1 - 1
                                }else{
                                     x[3] = as.numeric(x[4]) - mer1
                                     x[4] = as.numeric(x[3]) + mer1 - 1
                                }
                                 return(x)
                                }
            )
       )
up = makeGRangesFromDataFrame(data.frame(up), keep.extra.columns = TRUE)
down = t(apply(eej, 1, function(x){if (x[5] == "+"){
                                       x[3] = as.numeric(x[4])
                                       x[4] = as.numeric(x[3]) + mer1 - 1
                                  }else{
                                       x[3] = as.numeric(x[3]) - mer1 + 1
                                       x[4] = as.numeric(x[3]) + mer1 - 1
                                  }
                                   return(x)
                                  }
              )
         )
down = makeGRangesFromDataFrame(data.frame(down), keep.extra.columns = TRUE)
##  Retrieving the sequences of up and down sites surrounding exon-exon junctions.
if (bs.genome == TRUE){
    seq.target = BSgenome.Hsapiens.UCSC.hg38
}else{
    seq.target = FaFile(file.fa)
}
up.list = list()
down.list = list()
suppressWarnings(for (i in 0:mer2){
                      up.temp = up
                      down.temp = down
                      end(up.temp[as.character(strand(up.temp)) == "+", ]) =
                          end(up.temp[as.character(strand(up.temp)) == "+", ]) - i
                      end(down.temp[as.character(strand(down.temp)) == "+", ]) =
                          end(down.temp[as.character(strand(down.temp)) == "+", ]) - i
                      start(up.temp[as.character(strand(up.temp)) == "-", ]) =
                            start(up.temp[as.character(strand(up.temp)) == "-", ]) + i
                      start(down.temp[as.character(strand(down.temp)) == "-", ]) =
                            start(down.temp[as.character(strand(down.temp)) == "-", ]) + i
                      up.list[1 + i] = getSeq(x = seq.target, up.temp)
                      down.list[1 + i] = getSeq(x = seq.target, down.temp)
                     }
                )
##  Pairwise alignment of up and down sites surrounding exon-exon junctions.
al.glob.list = list()
for (i in 1:(mer1 - mer2 + 1)){
     al.glob.list[[i]] = pairwiseAlignment(pattern = up.list[[i]],
                                           subject = down.list[[i]],
                                           type = "global")
}
names(al.glob.list) = c(mer1:mer2)
##  Calculation of percent identity for up and down sites surrounding exon-exon junctions.
score.pid = lapply(al.glob.list, function(x){pid(unlist(x), type = "PID1")})
names(score.pid) = c(mer1:mer2)
##  Calculation of Levenshtein distances for up and down sites surrounding exon-exon junctions.
score.lev = list()
for (i in 1:(mer1 - mer2 + 1)){
     score.lev[[i]] = stringdist(up.list[[i]], down.list[[i]], method = "lv")
}
names(score.lev) = c(mer1:mer2)
##  Consolidation of results.
mer = mer1
pid.tr = double()
for (i in 0:(mer1 - mer2)){
mer = mer1 - i
if (mer > 8){
    pid.tr[i + 1] = 1 - (table(score.pid[[i + 1]] >= 100 - 1/mer * 100)[[1]]/nrow(eej))
}else{
    pid.tr[i + 1] = 1 - (table(score.pid[[i + 1]] == 100)[[1]]/nrow(eej))
}
}
mer = mer1
lev.tr = double()
for (i in 0:(mer1 - mer2)){
mer = mer1 - i
if (mer > 8){
    lev.tr[i + 1] = 1 - (table(score.lev[[i + 1]] < 2)[[1]]/nrow(eej))
}else{
    lev.tr[i + 1] = 1 - (table(score.lev[[i + 1]] == 0)[[1]]/nrow(eej))
}
}
exper.sites = list(pid.tr = pid.tr, pid = score.pid, lev.tr = lev.tr, lev = score.lev)
if (random == TRUE){
##  Generation of random sequences.
up.rand = list()
down.rand = list()
for (i in c(1:(mer1 - mer2 + 1))){
up.rand.temp = strsplit(x = as.character(up.list[[i]]), split = "")
up.rand.temp = sample(x = up.rand.temp, size = nrow(eej) * exc, replace = TRUE)
up.rand[[i]] = DNAStringSet(x = do.call(what = rbind,
                                        args = lapply(X = up.rand.temp,
                                                      FUN = function(x){paste(x, collapse = "")})))
down.rand.temp = strsplit(x = as.character(down.list[[i]]), split = "")
down.rand.temp = sample(x = down.rand.temp, size = nrow(eej) * exc, replace = TRUE)
down.rand[[i]] = DNAStringSet(x = do.call(what = rbind,
                                          args = lapply(X = down.rand.temp,
                                                        FUN = function(x){paste(x, collapse = "")})))
}
##  Pairwise alignment of random sequences.
rand.al.glob.list = list()
for (i in 1:(mer1 - mer2 + 1)){
     rand.al.glob.list[[i]] = pairwiseAlignment(pattern = up.rand[[i]],
                                                subject = down.rand[[i]],
                                                type = "global")
}
names(rand.al.glob.list) = c(mer1:mer2)
##  Calculation of percent identity for random sequences.
score.pid.rand = lapply(rand.al.glob.list, function(x){pid(unlist(x), type = "PID1")})
names(score.pid.rand) = c(mer1:mer2)
##  Calculation of Levenshtein distances for random sequences.
score.lev.rand = list()
for (i in 1:(mer1 - mer2 + 1)){
     score.lev.rand[[i]] = stringdist(up.rand[[i]], down.rand[[i]], method = "lv")
}
names(score.lev.rand) = c(mer1:mer2)
##  Consolidation of results.
mer = mer1
pid.tr.rand = double()
for (i in 0:(mer1 - mer2)){
mer = mer1 - i
if (mer > 8){
    pid.tr.rand[i + 1] = 1 - (table(score.pid.rand[[i + 1]] >= 100 - 1/mer * 100)[[1]]/(nrow(eej) * exc))
}else{
    pid.tr.rand[i + 1] = 1 - (table(score.pid.rand[[i + 1]] == 100)[[1]]/(nrow(eej) * exc))
}
}
mer = mer1
lev.tr.rand = double()
for (i in 0:(mer1 - mer2)){
mer = mer1 - i
if (mer > 8){
    lev.tr.rand[i + 1] = 1 - (table(score.lev.rand[[i + 1]] < 2)[[1]]/(nrow(eej) * exc))
}else{
    lev.tr.rand[i + 1] = 1 - (table(score.lev.rand[[i + 1]] == 0)[[1]]/(nrow(eej) * exc))
}
}
rand.sites = list(pid.tr = pid.tr.rand, pid = score.pid.rand,
                  lev.tr = lev.tr.rand, lev = score.lev.rand)
return(list(exper.sites = exper.sites, rand.sites = rand.sites))
}else{
return(exper.sites)
}
}
### A simple example of function use.
#   source("http://bio.bsu.by/genetics/files/eejmicrohomology.r")
eej.in = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, NCBI RefSeq EEJs.txt"
#eej.in = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, all EEJs.txt"
eej.in = "D:/Vasily Grinev/RUNX1-RUNX1T1 project, RUNX1-RUNX1T1, canonical EEJs.txt"
mer1 = 12
mer2 = 6
file.fa = "D:/Vasily Grinev/hg38_der8.fa"
exc = 10
res = eejMicroHomology(eej.in = eej.in,
                       mer1 = mer1,
                       mer2 = mer2,
                       bs.genome = FALSE,
                       file.fa = file.fa,
                       random = TRUE,
                       exc = exc)
