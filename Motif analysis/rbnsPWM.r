###################################################################################################
##  A complete R-wrap for calculation of splicing factors positional weight matrices             ##
##  using RNA Bind-n-Seq data.                                                                   ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### It is an R implementation of an algorithm developed by Daniel Dominguez with colleagues and
#   published in Molecular Cell (Dominguez et al., 2018, Molecular Cell 70, 854–867).
### Arguments of function:
##  d.work  - character string giving the path to and name of work directory.
##  f.input - character string giving the name of TXT file in tab-delimited format containing
#             pre-ranked kmers. This file should include two mandatory fields:
#             i) olig        - actual sequences of kmers;
#             ii) enrichment - enrichment scores.
##  kmer    - integer giving the length of oligomers (value of k). The current experimental
#             algorithm works only with k = 5. So, default value of kmer is 5.
##  bg.freq - background frequency of nucleotides. By default, the frequency is uniform.
rbnsPWM = function(d.work, f.input, kmer = 5, bg.freq = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)){
##  Loading of pre-ranked oligomers.
input = read.table(file = paste(d.work, f.input, sep = "/"),
                   sep = "\t",
                   header = TRUE,
                   quote = "\"",
                   as.is = TRUE)
#   Descending order.
input = input[order(-input[, 2]), ]
##  Sequential alignment.
list.mers = strsplit(x = input[, 1], split = "")
#   Development a list of founders.
pr.fndrs = list()
pr.fndrs[[1]] = list.mers[[1]]
for (i in 2:length(list.mers)){
    if (table(pr.fndrs[[1]] == list.mers[[i]])["FALSE"][[1]] == 1){
    }else if (table((c("", pr.fndrs[[1]]) ==
                     c(list.mers[[i]], ""))[2:5])["FALSE"][[1]] == 1 ||
              is.na(table((c("", pr.fndrs[[1]]) ==
                           c(list.mers[[i]], ""))[2:5])["FALSE"][[1]])){
    }else if (table((c(pr.fndrs[[1]], "") ==
                     c("", list.mers[[i]]))[2:5])["FALSE"][[1]] == 1 ||
              is.na(table((c(pr.fndrs[[1]], "") ==
                           c("", list.mers[[i]]))[2:5])["FALSE"][[1]])){
    }else if (length(unique((c("", "", pr.fndrs[[1]]) ==
                             c(list.mers[[i]], "", ""))[3:5])) == 1 &&
              unique((c("", "", pr.fndrs[[1]]) ==
                      c(list.mers[[i]], "", ""))[3:5]) == "TRUE"){
    }else if (length(unique((c(pr.fndrs[[1]], "", "") ==
                             c("", "", list.mers[[i]]))[3:5])) == 1 &&
              unique((c(pr.fndrs[[1]], "", "") ==
                      c("", "", list.mers[[i]]))[3:5]) == "TRUE"){
    }else{
         pr.fndrs[[length(pr.fndrs) + 1]] = list.mers[[i]]
    }
}
if (length(pr.fndrs) <= 2){
    fndrs = pr.fndrs
}else{
    fndrs = list()
    n = 1
    fndrs[[n]] = pr.fndrs[[1]]
    pr.fndrs = pr.fndrs[-1]
    while(length(pr.fndrs) > 0){
          index.mer = double()
          if (length(pr.fndrs) > 1){
              for (j in 2:length(pr.fndrs)){
                   if (table(pr.fndrs[[1]] == pr.fndrs[[j]])["FALSE"][[1]] == 1){
                       index.mer[[j]] = j
                  }else if (table((c("", pr.fndrs[[1]]) ==
                                   c(pr.fndrs[[j]], ""))[2:5])["FALSE"][[1]] == 1 ||
                            is.na(table((c("", pr.fndrs[[1]]) ==
                                         c(pr.fndrs[[j]], ""))[2:5])["FALSE"][[1]])){
                       index.mer[[j]] = j
                  }else if (table((c(pr.fndrs[[1]], "") ==
                                   c("", pr.fndrs[[j]]))[2:5])["FALSE"][[1]] == 1 ||
                            is.na(table((c(pr.fndrs[[1]], "") ==
                                         c("", pr.fndrs[[j]]))[2:5])["FALSE"][[1]])){
                       index.mer[[j]] = j
                  }else if (length(unique((c("", "", pr.fndrs[[1]]) ==
                                           c(pr.fndrs[[j]], "", ""))[3:5])) == 1 &&
                                   unique((c("", "", pr.fndrs[[1]]) ==
                                           c(pr.fndrs[[j]], "", ""))[3:5]) == "TRUE"){
                       index.mer[[j]] = j
                  }else if (length(unique((c(pr.fndrs[[1]], "", "") ==
                                           c("", "", pr.fndrs[[j]]))[3:5])) == 1 &&
                                   unique((c(pr.fndrs[[1]], "", "") ==
                                           c("", "", pr.fndrs[[j]]))[3:5]) == "TRUE"){
                       index.mer[[j]] = j
                  }
             }
              index.mer = index.mer[!is.na(index.mer)]
              if (length(index.mer) > 0){
                  pr.fndrs = pr.fndrs[-index.mer]
             }
          n = n + 1
          fndrs[[n]] = pr.fndrs[[1]]
          pr.fndrs = pr.fndrs[-1]
         }else if (length(pr.fndrs) == 1){
                   n = n + 1
                   fndrs[[n]] = pr.fndrs[[1]]
                   pr.fndrs = pr.fndrs[-1]
         }else if (length(pr.fndrs) == 0){
                   break
         }
    }
}
#   Sequential alignment.
LIST.ALIGNS = list()
sub.set = list.mers[-c(match(fndrs, list.mers))]
for (k in 1:length(fndrs)){
     list.align = list()
     for (l in 1:length(sub.set)){
          if (table(fndrs[[k]] == sub.set[[l]])["FALSE"][[1]] == 1){
              list.align[[l]] = rbind(fndrs[[k]], sub.set[[l]])
         }else if (table((c("", fndrs[[k]]) ==
                          c(sub.set[[l]], ""))[2:5])["FALSE"][[1]] == 1 ||
                   is.na(table((c("", fndrs[[k]]) ==
                                c(sub.set[[l]], ""))[2:5])["FALSE"][[1]])){
              list.align[[l]] = rbind(c("", fndrs[[k]]), c(sub.set[[l]], ""))
         }else if (table((c(fndrs[[k]], "") ==
                          c("", sub.set[[l]]))[2:5])["FALSE"][[1]] == 1 ||
                   is.na(table((c(fndrs[[k]], "") ==
                                c("", sub.set[[l]]))[2:5])["FALSE"][[1]])){
                   list.align[[l]] = rbind(c(fndrs[[k]], ""), c("", sub.set[[l]]))
         }else if (length(unique((c("", "", fndrs[[k]]) ==
                                  c(sub.set[[l]], "", ""))[3:5])) == 1 &&
                   unique((c("", "", fndrs[[k]]) ==
                           c(sub.set[[l]], "", ""))[3:5]) == "TRUE"){
                   list.align[[l]] = rbind(c("", "", fndrs[[k]]), c(sub.set[[l]], "", ""))
         }else if (length(unique((c(fndrs[[k]], "", "") ==
                                  c("", "", sub.set[[l]]))[3:5])) == 1 &&
                   unique((c(fndrs[[k]], "", "") ==
                           c("", "", sub.set[[l]]))[3:5]) == "TRUE"){
                   list.align[[l]] = rbind(c(fndrs[[k]], "", ""), c("", "", sub.set[[l]]))
         }
    }
     if (length(list.align) == 0){
         LIST.ALIGNS[[k]] = fndrs[[k]]
    }else{
     LIST.ALIGNS[[k]] = list.align
    }
}
ALIGNS = list()
for (m in 1:length(LIST.ALIGNS)){
     align = LIST.ALIGNS[[m]]
     align = align[lapply(X = align, FUN = length) > 0]
     if (class(align) == "character"){
         full.mat = align
    }else{
          if(length(align) == 1){
             full.mat = align[[1]]
         }
          if(length(align) > 1){
             full.mat = matrix(ncol = 9, nrow = length(align) + 1, "")
             full.mat[1, 3:7] = align[[1]][1, ][align[[1]][1, ] != ""]
             for (n in 1:length(align)){
                  if (ncol(align[[n]]) == 5){
                      full.mat[n + 1, 3:7] = align[[n]][2, ]
                 }
                  if (ncol(align[[n]]) == 6 & align[[n]][1, 1] == ""){
                      full.mat[n + 1, 2:6] = align[[n]][2, ][align[[n]][2, ] != ""]
                 }
                  if (ncol(align[[n]]) == 6 & align[[n]][1, ncol(align[[n]])] == ""){
                      full.mat[n + 1, 4:8] = align[[n]][2, ][align[[n]][2, ] != ""]
                 }
                  if (ncol(align[[n]]) == 7 & align[[n]][1, 1] == ""){
                      full.mat[n + 1, 1:5] = align[[n]][2, ][align[[n]][2, ] != ""]
                 }
                  if (ncol(align[[n]]) == 7 & align[[n]][1, ncol(align[[n]])] == ""){
                      full.mat[n + 1, 5:9] = align[[n]][2, ][align[[n]][2, ] != ""]
                 }
            }
         }
    }
     if (class(full.mat) == "character"){
         ALIGNS[[m]] = full.mat
    }else{
         ALIGNS[[m]] = full.mat[, colSums(full.mat != "") != 0]
    }
}
##  Assigning weights to each element of the matrix.
ALIGNS = lapply(X = ALIGNS,
                FUN = function(x){if (class(x) == "character"){
                                      c(x, input[paste(x, collapse = "") == input[, 1], 2])
                                 }else{
                                       cbind(x,
                                             apply(X = x,
                                                   MARGIN = 1,
                                                   FUN = function(y){input[paste(y,
                                                                                 collapse = "") ==
                                                                           input[, 1], 2]}))}})
ALIGNS = lapply(X = ALIGNS, FUN = function(x){if (class(x) == "matrix"){
                                                  cbind(x, x[1, ncol(x)])
                                             }else{
                                                  return(x)
                                             }})
mer = lapply(X = ALIGNS,
             FUN = function(x){if (class(x) == "character"){
                                   c(paste(x[1:(length(x) - 1)], collapse = ""), x[length(x)])
                              }else{
                                    t(apply(X = x,
                                            MARGIN = 1,
                                            FUN = function(y){y = cbind(paste(y[1:(length(y) - 2)],
                                                                              collapse = ""),
                                                                        y[length(y)])}))}})
mer = do.call(what = rbind, args = mer)
mer = aggregate(x = as.numeric(mer[, 2]), by = list(mer[, 1]), FUN = sum)
ALIGNS = lapply(X = ALIGNS,
                FUN = function(x){if (class(x) == "character"){
                                      return(x)
                                 }else{
                                      cbind(x[, -ncol(x)],
                                            apply(X = x,
                                                  MARGIN = 1,
                                                  FUN = function(y){y[length(y)] =
                                                                    mer[paste(y[1:(length(y) - 2)],
                                                                              collapse = "") ==
                                                                    mer[, 1], 2]}))}})
ALIGNS = lapply(X = ALIGNS,
                FUN = function(x){if (class(x) == "character"){
                                      return(x)
                                 }else{
                                      x[, ncol(x) - 1] =
                                      round(as.numeric(x[, ncol(x) - 1]) *
                                            as.numeric(x[1, ncol(x)])/as.numeric(x[, ncol(x)]),
                                            digits = 3)
                                      x = x[, -ncol(x)]
                                      return(x)}})
Pr.alignments = ALIGNS
#   Trimming of leading and/or lagging positions having more than 75% offset weight.
for (o in 1:length(ALIGNS)){
     if (class(ALIGNS[[o]]) == "matrix"){
         offset.weight = double()
         for (p in 1:(ncol(ALIGNS[[o]]) - 1)){
              offset.weight[[p]] = sum(as.numeric(ALIGNS[[o]][ALIGNS[[o]][, p] ==
                                       "", ncol(ALIGNS[[o]])]))/
                                   sum(as.numeric(ALIGNS[[o]][, ncol(ALIGNS[[o]])]))
        }
         offset.weight = c(offset.weight > 0.75, "FALSE")
         ALIGNS[[o]] = cbind(ALIGNS[[o]][, offset.weight == "FALSE"])
    }
}
#  Calculation of positional count matrices. Here the enrichment scores are used as counts.
PCMs = list()
for (s in 1:length(ALIGNS)){
     if (class(ALIGNS[[s]]) == "character"){
         mat = matrix(data = 0L, ncol = length(ALIGNS[[s]]) - 1, nrow = 4)
         rownames(mat) = c("A", "C", "G", "T")
         for (u in 1:4){
              mat[u, rownames(mat)[u] == ALIGNS[[s]][1:(length(ALIGNS[[s]]) - 1)]] =
                                         as.numeric(ALIGNS[[s]][length(ALIGNS[[s]])])
        }
    }else{
         mat = matrix(ncol = ncol(ALIGNS[[s]]) - 1, nrow = 4)
         rownames(mat) = c("A", "C", "G", "T")
         for (u in 1:ncol(mat)){
              mat[1, u] = sum(as.numeric(ALIGNS[[s]][ALIGNS[[s]][, u] == "A", ncol(ALIGNS[[s]])]))
              mat[2, u] = sum(as.numeric(ALIGNS[[s]][ALIGNS[[s]][, u] == "C", ncol(ALIGNS[[s]])]))
              mat[3, u] = sum(as.numeric(ALIGNS[[s]][ALIGNS[[s]][, u] == "G", ncol(ALIGNS[[s]])]))
              mat[4, u] = sum(as.numeric(ALIGNS[[s]][ALIGNS[[s]][, u] == "T", ncol(ALIGNS[[s]])]))
        }
    }
     PCMs[[s]] = mat
}
#  Calculation of pseudocount adjusted positional probability matrices.
PPMs = lapply(X = PCMs, FUN = function(y){y[y == 0] = min(y[y > 0]) * 0.001
                                          y = round(x = t(t(y)/apply(X = y,
                                                                     MARGIN = 2,
                                                                     FUN = sum)),
                                                    digits = 10)})
#   Calculation of positional weight matrices.
PWMs = lapply(X = PPMs, FUN = function(x, y){log2(x/(y = bg.freq))})
#   Calculation of information content profiles.
ICPs = lapply(X = PCMs,
              FUN = function(x){x = t(t(x)/apply(X = x, MARGIN = 2, FUN = sum))
                                2 + apply(X = x, MARGIN = 2, FUN = function(y){y = y[y > 0]
                                                                               sum(y * log2(y))})})
#   Calculation of information content matrices.
ICMs = list()
for (u in 1:length(PCMs)){
     x = t(t(PCMs[[u]])/apply(X = PCMs[[u]], MARGIN = 2, FUN = sum))
     ICMs[[u]] = t(ICPs[[u]] * t(x))
}
##  Returning the final object.
return(list("Primary alignments" = Pr.alignments,
            PCMs = PCMs,
            PPMs = PPMs,
            PWMs = PWMs,
            ICPs = ICPs,
            ICMs = ICMs))
}
