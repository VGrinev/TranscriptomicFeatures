###################################################################################################
##  R-wrap for quantile-based determination of appropriate threshold for splicing factor motif.  ##
##  (c) GNU GPL Vasily V. Grinev, 2017-2019. grinev_vv[at]bsu.by                                 ##
###################################################################################################
### Arguments of function:
##  d.work  - character string giving the path to and name of work directory.
##  f.motif - character string giving the name of TXT file in tab-delimited format containing
#             positional weight matix.
##  n       - integer, number of motif sequences to be assessed. Default value is 1e5.
##  bg.freq - background frequency of nucleotides. By default, the frequency is uniform.
##  quant   - integer, threshold value for quantile distribution of motif weights.
#             Default value is 0.99.
pwmSFsThr = function(d.work, f.motif, n = 1e5, bg.freq = c(0.25, 0.25, 0.25, 0.25), quant = 0.99){
##  Loading of motif to be analysed as a standard object of class "matrix".
pwm = data.matrix(frame = data.frame(read.table(file = paste(d.work, f.motif, sep = "/"),
                                                sep = "\t",
                                                header = FALSE,
                                                quote = "\"",
                                                row.names = 1,
                                                as.is = TRUE)),
                  rownames.force = TRUE)
##  Calculation of appropriate motif threshold.
m.length = length(pwm[1, ])
thr = double(3)
for (i in 1:3){
     weight = double(n)
     for (j in 1:n){
          r.motif = sample(x = c("A", "C", "G", "T"),
                           size = m.length,
                           replace = TRUE,
                           prob = bg.freq)
          weight[j] = sum(diag(x = pwm[r.motif, ]))
          names(weight)[j] = paste(r.motif, collapse = "")
    }
    thr[i] = quantile(x = weight, probs = quant)[[1]]
}
##  Returning the final object.
return(mean(thr))
}
