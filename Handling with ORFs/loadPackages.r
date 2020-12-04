### Loading of required packages.
##  Listing of required packages.
reqPkgs <- c("Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "data.table",
             "parallel", "Peptides", "Rcpp", "randomForest", "rtracklayer",
             "stringr")
##  Checking for installed required packages.
insPkgs <- reqPkgs[reqPkgs %in% installed.packages()[, 1]]
if (length(x=reqPkgs) > length(x=insPkgs)){
    if (length(x=reqPkgs) - length(x=insPkgs) == 1){
        stop(paste(paste("The following package has been not installed: ",
                         toString(x=paste0(reqPkgs[!reqPkgs %in% insPkgs])),
                         ".",
                         sep=""),
                   "Please, install this package before analysis."))
    }else{
        stop(paste(paste("The following packages have been not installed: ",
                         toString(x=paste0(reqPkgs[!reqPkgs %in% insPkgs])),
                         ".",
                         sep=""),
                   "Please, install these packages before analysis."))
    }
}else{
    print("All required packages are installed.")
}
##  Loading of required packages.
suppressMessages(expr=library(package=Biostrings))
suppressMessages(expr=library(package=BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(expr=library(package=data.table))
suppressMessages(expr=library(package=parallel))
suppressMessages(expr=library(package=Peptides))
suppressMessages(expr=library(package=Rcpp))
suppressMessages(expr=library(package=randomForest))
suppressMessages(expr=library(package=rtracklayer))
suppressMessages(expr=library(package=stringr))
##  Setting of C++ environment.
setwd("your_work_directory")
Rcpp::sourceCpp("getBaoMetrics.cpp")
Rcpp::sourceCpp("getCorrelationFactors.cpp")
