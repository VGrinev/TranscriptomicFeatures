#' Parse of json output of blastp
#' @description Parsing of the json output of blastp.
#' @param x name of input file in JSON format.
#' @param hits name of file in tab-delimited TXT format with consolidated
#'     blastp hits. It is typically output of function parseBlastPhits().
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return an object of class data frame with two fields: ID and alignment.
#' @author Vasily V. Grinev
#' @examples
#' x="TSX982GZ114-Alignment.json"
#' hits="Kasumi-1 StringTie-based proteome, blastp best hits.txt"
#' workDir="D:/Vasily Grinev/AML proteome"
#' res <- parseBlastPalignments(x=x, hits=hits, workDir=workDir)
#' @export
#' Last updated: December 20, 2022.

parseBlastPalignments <- function(x, hits, workDir=NULL){
  ### Loading the required package.
  #   This code was successfully tested with package rjson v.0.2.21.
  suppressMessages(expr=library(package=rjson))
  ### Setting the work directory.
  if (is.null(x=workDir)){
      workDir <- getwd()
  }
  ### Loading of the json output of blastp as object of class list.
  json <- fromJSON(file=paste(workDir, x, sep = "/"))
  ### Loading of the blastp best hits as data frame object.
  h <- read.table(file=paste(workDir, hits, sep = "/"),
                  sep="\t",
                  header=TRUE,
                  quote="\"",
                  as.is=TRUE)
  ### Parsing of the json output.
  alignments <- list()
  for (i in 1:lengths(x=json)[[1]]){
#      bl2seq <- json$BlastOutput2[[i]]$report$results$bl2seq
      bl2seq <- json$BlastOutput2[[i]]$report$results
      if (!is.null(x=bl2seq)){
          if (sort(x=lengths(x=bl2seq))[1] == 5){
              bl2seq <- bl2seq[lengths(x=bl2seq) == 5]
              query <- bl2seq[[1]]$query_title
              subject <- h[h$query == query, ]$subject
              hsps <- list()
              for (j in 1:length(x=bl2seq)){
                  if (length(x=grep(pattern=subject,
                                    x=as.data.frame(x=bl2seq[[j]])[1, ])) > 0){
                      hsps[[j]] <- cbind(c(query, subject, "midline"),
                                         c(bl2seq[[j]]$hits[[1]]$hsps[[1]]$qseq,
                                           bl2seq[[j]]$hits[[1]]$hsps[[1]]$hseq,
                                           bl2seq[[j]]$hits[[1]]$hsps[[1]]$midline))
                      hsps <- hsps[lengths(x=hsps) > 0]
                      hsps <- do.call(what=rbind, args=hsps)
                 }
              }
          alignments[[i]] <- hsps
          }
      }
  }
  alignments <- do.call(what=rbind, args=alignments)
  colnames(x=alignments) <- c("ID", "alignment")
  ### Returning of the final results.
  return(alignments)
}
