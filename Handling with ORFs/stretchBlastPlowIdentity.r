#' Identify of low identity protein alignment stretches
#' @description Identification of the low identity protein alignment stretches.
#' @param x name of tab-delimited TXT file with parsed protein alignments.
#'     It is typically output of function parseBlastPalignments().
#' @param hits name of file in tab-delimited TXT format with consolidated
#'     blastp hits. It is typically output of function parseBlastPhits().
#' @param thr integer value for minimal length of amino acids sequence with
#'     perfect match interrupting low identity alignment stretches.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that mean the current working directory.
#' @return an object of class data frame with coordinates of low identity
#'     stretches.
#' @author Vasily V. Grinev
#' @examples
#' x="alignments.txt"
#' hits="hits.txt"
#' thr=8
#' workDir="D:/Vasily Grinev"
#' res <- stretchBlastPlowIdentity(x=x, hits=hits, thr=8, workDir=workDir)
#' @export
#' Last updated: December 25, 2022.

stretchBlastPlowIdentity <- function(x, hits, thr=8, workDir=NULL){
  ### Setting the work directory.
  if (is.null(x=workDir)){
      workDir <- getwd()
  }
  ### Loading of the blastp alignments as data frame object.
  aligns <- read.table(file=paste(workDir, x, sep = "/"),
                       sep="\t",
                       header=TRUE,
                       quote="\"",
                       as.is=TRUE)
  ### Loading of the blastp hits as data frame object.
  hits <- read.table(file=paste(workDir, hits, sep = "/"),
                     sep="\t",
                     header=TRUE,
                     quote="\"",
                     as.is=TRUE)
  ### Parsing of the blastp alignments.
  query <- aligns[seq(from=1, to=nrow(x=aligns), by=3), ]
  midline <- aligns[seq(from=3, to=nrow(x=aligns), by=3), 2]
  coords <- list()
  for (i in 1:nrow(x=query)){
      N_ext <- hits[hits$query == query[i, ]$ID, ]$N
      MM <- grep(pattern=" ", x=strsplit(x=midline[i], split="")[[1]])
      gaps <- grep(pattern="-", x=strsplit(x=query[i, 2], split="")[[1]])
      if (length(x=gaps) > 0){
          MM <- MM[!MM %in% gaps]
          if (length(x=MM) > 0){
              adj <- tapply(X=gaps,
                            INDEX=cumsum(x=c(1, diff(x=gaps)) - 1 >= 8),
                            FUN=range)
              adj <- data.frame(do.call(what=rbind, args=adj))
              adj$dif <- adj$X2 - adj$X1 + 1
              for (j in 1:nrow(x=adj)){
                  MM[MM > adj$X2[j]] <- MM[MM > adj$X2[j]] - adj$dif[j]
                  adj[, 1:2] <- adj[, 1:2] - adj$dif[j]
              }
          }
      }
      C_ext <- hits[hits$query == query[i, ]$ID, ]$C
      y <- c(c(if (N_ext > 0){
                   1:N_ext
               }
               ),
             c(if (length(x=MM) > 0){
                   N_ext + MM
               }
               ),
             c(if (C_ext > 0){
                  (hits[hits$query == query[i, ]$ID, ]$q.end + 1):
                   hits[hits$query == query[i, ]$ID, ]$length
               }
               )
             )
      if (!is.null(x=y)){
          res <- tapply(X=y,
                        INDEX=cumsum(x=c(1, diff(x=y)) - 1 >= 8),
                        FUN=range)
          res <- cbind(query[i, ]$ID, do.call(what=rbind, args=res))
      }else{
          res <- cbind(query[i, ]$ID, 0, 0)
      }
      colnames(x=res) <- c("query", "start", "end")
      rownames(x=res) <- NULL
      res <- data.frame(res)
      res$start <- as.numeric(x=res$start)
      res$end <- as.numeric(x=res$end)
      coords[[i]] <- res
  }
  coords <- data.frame(do.call(what=rbind, args=coords))
  coords <- coords[coords$end - coords$start + 1 >=thr, ]
  rownames(x=coords) <- NULL
  ### Returning of the final results.
  return(coords)
}