###################################################################################################
##  R wrap for converting GTF/GFF annotations to a local SQLite database.                        ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
### Arguments of function:
##  d.work   - path to and name of work directory.
##  f.gtf    - the name of GTF/GFF file with reference annotations.
##  f.sqlite - the name of the SQLite database to be stored locally.
GTFtoSQLite = function(d.work, f.gtf, f.sqlite){
##  Loading of required auxiliary libraries.
#   This code was successfully tested with libraries GenomicFeatures v.1.36.4
#   and AnnotationDbi v.1.46.0.
suppressMessages(require(GenomicFeatures))
##  Make a TxDb object from annotations available as a GTF or GFF3 file.
TxDb = makeTxDbFromGFF(file = paste(d.work, f.gtf, sep = "/"))
##  Save the database to the specified file.
saveDb(x = TxDb, file = paste(d.work, f.sqlite, sep = "/"))
}
