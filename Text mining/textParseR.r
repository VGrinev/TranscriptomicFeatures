###################################################################################################
##  A simple R code for some parsing of text.                                                    ##
##  (c) GNU GPL Vasily V. Grinev, 2019. grinev_vv[at]bsu.by                                      ##
###################################################################################################
##  Installation of pdftools and stringr libraries from CRAN.
setRepositories()
install.packages("pdftools")
install.packages("stringr")
##  Loading of required auxiliary libraries.
suppressMessages(require(pdftools))
suppressMessages(require(stringr))
##  Setting of a work directory.
setwd("D:/Vasily Grinev")
##  Loading of a text from PDF file of interest.
#   This code was successfully tested with libraries pdftools v.2.2 and stringr v.1.4.0.
text.original = pdf_text(pdf = "Grinev V. V. et al. Manuscript_R1.pdf")
##  Inspection of text structure.
str(text.original)
text.original
##  Collection of pages to be parsed.
text.original = text.original[2:14]
##  Processing each page of text as a separate object.
words.all = list()
for (i in 1:length(x = text.original)){
     page = strsplit(x = text.original, split = "\r\n")[[i]]
     page = page[c(-1, -length(x = page))]
     textColumnL = list()
     for (j in 1:length(x = page)){
          textColumnL[[j]] = str_sub(string = page[j], start = 1L, end = 75L)
     }
     textColumnR = list()
     for (k in 1:length(x = page)){
          textColumnR[[k]] = str_sub(string = page[k], start = 76L, end = 140L)
     }
     text.primary = rbind(do.call(what = rbind, args = textColumnL),
                          do.call(what = rbind, args = textColumnR))
     text.primary = str_replace_all(string = text.primary,
                                    pattern = fixed("  "),
                                    replacement = "")
     text.primary = str_trim(string = text.primary, side = "both")
     text.primary = text.primary[text.primary != ""]
     text.primary = paste(text.primary, collapse = " ")
     text.primary = gsub(pattern = "- ", replacement = "", x = text.primary)
     text.primary = gsub(pattern = "[.]", replacement = "", x = text.primary)
     text.primary = gsub(pattern = ",", replacement = "", x = text.primary)
     text.primary = gsub(pattern = "[(]", replacement = "", x = text.primary)
     text.primary = gsub(pattern = ")", replacement = "", x = text.primary)
     text.primary = strsplit(x = text.primary, split = " ")[[1]]
     words.all[[i]] = text.primary
}
words.all = do.call(what = c, args = words.all)
##  Loading of abbreviations to be checked as a character vector.
abrevs = read.table(file = "Abbreviations.txt", sep = "\t", quote = "\"", as.is = TRUE)$V1
##  Some replacements if necessary.
abrevs[1] = sort(x = words.all)[662]
abrevs[2] = sort(x = words.all)[82]
abrevs[3] = sort(x = words.all)[113]
abrevs[4] = sort(x = words.all)[115]
##  Collection of statistics.
words.all = table(words.all)
words.targets = words.all[attr(words.all, "dimnames")$words.all %in% abrevs]
freq = list(number = length(x = unique(x = abrevs)),
            frequency = sum(words.targets)/sum(words.all))
##  Inspection of statistics.
freq
