##' add keyword in description line of FASTA file
##'
##' 
##' @title addKeywordFASTA
##' @param fasta fasta file
##' @param keyword keyword to be added
##' @param suffix logical, true add suffix, false add prefix
##' @importFrom Biostrings readBStringSet
##' @return NULL
##' @export
##' @author ygc
addKeywordFASTA <- function(fasta, keyword, suffix=TRUE) {
    ff <- readBStringSet(fasta)
    
    if (suffix == TRUE) {
        names(ff) <- paste(names(ff), keyword, sep="")
    } else {
        names(ff) <- paste(keyword, names(ff), sep="")
    }
    writeXStringSet(ff, fasta)
}
