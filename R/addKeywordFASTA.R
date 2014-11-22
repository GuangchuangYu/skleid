##' add keyword in description line of FASTA file
##'
##' 
##' @title addKeywordFASTA
##' @param fasta fasta file
##' @param keyword keyword to be added
##' @param suffix logical, true add suffix, false add prefix
##' @param outfile output file
##' @importFrom Biostrings readBStringSet
##' @return NULL
##' @export
##' @author ygc
addKeywordFASTA <- function(fasta, keyword, suffix=TRUE, outfile=NULL) {
    ff <- readBStringSet(fasta)
    
    if (suffix == TRUE) {
        names(ff) <- paste(names(ff), keyword, sep="")
    } else {
        names(ff) <- paste(keyword, names(ff), sep="")
    }
    if (is.null(outfile)) {
        writeXStringSet(ff, fasta)
    } else {
        writeXStringSet(ff, outfile)
    }
}
