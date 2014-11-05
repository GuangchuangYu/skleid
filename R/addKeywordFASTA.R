##' add keyword in description line of FASTA file
##'
##' 
##' @title addKeywordFASTA
##' @param fasta fasta file
##' @param keyword keyword to be added
##' @param suffix logical, true add suffix, false add prefix
##' @param type one of AA or DNA
##' @importFrom Biostrings readDNAStringSet
##' @importFrom Biostrings readAAStringSet
##' @return NULL
##' @export
##' @author ygc
addKeywordFASTA <- function(fasta, keyword, suffix=TRUE, type="DNA") {
    if (type == "DNA") {
        ff <- readDNAStringSet(fasta)
    } else if (type == "AA") {
        ff <- readAAStringSet(fasta)
    }
    if (suffix == TRUE) {
        names(ff) <- paste(names(ff), keyword, sep="")
    } else {
        names(ff) <- paste(keyword, names(ff), sep="")
    }
    writeXStringSet(ff, fasta)
}
