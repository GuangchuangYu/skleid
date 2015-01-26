##' add keyword in description line of FASTA file
##'
##' 
##' @title addKeywordFASTA
##' @param fasta fasta file
##' @param keyword keyword to be added
##' @param by one of suffix, prefix or mapping
##' @param outfile output file
##' @importFrom Biostrings readBStringSet
##' @return NULL
##' @export
##' @author ygc
addKeywordFASTA <- function(fasta, keyword, by, outfile=NULL) {
    by <- match.arg(by, c("suffix", "prefix", "mapping"))
    
    ff <- readBStringSet(fasta)
    
    if ( by == "suffix") {
        names(ff) <- paste(names(ff), keyword, sep="")
    } else if (by == "prefix") {
        names(ff) <- paste(keyword, names(ff), sep="")
    } else {
        kw <- read.csv(keyword, stringsAsFactor=FALSE)
        kw <- unique(kw)
        
        nn <- names(ff)
        ## idx <- sapply(kw[,1], function(x) {
        ##     y <- grep(x, nn)
        ##     if (length(y) == 0) {
        ##         return(FALSE)
        ##     }
        ##     return(TRUE)
        ## })

        ## kw <- unique(kw[idx,])

        
        for (i in 1:nrow(kw)) {
            nn <- gsub(kw[i,1], paste(kw[i, 1], kw[i, 2], sep="_"), nn)
        }

        names(ff) <- nn
    }
    
    if (is.null(outfile)) {
        writeXStringSet(ff, fasta)
    } else {
        writeXStringSet(ff, outfile)
    }
    invisible(ff)
}
