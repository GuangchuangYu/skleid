##' filter out read numbers below specific percentage
##'
##' 
##' @title readNumFilter
##' @param files files
##' @param percentage percentage 
##' @return NULL
##' @author ygc
##' @export
readNumFilter <- function(files, percentage=2) {
    ## remove contig with criteria of reads number less than 2% (by default) 
    for ( ff in files ) {
        readNumFilter_internal(ff, percentage)
    }
}

##' @importFrom Biostrings readDNAStringSet
readNumFilter_internal <- function(file, percentage=2) {
    y <- readDNAStringSet(file)
    if (length(grep("numreads=", names(y)))) {
        ## num of reads
        nr <- sub(".*numreads=", "", names(y))
        nr <- as.numeric(nr)
        i <- nr/sum(nr) * 100 > percentage
        if (sum(i) >= 1) {
            writeXStringSet(y[i], filepath=file)
        }
        cat(file, "\t", "processed...\n")
    } else {
        cat(file, "\t", "omitted...\n")
    }
}
