##' download genbank (.gb) file by GI number
##'
##' 
##' @title download_genbank_gi 
##' @param gi gi number
##' @param db supported db, currently 'nuccore'
##' @param outfile output file, by default, gi_number.gb
##' @return NULL
##' @export
##' @author Guangchuang Yu
download_genbank_gi <- function(gi, db="nuccore", outfile=NULL) {
    if (db != "nuccore") {
        stop("currently, only nuccore is supported...")
    }
    url <- paste0("www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=",
                  db, '&dopt=genbank&sort=&val=', gi)

    if (is.null(outfile)) {
        outfile <- paste0(gi, ".gb")
    }
    
    download.file(url=url, destfile = outfile, method="curl")
}

##' download genbank (.gb) file by accession number
##'
##' 
##' @title download_genbank_acc
##' @param acc accession number
##' @param db supported db, currently 'nuccore'
##' @param outfile output file, by default, acc_number.gb
##' @return NULL
##' @export
##' @author Guangchuang Yu
download_genbank_acc <- function(acc, db="nuccore", outfile=NULL) {
    if (is.null(outfile)) {
        outfile <- paste0(acc, ".gb")
    }
    gi <- acc2gi(acc)
    download_genbank_gi(gi, db, outfile)
}

acc2gi <- function(acc) {
    url <- paste0("http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=", acc)
    x <- readLines(url)
    i <- grep('<Id>', x)
    if (length(i) == 0) {
        stop("GI number not found...")
    }
    id <- x[i]
    gi <- gsub('</*Id>', '', id)
    return(gi)
}
