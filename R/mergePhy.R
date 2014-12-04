##' merge two phylip files into one. 
##'
##' 
##' @title mergePhy
##' @param p1 phylip file 1
##' @param p2 phylip file 2
##' @param outfile output file
##' @return NULL
##' @export
##' @author ygc
mergePhy <- function(p1, p2, outfile) {
    p1.info <- getPhyInfo(p1)
    p2.info <- getPhyInfo(p2)
    
    if ( ! p1.info$width == p2.info$width ) {
        stop("sequence width not consistent in two phylip files...")
    }
    
    pp1 <- readLines(p1)
    pp2 <- readLines(p2)
    if ( length(pp1) != (p1.info$num +1) || length(pp2) != (p2.info$num + 1) ) {
        stop("only collapse version (sequence in one line) of phylip is supported...")
    }

    num <- p1.info$num + p2.info$num
    hh <- paste(num, "\t", p1.info$width)
    out <- file(outfile, "w")
    writeLines(hh, out)
    for (i in 2:length(pp1)) {
        writeLines(pp1[i], out)
    }
    for (i in 2:length(pp2)) {
        writeLines(pp2[i], out)
    }
    close(out)
}

getPhyInfo <- function(phy) {
    line1 <- readLines(phy, n=1)
    res <- strsplit(line1, split="\\s")[[1]]
    res <- res[res != ""]

    return(list(num=as.numeric(res[1]), width=as.numeric(res[2])))
}
