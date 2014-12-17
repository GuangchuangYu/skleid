##' read aligned sequences in phylip format
##'
##' 
##' @title read.phylip
##' @param file phylip file
##' @importFrom Biostrings BStringSet
##' @return BStringSet object
##' @export
##' @author ygc
read.phylip <- function(file) {
    info <- getPhyInfo(file)
    phy <- readLines(file)

    type <- "standard"
    if (length(phy) == (info$num + 1)) {
        type <- "oneLine"
    }
    if (type == "standard") {
        stop("not implemented yet...")
    }
    if (type == "oneLine") {
        seq <- lapply(2:length(phy), function(i) unlist(strsplit(phy[i], "\\s+")))
        seq2 <- sapply(seq, function(x) x[2])
        names(seq2) <- sapply(seq, function(x) x[1])
        res <- BStringSet(seq2)
    }
    return(res)
}

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

##' convert fasta to phylip format
##'
##' 
##' @title fas2phy 
##' @param fas aligned sequences in fasta format
##' @param type one of DNA, RNA or AA
##' @param outfile output file
##' @param collapse collapse to 1 line or not
##' @return NULL
##' @importFrom Biostrings readDNAMultipleAlignment
##' @importFrom Biostrings readRNAMultipleAlignment
##' @importFrom Biostrings readAAMultipleAlignment
##' @importFrom Biostrings write.phylip
##' @importFrom Biostrings width
##' @export
##' @author ygc
fas2phy <- function(fas, type="DNA", outfile="out.phy", collapse=FALSE) {
    if (type == "DNA") {
        x <- readDNAMultipleAlignment(fas)
    } else if (type == "RNA") {
        x <- readRNAMultipleAlignment(fas)
    } else if (type == "AA") {
        x <- readAAMultipleAlignment(fas)
    } else {
        stop("type should be one of DNA, RNA or AA...\n")
    }
    
    if (collapse == TRUE) {
        out <- file(outfile, "w")
        dna <- x@unmasked
        header <- paste(length(dna), "\t", width(dna[1]))
        writeLines(header, out)
        nn <- max(nchar(names(dna)))
        for (i in 1:length(dna)) {
            n <- names(dna[i])
            sep.blank <- paste(rep(" ", nn-nchar(n)+4), sep="", collapse="")
            line <- paste(n, sep.blank, toString(dna[i]), sep="")
            writeLines(line, out)
        }
        close(out)
    } else {
        write.phylip(x, outfile)
    }
}

getPhyInfo <- function(phy) {
    line1 <- readLines(phy, n=1)
    res <- strsplit(line1, split="\\s")[[1]]
    res <- res[res != ""]

    return(list(num=as.numeric(res[1]), width=as.numeric(res[2])))
}
