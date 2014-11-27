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
        file(outfile, "w")
        dna <- x@unmasked
        header <- paste(length(dna), "\t", width(dna))
        writeLines(header, outfile)
        nn <- max(nchar(names(dna)))
        for (i in 1:length(dna)) {
            n <- names(dna[i])
            sep.blank <- paste(rep(" ", nn-nchar(n)+4), sep="", collapse="")
            line <- paste(n, sep.blank, toString(dna[i]), sep="")
            writeLines(line, outfile)
        }
        close(outfile)
    } else {
        write.phylip(x, outfile)
    }
}
