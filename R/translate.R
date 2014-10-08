##' translate DNA sequence to protein sequence
##'
##' 
##' @title translateSeq
##' @param seq DNA sequence
##' @param start translate start position
##' @param print TRUE OR FALSE
##' @return character 
##' @author ygc
##' @importFrom Biostrings translate
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings toString
##' @export
translateSeq <- function(seq, start=1, print=TRUE) {
    s <- DNAStringSet(seq, start=start)
    p <- suppressWarnings(translate(s))
    res <- toString(p)
    if (print) 
        printFormatSeq(res)
    invisible(res)
}

##' auto determine the best translate start postion (offset 1,2 or 3),
##' and perform translation
##'
##' 
##' @title autoTranslate
##' @param seq DNA sequence 
##' @return list
##' @author ygc
##' @export
autoTranslate <- function(seq) {
    s1 <- suppressWarnings(translateSeq(seq, 1, print=FALSE))
    s2 <- suppressWarnings(translateSeq(seq, 2, print=FALSE))
    s3 <- suppressWarnings(translateSeq(seq, 3, print=FALSE))
    m1 <- max(nchar(unlist(strsplit(s1, "\\*"))))
    m2 <- max(nchar(unlist(strsplit(s2, "\\*"))))
    m3 <- max(nchar(unlist(strsplit(s3, "\\*"))))

    m <- m1
    s <- s1
    start <- 1
    if (m < m2) {
        s <- s2
        start <- 2
    }
    if (m < m3) {
        s <- s3
        start <- 3
    }

    peptides <- unlist(strsplit(s, "\\*"))
    nn <- nchar(peptides)
    i <- which.max(nn)

    AAstart <- regexpr(peptides[i], s)-1 +start -1 +regexpr("M", peptides[i]) -1
    AAstop <- regexpr(peptides[i], s) + nchar(peptides[i])

    startcodon <- AAstart*3-start+2
    stopcodon <- AAstop*3-(start-1)
    orf <- substring(seq, startcodon, stopcodon)
    pp <- peptides[i]
    pep <- substring(pp, AAstart+1)
    
    cat("start translate from position ", start, "\n")
    cat("longest peptide", "(", nchar(pep), 'AA)', "is\n", "\n")
    printFormatSeq(pep)
    
    cat("\nstart codon position is ", startcodon, "\n")
    cat("stop codon position is ", stopcodon, "\n")

    cat("\ncorresponding reading frame", "(", nchar(orf), 'NT)', "is\n", "\n")
    printFormatSeq(orf)
    
    result <- list(orf=orf, peptide=pep, start=startcodon, stop=stopcodon)
    invisible(result)
}
