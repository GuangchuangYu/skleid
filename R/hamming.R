##' @rdname hamming
##' @return numeric
##' @method hamming BStringSet
##' @export
##' @author ygc
hamming.BStringSet <- function(seq1, seq2, indel=FALSE, ...) {
    s1 <- toString(seq1)
    s2 <- toString(seq2)
    hamming(s1, s2, indel, ...)
}

##' @rdname hamming
##' @method hamming character
##' @export
##' @author ygc
hamming.character <- function(seq1, seq2, indel=FALSE, ...) {
    if (length(seq1) == 1 & length(seq2) == 1) {
        n1 <- nchar(seq1)
        n2 <- nchar(seq2)
        if (n1 != n2) {
            stop("length of two sequences should be consistent...")
        }
        seq1 <- substring(seq1, 1:n1, 1:n1)
        seq2 <- substring(seq2, 1:n2, 1:n2)
    }

    if ( length(seq1) > 1 & length(seq2) > 1) {
        class(seq1) <- "characterVector"
        hamming(seq1, seq2, indel, ...)
    } else {
        stop("seq1 and seq2 should be character of equal length...")
    }
}

##' @rdname hamming
##' @method hamming characterVector
##' @export
##' @author ygc
hamming.characterVector <- function(seq1, seq2, indel, ...) {
    ii <- which(seq1 != seq2)
    
    if (indel == FALSE) {
        message("--> indel will not count...")
        ii <- ii[seq1[ii] != '-' & seq2[ii] != '-']
    }
    
    return(length(ii))
}

