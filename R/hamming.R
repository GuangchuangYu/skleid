##' @rdname hamming
##' @return numeric
##' @method hamming BStringSet
##' @export
##' @author ygc
hamming.BStringSet <- function(seq1, seq2=NULL, indel=FALSE, ...) {
    seqs <- seq1
    if (!is.null(seq2)) {
        warning("seq2 will be omitted...")
    }
    
    n <- length(seqs) 
    if (n < 2) {
        stop("at least 2 sequences was expected...")
    }
    if ( n == 2 ) {
        seq1 <- seqs[1,]
        seq2 <- seqs[2,]
        s1 <- toString(seq1)
        s2 <- toString(seq2)
        res <- hamming(s1, s2, indel, ...)
    } else {
        res <- matrix(NA, nrow=n, ncol=n)
        colnames(res) <- rownames(res) <- names(seqs)
        seqs2 <- lapply(seqs, function(x) {
            y <- toString(x)
            substring(y, 1:nchar(y), 1:nchar(y))
        })
        cnt <- 1
        pb <- txtProgressBar(min=0, max=sum(1:n), style=3)
        for (i in 1:n) {
            for (j in 1:i) {
                setTxtProgressBar(pb, cnt)
                cnt <- cnt + 1
                if (i == j) {
                    res[i,j] <- 0
                } else {
                    res[i,j] <- suppressMessages(hamming.characterVector(seqs2[[i]], seqs2[[j]], indel, ...))
                    res[j, i] <- res[i, j]
                }
            }
        }
        close(pb)
    }
    return(res)
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

