##' reverse complement of DNA string
##'
##' 
##' @title revcomp
##' @param seq DNA sequence string
##' @return character
##' @author ygc
##' @export
revcomp <- function(seq) {
	## complement
	seq2 <- chartr("ATCGMKRYWSVBHDXN", "TAGCKMYRWSBVDHXN", seq)
	## vectorize
	seq2 <- substring(seq2, 1:nchar(seq2), 1:nchar(seq2))
	## reverse
	x <- paste(rev(seq2), collapse="")
	return(x)	
}
