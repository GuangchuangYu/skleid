##' remove sequence(s) in alignment by the specific index.
##'
##' 
##' @title removeSeq
##' @param aln alignment
##' @param idx index
##' @return alignment 
##' @author ygc
##' @export
removeSeq <- function(aln, idx) {
    aln$seqs <- aln$seqs[-idx,]
    aln$num <- aln$num - length(idx)
    return(aln)
}
