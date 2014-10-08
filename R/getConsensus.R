##' determine the consensus sequence
##'
##' 
##' @title getConsensus
##' @param seq.df sequence data.frame, sequence is stored in row
##' @return consensus sequence in character vector
##' @author ygc
##' @export
getConsensus <- function(seq.df) {
    cs <- apply(seq.df, 2, function(i) {
        nt <- unique(i)
        nt <- nt[nt!="-"]
        nt <- toupper(nt)
        nt <- sort(nt)
        if (length(nt) == 1) {
            return(nt)
        } else if (length(nt) == 2) {
            if ( all(nt == c("A", "C")) ) {
                return("M")
            } else if (all(nt == c("A", "G"))) {
                return("R")
            } else if (all(nt == c("A", "T"))) {
                return("W")
            } else if (all(nt == c("C", "G"))) {
                return("S")
            } else if (all(nt == c("C", "T"))) {
                return("Y")
            } else if (all(nt == c("G", "T"))) {
                return("K")
            }
        } else if ( length(nt) == 3) {
            if (all(nt == c("A", "C", "G"))) {
                return("V")
            } else if (all(nt== c("A", "C", "T"))) {
                return("H")
            } else if (all(nt == c("A", "G", "T"))) {
                return("D")
            } else if (all(nt == c("C", "G", "T"))) {
                return("B")
            }
        } else if (length(nt) == 4) {
            ## length of 4
            return("N")
        } else {
            ## length of 0, only -
            return("-")
        }
            
    })
    cs <- unlist(cs)
    cs <- cs[cs != "-"]
    return(cs)
}
