##' write consensus sequence to fasta file
##'
##' 
##' @title writeConsensus
##' @param aln alignment
##' @param output output fasta file 
##' @return NULL
##' @author ygc
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings writeXStringSet
##' @export
writeConsensus <- function(aln, output) {
    aln$seqs <- aln$seqs[-1,]
    seq.df <- aln2seqDF(aln)
    idx <- apply(seq.df, 2, function(i) all(i=="-"))
    seq.df <- seq.df[,!idx]

    cs <- getConsensus(seq.df)
    cs2 <- DNAStringSet(paste(cs, collapse=""))
    names(cs2) <- output
    writeXStringSet(cs2, output)
}



