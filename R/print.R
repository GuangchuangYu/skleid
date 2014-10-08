##' print aligned sequences
##'
##' 
##' @title printAlignedSeq
##' @param aln alignment
##' @param window window
##' @return NULL
##' @author ygc
##' @export
printAlignedSeq <- function(aln, window=80) {
    n <- aln$length
    idx <- seq(1, n, by=window)
    if (idx[length(idx)] < n) {
        idx <- c(idx, n)
    }
    seq.df <- aln2seqDF(aln)
    nn <- paste(substring(aln$seqs[,1], 1, 20), "...", sep="")
 
    ss2 <- sapply(1:(length(idx)-1), function(i) {
        ss <- seq.df[, idx[i]:idx[i+1]]
        for (k in 1:ncol(ss)) {
            if (length(unique(ss[,k])) == 1) {
                ss[,k] = "."
            }
        }
        apply(ss, 1, paste, collapse="")
    })
    for (i in 1:ncol(ss2)) {
        for (j in 1:nrow(ss2)) {
            cat(nn[j], "  ", ss2[j, i], "\n")
        }
        cat("\n")
    }

    
    p <- which(apply(seq.df, 2, function(i) length(unique(i)) != 1))
    if (length(p) > 0) {
        x <- seq.df[, p]
        colnames(x) <- p
        rownames(x) <- nn
        cat("Ambiguous bases:\n")
        print(x)
    }
}

##' print the consensus sequence
##'
##' 
##' @title printConsensus
##' @param aln alignment
##' @param window window of ambiguous base
##' @return character
##' @author ygc
##' @export
printConsensus <- function(aln, window=10) {
    seq.df <- aln2seqDF(aln)
    cs <- getConsensus(seq.df[-1,])
    ii <- cs %in% c("A", "C", "G", "T")

    print("Consensus sequence:")
    cs2 <- paste(cs, collapse="")

    printFormatSeq(cs2)
    
    if (all(ii) == FALSE) {
        idx <- which(!ii)
        print("Ambiguity bases:")
        for (i in idx) {
            cat("\n\tposition of ", i, " is ", cs[i], "\n")
            pre <- ifelse((i-window) > 0, i-window, 0)
            pro <- ifelse((i+window) < ncol(seq.df), i+window, ncol(seq.df))
            aa <- apply(seq.df[,pre:pro], 1, paste, collapse="")
            adf <- data.frame(seq=aa)
            names(adf) <- paste(paste(rep(" ", i-pre), collapse=""), "|", paste(rep(" ", pro-i), collapse=""), sep="")
            rownames(adf) <- paste(substring(aln$seqs[,1], 1, 20), "...", "")
            cat("\n")
            print(adf)
            cat("\n")
        }
    }
    invisible(cs2)
}

##' print nucleotide code
##'
##' 
##' @title printNTcode
##' @return NULL
##' @author ygc
##' @export
printNTcode <- function() {
    cat("reference:\thttp://www.bios.niu.edu/johns/bioinform/nucleotide_ambiguity_codes.htm\n\n")
    cat("Symbol\tNucleotides\tName\t\t\t\tComplement\n")
    cat("A\tA\t\tadenosine\t\t\tT\n")
    cat("C\tC\t\tcytosine\t\t\tG\n")
    cat("G\tG\t\tguanine\t\t\t\tC\n")
    cat("T\tT\t\tthymidine\t\t\tA\n")
    cat("M\tAC\t\tamino (C & A)\t\t\tK\n")
    cat("R\tAG\t\tpurine (A & G)\t\t\tY\n")
    cat("W\tAT\t\tweak (A & T)\t\t\tW\n")
    cat("S\tCG\t\tstrong (G & C)\t\t\tS\n")
    cat("Y\tCT\t\tpyrimidine (C & T)\t\tR\n")
    cat("K\tGT\t\tketo (T & G)\t\t\tM\n")
    cat("V\tACG\t\tnot T\t\t\t\tB\n")
    cat("H\tACT\t\tnot G\t\t\t\tD\n")
    cat("D\tAGT\t\tnot C\t\t\t\tH\n")
    cat("B\tCGT\t\tnot A\t\t\t\tV\n")
    cat("N\tACGT\t\tany nucleotide\t\t\tX/N\n")
    cat("X\tACGT\t\tany nucleotide or unknown\tX/N\n")
}
