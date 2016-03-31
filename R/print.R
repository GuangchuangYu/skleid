shorten_name <- function(nn) {
    ## in old version of mira sequence
    ##
    ## header is
    ## >PB2_EPI529527_A_BRISBANE_165_2013(H3N2)_BB_1 PB2_EPI529527_A_BRISBANE_165_2013(H3N2)_BB
    ##
    ## which contains _BB
    ##
    ##
    
    ## idx_mira <- grep("_BB", nn)


    ## new version is
    ##
    ## >PA_JX485431_A_duck_Shanghai_28-1_2009(H4N2),_1..1034__length=1036___numreads=46
    ##
    ## which is similar to ref with numreads recorded
    ## numreads also recorded in 454 fas file, which name by contig00X
    ##
    ## 
    ##
    ##
    ## there may exists similar name:
    ## >PA_JX485431_A_duck_Shanghai_28-1_2009(H4N2),_1144..2233__length=1102___numreads=168
    ##
    ## so that previous solution can't solve this issue.
    ##
    ## The new solution is to keep the start..end position and put it at the beginning of the name
    ##

    ii <- grep("\\d+\\.\\.\\d+", nn)
    if (length(ii)) {
        position <- gsub(".*_(\\d+\\.\\.\\d+)_.*", "\\1", nn[ii])
        nn[ii] <- paste(position, nn[ii], sep="_")
    }

    ## idx_mira <- which(!grepl("contig", nn) & grepl("numreads", nn))
    ## if (length(idx_mira) > 0) {
    ##     nn[idx_mira] <- paste0("BB_", nn[idx_mira])
    ## }
    
    nn <- paste(substring(nn, 1, 20), "...", sep="")
    return(nn)
}

##' print aligned sequences
##'
##' 
##' @title printAlignedSeq
##' @param aln alignment
##' @param window window
##' @param printAmbiguous logical
##' @importMethodsFrom Biostrings width
##' @return NULL
##' @author ygc
##' @export
printAlignedSeq <- function(aln, window=80, printAmbiguous=FALSE) {
    ## n <- aln$length
    n <- unique(width(aln))
    idx <- seq(1, n, by=window)
    if (idx[length(idx)] < n) {
        idx <- c(idx, n)
    }
    seq.df <- aln2seqDF(aln)
    ## nn <- paste(substring(aln$seqs[,1], 1, 20), "...", sep="")
    nn <- shorten_name(names(aln))
     
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

    if (printAmbiguous) {
        p <- which(apply(seq.df, 2, function(i) {
            x <- unique(i)
            res <- length(x) != 1
            if (res) {
                if (length(x) == 2 && "-" %in% x) {
                    res <- FALSE
                }
            }
            return(res)
        }))
        
        if (length(p) > 0) {
            if (length(p) == 1) {
                x <- data.frame(a=seq.df[,p])
                colnames(x) <- p
                rownames(x) <- nn
                cat("Ambiguous bases:\n")
                print(x)
            } else {
                x <- seq.df[, p]
                colnames(x) <- p
                rownames(x) <- nn
                cat("Ambiguous bases:\n")
                options(width=70)
                print(x)
            }
        }
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
    cs <- getConsensus(seq.df[-1,, drop=FALSE])
   
    idx <- attr(cs, "index")
    if (length(idx) > 0) {
        seq.df <- seq.df[, idx]
    }
    
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
            ## rownames(adf) <- paste(substring(aln$seqs[,1], 1, 20), "...", "")
            ## rownames(adf) <- paste(substring(names(aln), 1, 20), "...", "")
            rownames(adf) <- shorten_name(names(aln))
            
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
