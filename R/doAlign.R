##' align 454 contigs and mira sequence to reference sequence.
##'
##' 
##' @title doAlign
##' @param x file names of 454 contig, reference and mira sequences
##' @return fasta S3 object
##' @author ygc
##' @importFrom muscle read.fasta
##' @importFrom muscle muscle
##' @export
doAlign <- function(x) {
    ## suppose x[1] is 454 contigs
    ## x[2] is reference sequence
    ## x[3] is mira sequence, that is read assembly that mapping to reference.
    f454 <- read.fasta(x[1])

    fref <- read.fasta(x[2])
    fref$seqs <- fref$seqs[1,]
    fref$num <- 1
    
    fmira <- read.fasta(x[3])
    
    for (i in 1:nrow(f454$seqs)) {
        fa <- list(seqs=rbind(f454$seqs[i,], fref$seqs), num=fref$num+1)
        class(fa) <- "fasta"
        aln <- muscle(fa, quiet = TRUE)

        fa2 <- fa
        fa2$seqs[1,2] <- revcomp(fa2$seqs[1,2])
        aln2 <- muscle(fa2, quiet = TRUE)
        if (identityRatio(aln) < identityRatio(aln2)) {
            f454$seqs[i,1] <- paste("RC", f454$seqs[i,1], sep="_")
            f454$seqs[i,2] <- revcomp(f454$seqs[i,2])
        }
    }

    fasta <- list(seqs=rbind(fref$seqs, fmira$seqs, f454$seqs),
                  num=fref$num+fmira$num+f454$num)
    class(fasta) <- "fasta"
    
    if (f454$num == 1) {
        res2 <- muscle(fasta, quiet = TRUE)
        ii <- getIdx(fasta$seqs[,1], res2$seqs[,1])
        res2$seqs <- res2$seqs[ii,]
        return(res2)
    }
    
    res <- lapply(1:f454$num, function(i) {
        fa <- list(seqs=rbind(fref$seqs, fmira$seqs, f454$seqs[i,]),
                   num=3)
        class(fa) <- "fasta"
        muscle(fa, quiet = TRUE)
    })
    
    if (length(unique(sapply(res, function(x) x$length))) == 1) {
        jj <- sapply(res, function(x) getIdx(fasta$seqs[1,1], x$seqs[,1]))

        r1 <- res[[1]]$seqs[jj[1],2]
        flag <- FALSE
        for (i in 2:length(jj)) { ## jj should >= 2
            if (r1 != res[[i]]$seqs[jj[i],2]) {
                flag <- TRUE
                break
            }
        }
        if (flag == TRUE) {
            res2 <- muscle(fasta, quiet = TRUE)
        } else {
            seqs <- lapply(res, function(x) x$seqs)
            seqs <- do.call("rbind", seqs)
            seqs <- unique(seqs)
            sn <- seqs[,1]
            if (length(sn) != length(unique(sn))) {
                k <- sapply(unique(sn), function(i) which(i == sn)[1])
                seqs <- seqs[k,]
            }
            res2 <- list(seqs=seqs, num=nrow(seqs), length=res[[1]]$length)
        }
    } else {
        res2 <- muscle(fasta, quiet = TRUE)
    }
    
    ii <- getIdx(fasta$seqs[,1], res2$seqs[,1])
    res2$seqs <- res2$seqs[ii,]
    return(res2)
}


##' align multiple sequence stored in several fasta files.
##'
##' 
##' @title doAlign2
##' @param files fasta files
##' @return alignment 
##' @author ygc
##' @importFrom Biostrings readDNAStringSet
##' @importFrom muscle muscle
##' @export
doAlign2 <- function(files) {
    seqs <- lapply(files, readDNAStringSet)
    str <- unlist(lapply(seqs, toString))
    nn <- unlist(lapply(seqs, names))
    fa <- list(seqs=data.frame(V1=nn, V2=str), num=length(nn))
    class(fa) <- "fasta"
    aln <- muscle(fa)
    return(aln)
}

##' write aligned sequence to fasta file
##'
##' 
##' @title writeAlignedSeq
##' @param aln alignment
##' @param output out file
##' @return NULL
##' @author ygc
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings writeXStringSet
##' @export
writeAlignedSeq <- function(aln, output) {
    fa <- DNAStringSet(aln$seqs[,2])
    names(fa) <- aln$seqs[,1]
    writeXStringSet(fa, output)
}
