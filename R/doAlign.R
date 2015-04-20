##' align 454 contigs and mira sequence to reference sequence.
##'
##' 
##' @title doAlign
##' @param x file names of 454 contig, reference and mira sequences
##' @return fasta S3 object
##' @author ygc
##' @importFrom muscle muscle
##' @importFrom Biostrings readDNAStringSet
##' @importFrom Biostrings reverseComplement
##' @importFrom magrittr %<>%
##' @export
doAlign <- function(x) {
    ## suppose x[1] is 454 contigs
    ## x[2] is reference sequence
    ## x[3] is mira sequence, that is read assembly that mapping to reference.

    ## f454 <- read.fasta(x[1])
    f454 <- readDNAStringSet(x[1])
    
    fref <- readDNAStringSet(x[2])
    fref <- fref[1]
    ## fref$seqs <- fref$seqs[1,]
    ## fref$num <- 1
    
    ## fmira <- read.fasta(x[3])
    fmira <- readDNAStringSet(x[3])
    
    ## for (i in 1:nrow(f454$seqs)) {
    for (i in 1:length(f454)) {
        ## fa <- list(seqs=rbind(f454$seqs[i,], fref$seqs), num=fref$num+1)
        ## class(fa) <- "fasta"
        ## aln <- muscle(fa, quiet = TRUE)
        fa <- c(f454[i], fref)
        aln <- muscle(fa, quiet=TRUE)
        
        fa2 <- fa
        fa2[1] <- reverseComplement(fa2[1])

        ## fa2$seqs[1,2] <- revcomp(fa2$seqs[1,2])
        aln2 <- muscle(fa2, quiet = TRUE)
        if (identityRatio(aln) < identityRatio(aln2)) {
            ## f454$seqs[i,1] <- paste("RC", f454$seqs[i,1], sep="_")
            ## f454$seqs[i,2] <- revcomp(f454$seqs[i,2])
            f454[i] <- reverseComplement(f454[i])
            names(f454)[i] <- paste("RC", names(f454)[i], sep="_")
        }
    }

    ## fasta <- list(seqs=rbind(fref$seqs, fmira$seqs, f454$seqs),
    ##               num=fref$num+fmira$num+f454$num)
    ## class(fasta) <- "fasta"

    fasta <- c(fref, fmira, f454)
    
    ## if (f454$num == 1) {
    if (length(f454) == 1) {
        res2 <- muscle(fasta, quiet = TRUE)
        ## ii <- getIdx(fasta$seqs[,1], res2$seqs[,1])
        ## res2$seqs <- res2$seqs[ii,]
        result <- DNAStringSet(res2)
        ii <- getIdx(names(fasta), names(result))
        result <- result[ii]
        return(result)
    }

    ## res <- lapply(1:f454$num, function(i) {
    res <- lapply(1:length(f454), function(i) {
        ## fa <- list(seqs=rbind(fref$seqs, fmira$seqs, f454$seqs[i,]),
        ##            num=3)
        ## class(fa) <- "fasta"
        fa <- c(fref, fmira, f454[i])
        xx <- muscle(fa, quiet = TRUE)
        yy <- DNAStringSet(xx)
        ## jj <- getIdx(fasta$seqs[1,1], xx$seqs[,1])
        jj <- getIdx(names(fasta)[1], names(yy))
        ## ii <- 1:nrow(xx$seqs)
        ii <- 1:length(yy)
        yy[c(jj, ii[-jj])]
        ## xx$seqs <- xx$seqs[c(jj, ii[-jj]),]
        ## return(xx)
        return(yy)
    })

    ## if (length(unique(sapply(res, function(x) x$length))) != 1) {
    if (length(unique(sapply(res, width)[1,])) != 1) {
        ## yy <- list(seqs=do.call("rbind", lapply(res, function(x) x$seqs[1,])),
        ##            num=length(x))
        ## yy$seqs[,2] %<>% gsub("-", "X", .)
        ## class(yy) <- "fasta"
        yy <- toString(res[[1]][1])
        for (i in 2:length(res)) {
            yy <- c(yy, toString(res[[i]][1]))
        }
        yy %<>% gsub("-", "N", .)
        yy <- DNAStringSet(yy)
        yres <- muscle(yy, quiet = TRUE)
        yres <- DNAStringSet(yres)
        for(j in seq_along(res)) {
            res_item <- res[[j]]
            ## res_seq <- res_item$seqs[,2]
            res_seq <- sapply(res_item, toString)
            
            ## seq <- yres$seqs[j, 2]
            seq <- toString(yres[j])
            
            idx <- gregexpr("-+", seq)
            idx <- idx[[1]]
            
            if (length(idx) == 1 && idx < 1) {
                next
            }
            
            for (i in seq_along(idx)) {
                for (k in seq_along(res_seq)) {
                    if (idx[i] == 1) {
                        res_seq[k] <- paste0(
                            paste0(rep("-", attr(idx, "match.length")[i]), collapse = ""),
                            res_seq[k], collapse = "")
                    } else if (idx[i] > nchar(res_seq[k])) {
                        res_seq[k] <- paste0(res_seq[k],
                                               paste0(rep("-", attr(idx, "match.length")[i]), collapse = ""),
                                               collapse = "")
                    } else {
                        res_seq[k] <- paste0(substring(res_seq[k], 1, idx[i]-1),
                                             paste0(rep("-", attr(idx, "match.length")[i]), collapse = ""),
                                             substring(res_seq[k], idx[i]), collapse = "")
                    }
                }
            }

            ## res_item$seqs[,2] <- res_seq
            ## res[[j]] <- res_item
            res[[j]] <- DNAStringSet(res_seq)
            ## res[[j]]$length <- nchar(res_seq[1])
        }
    }
    
    ## if (length(unique(sapply(res, function(x) x$length))) == 1) {
    if (length(unique(sapply(res, width)[1,])) == 1) {
        ## jj <- sapply(res, function(x) getIdx(fasta$seqs[1,1], x$seqs[,1]))
        jj <- sapply(res, function(x) getIdx(names(fasta)[1], names(x)))

        ## r1 <- res[[1]]$seqs[jj[1],2]
        r1 <- res[[1]][jj[1]]
        flag <- FALSE
        for (i in 2:length(jj)) { ## jj should >= 2
            ## if (r1 != res[[i]]$seqs[jj[i],2]) {
            if (r1 != res[[i]][jj[i]]) {
                flag <- TRUE
                break
            }
        }
        if (flag == TRUE) {
            res2 <- muscle(fasta, quiet = TRUE)
            result <- DNAStringSet(res2)
        } else {
            ## seqs <- lapply(res, function(x) x$seqs)
            ## seqs <- do.call("rbind", seqs)
            ## seqs <- unique(seqs)
            seqs <- res[[1]]
            for (i in 2:length(res)) {
                seqs <- c(seqs, res[[i]])
            }
            seqs <- unique(seqs)
            
            ## sn <- seqs[,1]
            sn <- names(seqs)
         
            if (length(sn) != length(unique(sn))) {
                k <- sapply(unique(sn), function(i) which(i == sn)[1])
                ## seqs <- seqs[k,]
                seqs <- seqs[k]
            }
            ## res2 <- list(seqs=seqs, num=nrow(seqs), length=res[[1]]$length)
            result <- seqs
        }
    } else {
        res2 <- muscle(fasta, quiet = TRUE)
        result <- DNAStringSet(res2)
    }
    
    ## ii <- getIdx(fasta$seqs[,1], res2$seqs[,1])
    ii <- getIdx(names(fasta), names(result))
    ## res2$seqs <- res2$seqs[ii,]
    result <- result[ii]
    return(result)
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
    fa <- seqs[[1]]
    if (length(seqs) > 1) {
        for (i in 2:length(seqs)) {
            fa <- c(fa, seqs[[i]])
        }
    }
    
    ## str <- unlist(lapply(seqs, toString))
    ## nn <- unlist(lapply(seqs, names))
    ## fa <- list(seqs=data.frame(V1=nn, V2=str), num=length(nn))
    ## class(fa) <- "fasta"
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
    ## fa <- DNAStringSet(aln$seqs[,2])
    ## names(fa) <- aln$seqs[,1]
    ## writeXStringSet(fa, output)
    writeXStringSet(DNAStringSet(aln), output)
}
