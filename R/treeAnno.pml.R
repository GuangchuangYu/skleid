##' tree annotation of sequence substitution by comparing to parent node
##'
##' 
##' @title treeAnno.pml
##' @param pmlTree tree in pml format, output of optim.pml (phangorn)
##' @param outTree output tree file
##' @param translate logical
##' @param removeGap logical
##' @param plot logical
##' @return tree
##' @importFrom ape read.tree
##' @importFrom ape reorder.phylo
##' @importFrom ape plot.phylo
##' @importFrom ape write.tree
##' @importFrom phangorn ancestral.pml
##' @export
##' @author Yu Guangchuang
treeAnno.pml <- function(pmlTree, outTree="out.nwk", translate=TRUE, removeGap=TRUE, plot=FALSE) {
    ## pmlTree is the output of optim.pml
    ##
    ## infile="sH2N2-98-ref.fas"
    ## tre <- read.tree("H2.nwk")
    ## dat <- read.phyDat(infile,format="fasta")
    ## fit <- pml(tre, dat, k=4)
    ## fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
    ##                  optInv=T, optGamma=T, optEdge=TRUE,
    ##                  optRooted=FALSE, model = "GTR")

    ## anno <- ancestral.pml(pmlTree, "ml")
    ## anno <- matrix2vector.phyDat(anno)
    ## names(anno) already sorted to match 1:(2n-2) that match the phylo class.

    anno <- pmlToSeqString(pmlTree, includeAncestor=TRUE)
    
    
    tr <- pmlTree$tree
    tr <- reorder.phylo(tr, "postorder")
    if (is.null(tr$node.label)) {
        n <- length(tr$tip.label)
        nl <- (n+1):(2*n-2)
        tr$node.label <- as.character(nl)
    }
    
    nodeName <- getNodeName(tr)
    root = tr$edge[nrow(tr$edge),1]
    cn <- getChild(tr, root)
    
    addLabel <- function(tr, node, anno) {
        nodeName <- getNodeName(tr)
        root <- tr$edge[nrow(tr$edge),1]
        cc <- node ## current node
        nodename <- nodeName[cc]
        parent <- getParent(tr, cc)
        parent.nodename <- nodeName[parent]

        labs <- getSubsLabel(anno, parent, cc, translate, removeGap)
        
        if ( ! is.null(labs) ) {
            xx <- paste("[", labs, "]", sep="", collapse="")
            cat(parent.nodename, "->", nodename, ": ", xx, "\n") 
            x <- paste(nodename, xx) 
            if (nodename %in% tr$node.label) {
                tr$node.label[tr$node.label==nodename] <<- x
            } else {
                tr$tip.label[tr$tip.label==nodename] <<-x
            }
        }
    
        children <- getChild(tr, cc)
        for (child in children) {
            addLabel(tr, child, anno)
        }
    }
    
    for (cc in cn) {
        addLabel(tr, cc, anno)
    }
    
    write.tree(tr, file=outTree)
    if (plot) {
        plot.phylo(tr, show.node.label=TRUE)
    }
    
    cat("-> done... \n")
    cat("------------\n")
    mooseSay()

    invisible(tr)
}

##' convert pml object to fasta file or XStringSet object
##'
##' 
##' @title pmlToSeq 
##' @param pml pml object
##' @param includeAncestor logical 
##' @param output out fasta file
##' @return XStringSet
##' @importFrom phangorn ancestral.pml
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings writeXStringSet
##' @export
##' @author ygc
pmlToSeq <- function(pml, includeAncestor=TRUE, output="seq.fa") {
    res <- pmlToSeqString(pml, includeAncestor)
    seqobj <- DNAStringSet(res)
    writeXStringSet(seqobj, output)
    invisible(seqobj)
}

pmlToSeqString <- function(pml, includeAncestor=TRUE) {
    if (includeAncestor == FALSE) {
        phyDat <- pml$data
    } else {
        phyDat <- ancestral.pml(pml, "ml")
    }
    
    phyDat <- matrix2vector.phyDat(phyDat)
    ## defined by phangorn
    labels <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", 
                "y", "k", "v", "h", "d", "b", "n", "?", "-")
    labels <- toupper(labels)

    index <- attr(phyDat, "index")
    
    result <- do.call(rbind, phyDat)
    result <- result[, index, drop=FALSE]

    res <- apply(result, 2, function(i) labels[i])
    res <- apply(res, 1, paste, collapse="")
    names(res) <- rownames(result)
    return(res)
}


seq2codon <- function(x) {
    substring(x, first=seq(1, nchar(x)-2, 3), last=seq(3, nchar(x), 3))
}

##' @importFrom Biostrings GENETIC_CODE
codon2AA <- function(codon) {
    aa <- GENETIC_CODE[codon]
    aa[is.na(aa)] <- "X"
    return(aa)
}

##' @importFrom magrittr %>%
getSubsLabel <- function(seqs, A, B, translate, removeGap) {
    seqA <- seqs[A]
    seqB <- seqs[B]

    if (translate == TRUE) {
        AA <- seqA %>% seq2codon %>% codon2AA
        BB <- seqB %>% seq2codon %>% codon2AA
        ## codonA <- seq2codon(seqA)
        ## codonB <- seq2codon(seqB)
        
        ## AA <- codon2AA(codonA)
        ## BB <- codon2AA(codonB)
    } else {
        n <- nchar(seqA) ## should equals to nchar(seqB)
        AA <- substring(seqA, 1:n, 1:n)
        BB <- substring(seqB, 1:n, 1:n)
    }
    
    ii <- which(AA != BB)

    if (removeGap == TRUE) {
        if (length(ii) > 0 && translate == TRUE) {
            ii <- ii[AA[ii] != "X" & BB[ii] != "X"]
        }

        if (length(ii) > 0 && translate == FALSE) {
            ii <- ii[AA[ii] != "-" & BB[ii] != "-"]
        }
    }
    
    if (length(ii) == 0) {
        return(NULL)
    }
    
    res <- paste(AA[ii], ii, BB[ii], sep="", collapse="/")
    return(res)
}

matrix2vector.phyDat <- function(x) {
    index <- attr(x, "index")
    res <- lapply(x, matrix2vector.phyDat.item)
    names(res) <- names(x)
    attr(res, "index") <- index
    class(res) <- "phyDat"
    return(res)
}

matrix2vector.phyDat.item <- function(y) {
    ii <- apply(y, 1, function(xx) {
        ## return index of a c g and t, if it has highest probability
        ## otherwise return index of -
        jj <- which(xx == max(xx))
        if ( length(jj) > 1) {
            if (length(jj) < 4) {
                warning("ambiguous found...\n")
            } else {
                ## cat("insertion found...\n")
            }
            ## 18 is the gap(-) index of base character defined in phangorn
            ## c("a", "c", "g", "t", "u", "m", "r", "w", "s", 
	    ##   "y", "k", "v", "h", "d", "b", "n", "?", "-")
            18
        } else {
            jj
        }
    })
    unlist(ii)
}

