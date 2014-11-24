##' tree annotation of sequence substitution by comparing to parent node
##'
##' 
##' @title treeAnno.pml
##' @param pmlTree tree in pml format, output of optim.pml (phangorn)
##' @param outTree output tree file
##' @param plot logical
##' @return tree
##' @importFrom ape read.tree
##' @importFrom ape reorder.phylo
##' @importFrom ape plot.phylo
##' @importFrom ape write.tree
##' @importFrom phangorn ancestral.pml
##' @export
##' @author Yu Guangchuang
treeAnno.pml <- function(pmlTree, outTree="out.nwk", plot=FALSE) {
    ## pmlTree is the output of optim.pml
    ##
    ## infile="sH2N2-98-ref.fas"
    ## tre <- read.tree("H2.nwk")
    ## dat <- read.phyDat(infile,format="fasta")
    ## fit <- pml(tre, dat, k=4)
    ## fit <- optim.pml(fit, optNni=FALSE, optBf=T, optQ=T,
    ##                  optInv=T, optGamma=T, optEdge=TRUE,
    ##                  optRooted=FALSE, model = "GTR")

    anno <- ancestral.pml(pmlTree, "ml")
    anno <- matrix2vector.phyDat(anno)
    ## names(anno) already sorted to match 1:(2n-2) that match the phylo class.
    
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

        labs <- getSubsLabel(anno, parent, cc)
        
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

    seq <- lapply(phyDat, function(i) labels[i])
    seq <- lapply(seq, paste, collapse="")
    seqobj <- DNAStringSet(unlist(seq))
    writeXStringSet(seqobj, output)
    invisible(seqobj)
}


getSubsLabel <- function(phyDat, A, B) {
    seqA <- phyDat[[A]]
    seqB <- phyDat[[B]]
    ii <- which(seqA != seqB)
    if (length(ii) == 0) {
        return(NULL)
    }
        
    ## defined by phangorn
    labels <- c("a", "c", "g", "t", "u", "m", "r", "w", "s", 
                "y", "k", "v", "h", "d", "b", "n", "?", "-")
    labels <- toupper(labels)
    res <- paste(labels[seqA[ii]], ii, labels[seqB[ii]], sep="", collapse="/")
    return(res)
}

matrix2vector.phyDat <- function(x) {
    res <- lapply(x, matrix2vector.phyDat.item)
    names(res) <- names(x)
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

