##' annotating phylogenetic tree with functional site mutation
##'
##' 
##' @title treeAnno
##' @param inTree input tree in nwk format
##' @param outTree output tree in nwk format
##' @param anno site annotation, e.g. N-glyco sites 
##' @param plot boolean variable, plot or not.
##' @importFrom ape read.tree
##' @importFrom ape reorder.phylo
##' @importFrom ape plot.phylo
##' @importFrom ape write.tree
##' @export
##' @return tree
##' @author Guangchuang Yu
treeAnno <- function(inTree, outTree="out.nwk", anno, plot=FALSE) {
    tr <- read.tree(inTree)
    tr <- reorder.phylo(tr, "postorder")
    
    nodeName <- c(tr$tip.label, tr$node.label)
    root = tr$edge[nrow(tr$edge),1]
    cn <- getChild(tr, root)

    addLabel <- function(tr, node, nn) {
        nodeName <- c(tr$tip.label, tr$node.label)
        root <- tr$edge[nrow(tr$edge),1]
        cc <- node ## current node
        nodename <- nodeName[cc]
        pos <- getPos(nn, nodename)
        x <- nodename
        if (is.null(pos)) {
            x <- paste(nodename, "[*]", sep="",collpase="")
        } else if (!is.null(pos) && getParent(tr, cc) == root) {
            x <- paste(c(nodename, " [", paste(pos, collapse="/"), "]"),
                       sep="", collapse="")
        } else {
            parent <- getParent(tr, cc)
            parent.nodename <- nodeName[parent]
            parent.pos <- getPos(nn, parent.nodename)
            if (!is.null(parent.pos)) {
                plus <- pos[!pos %in% parent.pos]
                if (length(plus) != 0) {
                    gain <- paste(paste("+", plus, sep=""), collapse="/")
                } else {
                    gain <- ""
                }
                neg <- parent.pos[! parent.pos %in% pos]
                if (length(neg)!=0) {
                    loss <- paste(paste("-", neg, sep=""), collapse="/")
                } else {
                    loss <- ""
                }
                if (length(plus) == 0 && length(neg) == 0) {
                    x <- paste(nodename, "[-]", collapse="")
                } else if ( length(plus) == 0) {
                    x <- paste(nodename, " [", loss, "]", sep="", collaspe="")
                } else if (length(neg) == 0 ) {
                    x <- paste(nodename, " [",  gain, "]", sep="", collaspe="")
                } else {
                    x <- paste(nodename, " [",  gain, "|",loss, "]", sep="", collaspe="")
                }    
            } else {
                x <- paste(c(nodename, " [", paste(pos, collapse="/"), "]"),
                           sep="", collapse="")
            }
        }
        
        if (nodename %in% tr$node.label) {
            tr$node.label[tr$node.label==nodename] <<- x
        } else {
            tr$tip.label[tr$tip.label==nodename] <<-x
        }
        
        children <- getChild(tr, cc)
        for (child in children) {
            addLabel(tr, child, nn)
        }
    }
    
    for (cc in cn) {
        addLabel(tr, cc, anno)
    }
    
    write.tree(tr, file=outTree)
    if (plot)
        plot.phylo(tr, show.node.label=TRUE)
    invisible(tr)
}


getChild <- function(tr, node) {
    ## tree should be in postorder
    tr$edge[tr$edge[,1] == node,2]
}

getParent <- function(tr, node) {
    ## tree should be in postorder
    tr$edge[tr$edge[,2] == node,1]
}

getPos <- function(anno, nodename) {
    pos <- anno[tolower(anno$name)== tolower(nodename), "position"]
    return(sort(pos))
}
