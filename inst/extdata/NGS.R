##' Class "NGS"
##' This class stores information for N-glycosylation sites
##'
##'
##' @name NGS-class
##' @aliases NGS-class
##'       show,NGS-method toCodon,NGS-method
##'       toAA,NGS-method toPosition,NGS-method
##'       getPosition,NGS-method getSite,NGS-method
##'       toString,NGS-method setNGS_,NGS-method
##'
##' @slot stringSet BStringSet
##' @slot start start position
##' @slot characterSet sequence character vector
##' @slot codonSet codon of sequences
##' @slot AASet translate codonSet to AASet
##' @slot position NGS position
##' @exportClass NGS
##' @author GuangchuangYu \url{http://ygc.name}
##' @keywords classes
setClass("NGS",
         representation = representation(
             stringSet = "BStringSet",
             start = "numeric",
             characterSet = "character",
             codonSet = "matrix",
             AASet = "matrix",
             position = "list"
             ),
         prototype = prototype(start = 1)
         )


seq2codon <- ggtree:::seq2codon
codon2AA <- ggtree:::codon2AA

##' read fasta sequences and construct NGS object
##'
##' 
##' @title toNGS 
##' @param fasta fasta file 
##' @param start start position, for translation
##' @return NGS object
##' @importFrom magrittr %>%
##' @export
##' @author Yu Guangchuang
toNGS <- function(fasta, start=1) {
    new("NGS",
        stringSet = readBStringSet(fasta),
        start = start) %>%
            setNGS_
}

##' setNGS_ method
##'
##'
##' @name setNGS_
##' @rdname setNGS_-methods
##' @aliases setNGS_,NGS,ANY-method
##' @title setNGS_ method
##' @param x NGS object
##' @param ... additional parameter
##' @return NGS object
##' @exportMethod toCodon
##' @author Guangchuang Yu
##' @usage setNGS_(x, ...)
setMethod("setNGS_", signature(x="NGS"),
          function(x, ...) {
              x@characterSet <- toString(x)
              x@codonSet <- toCodon(x)
              x@AASet <- toAA(x)
              x@position <- toPosition(x)
              return(x)
          })


##' show method for \code{NGS} instance
##'
##'
##' @name show
##' @rdname show-methods
##'
##' @title show method
##' @param object A \code{NGS} instance
##' @return print info
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuanbg Yu
setMethod("show", signature(object = "NGS"),
          function(object) {
              cat("NGS class:\n")
              cat("\t", length(object@characterSet), "sequences",
                  "of length", nchar(object@characterSet[1]), "\n")
              cat("\t use _getPosition_ to access the NGS positions\n")
              cat("\t use _getSite_ to access the NGS sites\n")
          })

          
##' toCodon method
##'
##'
##' @name toCodon
##' @rdname toCodon-methods
##' @aliases toCodon,NGS,ANY-method
##' @title toCodon method
##' @param x NGS object
##' @param ... additional parameter
##' @return matrix
##' @exportMethod toCodon
##' @author Guangchuang Yu
##' @usage toCodon(x, ...)
setMethod("toCodon", signature(x="NGS"),
          function(x, ...) {
              seqs <- x@characterSet
              codon <- t(sapply(seqs, seq2codon))
              rownames(codon) <- NULL
              return(codon)
          })


##' toAA method
##'
##'
##' @name toAA
##' @rdname toAA-methods
##' @aliases toAA,NGS,ANY-method
##' @title toAA method
##' @inheritParams toCodon
##' @return matrix
##' @exportMethod toAA
##' @author Guangchuang Yu
##' @usage toAA(x, ...)
setMethod("toAA", signature(x="NGS"),
          function(x, ...) {
              codon <- x@codonSet
              matrix(codon2AA(codon), ncol=ncol(codon))  
          })


##' toString method
##'
##'
##' @name toString
##' @rdname toString-methods
##' @aliases toString,NGS,ANY-method
##' @title toString method
##' @inheritParams toCodon
##' @return matrix
##' @importFrom Biostrings toString
##' @exportMethod toString
##' @author Guangchuang Yu
##' @usage toString(x, ...)
setMethod("toString", signature(x = "NGS"),
          function(x, ...) {
              toString.NGS(x, ...)
          })


##' toPosition method
##'
##'
##' @name toPosition
##' @rdname toPosition-methods
##' @aliases toPosition,NGS,ANY-method
##' @title toPosition method
##' @inheritParams toCodon
##' @return list
##' @exportMethod toPosition
##' @author Guangchuang Yu
##' @usage toPosition(x, ...)
setMethod("toPosition", signature(x="NGS"),
          function(x, ...) {
              AASet <- x@AASet
              ngs <- apply(AASet, 1, position.NGS)
              names(ngs) <- names(x@stringSet)
              return(ngs)
          })

##' getPosition method
##'
##'
##' @name getPosition
##' @rdname getPosition-methods
##' @aliases getPosition,NGS,ANY-method
##' @title getPosition method
##' @param x NGS object
## @param to output format
##' @param ... additional parameters
##' @return NGS positions
##' @exportMethod getPosition
##' @author Guangchuang Yu
##' @usage getPosition(x, ...)
setMethod("getPosition", signature(x="NGS"),
          function(x, to="data.frame", ...) {
              pos <- x@position
              if (to == "data.frame") {
                  pos <- to.data.frame.NGS(pos)
                  colnames(pos) <- c("seq", "pos")
              }
              return(pos)
          })


##' getSite method
##'
##'
##' @name getSite
##' @rdname getSite-methods
##' @aliases getSite,NGS,ANY-method
##' @title getSite method
##' @param x NGS object
## @param to output format
##' @param ... additional parameters
##' @return NGS site
##' @exportMethod getSite
##' @author Guangchuang Yu
##' @usage getSite(x, ...)
setMethod("getSite", signature(x="NGS"),
          function(x, to="data.frame", ...) {
              pos <- getPosition(x, to="list")
              AASet <- x@AASet
              sites <- sapply(1:length(pos), function(i) {
                  aa <- AASet[i,]
                  idx <- pos[[i]]
                  site <- matrix(aa[sapply(idx, function(i) i:(i+2))],
                                 byrow=T, ncol=3)
                  site <- apply(site, 1, paste, collapse="")
                  return(site)
              })
              names(sites) <- names(pos)
              if (to == "data.frame") {
                  sites <- to.data.frame.NGS(sites)
                  colnames(sites) <- c("seq", "site")
              }
              return(sites)
          })



position.NGS <- function(aa) {
    idx <- which(aa == "N")
    idx <- idx[which(aa[idx+1] != "P")]
    idx <- idx[which(aa[idx+2] == "S" | aa[idx+2] == "T")]
    return(idx)
}

toString.NGS <- function(x, ...) {
    start <- x@start
    ss <- x@stringSet
    n <- length(ss)
    seqs <- sapply(1:n, function(i) {
        s <- toString(ss[i])
        substring(s, start, nchar(s))
    })
    return(seqs)
}

to.data.frame.NGS <- function(LIST) {
    nn <- rep(names(LIST), times=sapply(LIST, length))
    df <- data.frame(V1=nn, V2=unlist(LIST))
    row.names(df) <- NULL
    return(df)
}









