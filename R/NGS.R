##' Class "NGS"
##' This class stores information for N-glycosylation sites
##'
##'
##' @name NGS-class
##' @aliases NGS-class
##'       show,NGS-method toCodon,NGS-method
##'       toAA,NGS-method toPosition,NGS-method
##'       getPosition,NGS-method getSite,NGS-method
##'       toString,NGS-method
##'
##' @docType class
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

##' read fasta sequences and construct NGS object
##'
##' 
##' @title toNGS 
##' @param fasta fasta file 
##' @param start start position, for translation
##' @return NGS object
##' @export
##' @author Yu Guangchuang
toNGS <- function(fasta, start=1) {
    fa <- readBStringSet(fasta)
    res <- new("NGS", stringSet = fa,
               start = start)
    res@characterSet <- toString(res)
    res@codonSet <- toCodon(res)
    res@AASet <- toAA(res)
    res@position <- toPosition(res)
    return(res)
}

##' show method for \code{NGS} instance
##'
##'
##' @name show
##' @docType methods
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
              cat("\t use _getPos_ to access the NGS positions\n")
              cat("\t use _getSite_ to access the NGS sites\n")
          })

          
##' toCodon method
##'
##'
##' @docType methods
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
##' @docType methods
##' @name toAA
##' @rdname toAA-methods
##' @aliases toAA,NGS,ANY-method
##' @title toAA method
##' @param x NGS object
##' @param ... additional parameters
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
##' @docType methods
##' @name toString
##' @rdname toString-methods
##' @aliases toString,NGS,ANY-method
##' @title toString method
##' @param x NGS object
##' @param ... additional parameters
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
##' @docType methods
##' @name toPosition
##' @rdname toPosition-methods
##' @aliases toPosition,NGS,ANY-method
##' @title toPosition method
##' @param x NGS object
##' @param ... additional parameters
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
##' @docType methods
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
##' @docType methods
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









