##' @title hamming
##' @param seq1 sequence 1
##' @param seq2 sequence 2
##' @param indel include indel or not
##' @param ... additional parameter
##' @rdname hamming
##' @export
hamming <- function(seq1, seq2, indel=FALSE, ...) {
    UseMethod("hamming")
}

##' @docType methods
##' @name toCodon
##' @rdname toCodon-methods
##' @title toCodon method
##' @export
setGeneric("toCodon", function(x, ...) standardGeneric("toCodon"))

##' @docType methods
##' @name toAA
##' @rdname toAA-methods
##' @title toAA method
##' @export
setGeneric("toAA", function(x, ...) standardGeneric("toAA"))

##' @docType methods
##' @name toPosition
##' @rdname toPosition-methods
##' @title toPosition method
##' @export
setGeneric("toPosition", function(x, ...) standardGeneric("toPosition"))

##' @docType methods
##' @name getPosition
##' @rdname getPosition-methods
##' @title getPosition method
##' @export
setGeneric("getPosition", function(x, ...) standardGeneric("getPosition"))

##' @docType methods
##' @name getSite
##' @rdname getSite-methods
##' @title getSite method
##' @export
setGeneric("getSite", function(x, ...) standardGeneric("getSite"))

##' @docType methods
##' @name setNGS_
##' @rdname setNGS_-methods
##' @title setNGS_ method
##' @export
setGeneric("setNGS_", function(x, ...) standardGeneric("setNGS_"))

