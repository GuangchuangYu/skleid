##' summarize the result
##'
##' 
##' @title summarize
##' @param out.folder out folder
##' @param name.file name file 
##' @return NULL
##' @export
##' @author ygc
summarize <- function(out.folder, name.file) {
    nameMap <- read.delim(name.file, header=F, stringsAsFactor=FALSE)
    colnames(nameMap) <- c("ID", "name")
    cs <- list.files(path=out.folder, pattern=".fas$")
    cs2 <- paste(out.folder, cs, sep="/")

    subtype <- getSubtype(nameMap, cs)
    nameMap$subtype <- subtype

    kw <- c("HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2")
    asite <- sapply(kw, getAmbiguousSite, nameMap=nameMap, cs2=cs2)

    res <- cbind(nameMap, asite)
    write.csv(res, file="summary.csv", row.names=F)
}

getAmbiguousSite <- function(nameMap, cs2, keyword) {
    if (keyword == "HA") {
        keyword = "H\\d"
    } else if (keyword == "NA") {
        keyword = "N\\d"
    }

    result <- sapply(1:nrow(nameMap), function(i) {
        strain <- cs2[grep(nameMap[i,2], cs2)]
        if (length(strain) == 0) {
            return("")
        }
        gg <- strain[grep(paste(keyword, "[_mixed]*\\.fas$", sep=""), strain)]
        if (length(gg) == 0) {
            return("")
        }
        x <- toupper(toString(readBStringSet(gg)))
        x <- substring(x, 1:nchar(x), 1:nchar(x))
        ii <- which(! x %in% c("A", "C", "G", "T"))
        if (length(ii) > 0) {
            ii <- ii[ii > 10 && ii < (length(x)-10)]
        }
        if (length(ii) == 0) {
            return("full")
        }
        res <- paste(ii, x[ii], sep="")
        if (length(ii) > 1) {
            res <- paste(res, sep=",", collapse = ",")
        }
        return(res)
    })
    return(unlist(result))
}

getSubtype <- function(nameMap, cs) {
    res <- sapply(1:nrow(nameMap), function(i) {
        strain <- cs[grep(nameMap[i,2], cs)]
        if (length(strain) == 0) {
            return("")
        }
        hh <- strain[grep("H\\d[_mixed]*\\.fas$", strain)]
        nn <- strain[grep("N\\d[_mixed]*\\.fas$", strain)]
        ht <- gsub(".*(H\\d)[_mixed]*\\.fas", "\\1", hh)
        nt <- gsub(".*(N\\d)[_mixed]*\\.fas", "\\1", nn)
        subtype <- paste(ht, nt, sep="")
        if (length(subtype) == 0)
            return("")
        return(subtype)
    })
    return(unlist(res))
}
