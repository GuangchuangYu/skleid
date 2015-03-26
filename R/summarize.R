##' summarize the result
##'
##' 
##' @title summarize
##' @param out.folder out folder
##' @param name.file name file 
##' @return NULL
##' @importFrom magrittr %>%
##' @export
##' @author ygc
summarize <- function(out.folder, name.file) {
    nameMap <- read.delim(name.file, header=F, stringsAsFactor=FALSE)
    colnames(nameMap) <- c("ID", "name")
    cs <- list.files(path=out.folder, pattern=".fas$")
    cs2 <- paste(out.folder, cs, sep="/")
    cat("-> determing subtypes\t\t\t", format(Sys.time(), "%Y-%m-%d %X"), "\n")
    subtype <- getSubtype(nameMap, cs)

    kw <- c("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")
    asite <- sapply(kw, getAmbiguousSite, nameMap=nameMap, cs2=cs2)

    ## isMixed <- apply(asite, 1, function(i) any(i == ""))

    Mixed <- rep("unmixed", nrow(asite))

    ## mix <- getFiles("mixed")
    ## mix <- gsub("mixed/S\\d+_", "", mix)
    ## mix <- gsub("_S\\d+.*", "", mix) %>% unique
    mixID <- gsub("^S", "", getMixedStrain() %>% unique) %>% as.numeric    
    ## Mixed[isMixed] <- "mixed"
    Mixed[match(mixID, nameMap[,1])] <- "mixed"
    
    nameMap$Mixed <- Mixed
    nameMap$subtype <- subtype
    
    res <- cbind(nameMap, asite)
    write.csv(res, file="summary.csv", row.names=F)
    cat("-> output info to summary.csv\t\t", format(Sys.time(), "%Y-%m-%d %X"), "\n")
    cat(">> done... \n")
    cat("------------\n")
    says()
}

getAmbiguousSite <- function(nameMap, cs2, keyword) {
    cat("-> get ambiguous site info from", keyword, "\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
    if (keyword == "HA") {
        keyword = "H\\d+"
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

        jj <- getAmbiguous.index.base(gg)
        ii <- jj$index
        if (length(ii) == 0) {
            return("full")
        }
        res <- paste(ii, jj$base, sep="")
        if (length(ii) > 1) {
            res <- paste(res, sep=",", collapse = ",")
        }
        return(res)
    })
    return(unlist(result))
}

getAmbiguous.index.base <- function(consensusSeqFile) {
    x <- toupper(toString(readBStringSet(consensusSeqFile)))
    x <- substring(x, 1:nchar(x), 1:nchar(x))
    ii <- which(! x %in% c("A", "C", "G", "T"))
    if (length(ii) > 0) {
        ii <- ii[ii > 10 & ii < (length(x)-10)]
    }

    result <- list(index=ii, base=x[ii])
    return(result)
}


getSubtype <- function(nameMap, cs) {
    res <- sapply(1:nrow(nameMap), function(i) {
        strain <- cs[grep(nameMap[i,2], cs)]
        if (length(strain) == 0) {
            return("")
        }
        hh <- strain[grep("H\\d+[_mixed]*\\.fas$", strain)]
        nn <- strain[grep("N\\d[_mixed]*\\.fas$", strain)]
        ht <- gsub(".*(H\\d+)[_mixed]*\\.fas", "\\1", hh)
        nt <- gsub(".*(N\\d)[_mixed]*\\.fas", "\\1", nn)
        ht <- paste(ht, collapse = "|")
        nt <- paste(nt, collapse = "|")
        if (ht == "") {
            subtype <- nt
        } else if (nt == "") {
            subtype <- ht 
        } else {
            subtype <- paste(ht, nt, sep="|")
        }
        if (length(subtype) == 0)
            return("")
         return(subtype)
    })
    return(unlist(res))
}
