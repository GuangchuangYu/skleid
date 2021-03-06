##' auto update skleid package if needed
##'
##' 
##' @title update_skleid
##' @return NULL
##' @importFrom RCurl getURL
##' @importFrom utils packageVersion
##' @export
##' @author ygc
update_skleid <- function() {
    x <- getURL("https://raw.githubusercontent.com/GuangchuangYu/skleid/master/DESCRIPTION",
                .opts = list(ssl.verifypeer = FALSE))
    vv <- gsub(".*\nVersion: (\\d+\\.\\d+\\.\\d+)\n.*", "\\1", x)
    if (as.character(packageVersion("skleid")) != vv) {
        cat("-> new version (", vv, ") of skleid is available\n")
        cat("-> press ENTER to update the package...\n")
        pause()
        detach("package:skleid", character.only=TRUE)

        ## source("http://ygc.name/get_skleid.R")
        source_github("https://raw.githubusercontent.com/GuangchuangYu/skleid/master/inst/extdata/get_skleid.R")
        
        ## ## 
        ## flag <- FALSE
        ## for (lp in .libPaths()) {
        ##     libs <- list.files(lp)
        ##     if (length(grep("skleid", libs)) == 1) {
        ##         source(paste(lp, "skleid/extdata/get_skleid.R", sep="/"))
        ##         flag <- TRUE
        ##         break
        ##     }
        ## }
        ## if (flag == FALSE) {
        ##     source("http://ygc.name/get_skleid.R")
        ## }
    }
         
    library("skleid")
    
    cat("  __________________________________________\n")
    cat("/ skleid (version=", vv, ") is up to date... \\\n") 
    cat("\\ Have fun with SKLEID...                   /\n")
    cat("  ------------------------------------------\n")
    cheeseSay()
}


source_github <- function(url, global=TRUE) {
    script <- getURL(url, .opts = list(ssl.verifypeer = FALSE))
    if (global) {
        eval(parse(text = script), envir=.GlobalEnv)
    } else {
        eval(parse(text = script))
    }
}
