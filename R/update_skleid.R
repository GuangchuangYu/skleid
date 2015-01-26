##' auto update skleid package if needed
##'
##' 
##' @title update_skleid
##' @return NULL
##' @importFrom RCurl getURL
##' @importFrom utils packageVersion
##' @importFrom devtools install_github
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
        flag <- FALSE
        ## source("http://ygc.name/get_skleid.R")
        for (lp in .libPaths()) {
            libs <- list.files(lp)
            if (length(grep("skleid", libs)) == 1) {
                source(paste(lp, "skleid/extdata/get_skleid.R", sep="/"))
                flag <- TRUE
                break
            }
        }
        if (flag == FALSE) {
            source("http://ygc.name/get_skleid.R")
        }
    }

    y <- getURL("https://raw.githubusercontent.com/GuangchuangYu/ggtree/master/DESCRIPTION",
                .opts = list(ssl.verifypeer = FALSE))
    v <- gsub(".*\nVersion: (\\d+\\.\\d+\\.\\d+)\n.*", "\\1", y)
    if (as.character(packageVersion("ggtree")) != v) {
        cat("-> new version (", v, ") of ggtree is available\n")
        cat("-> press ENTER to update the package...\n")
        pause()
        install_github("GuangchuangYu/ggtree")
    }
         
    library("skleid")
    
    cat("  __________________________________________\n")
    cat("/ skleid (version=", vv, ") is up to date... \\\n") 
    cat("\\ Have fun with SKLEID...                   /\n")
    cat("  ------------------------------------------\n")
    cheeseSay()
}
