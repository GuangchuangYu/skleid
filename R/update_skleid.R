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
        cat("new version (", vv, ") of skleid found...\n")
        ## cat("press ENTER to update the package...\n")
        pause()
        source("http://ygc.name/get_skleid.R")
    } else {
        cat("skleid package (version=", vv, ") is up to date...\n")
    }

    cat("\n\nHave fun with SKLEID...\n")
    cat("___________________________\n")
    cheeseSay()
}
