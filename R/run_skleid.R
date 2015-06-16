##' launch shiny app
##'
##' 
##' @title run_skleid
##' @return NULL
##' @importFrom shiny runApp
##' @import shinyFiles
##' @export
##' @author Guangchuang Yu
run_skleid <- function() {
    require("shinyFiles", character.only=TRUE)
    dir <- system.file("app", package="skleid")
    options(shiny.launch.browser = TRUE)
    runApp(dir)
}


##' generate desktop shortcut
##'
##' 
##' @title generate_run_skleid_app
##' @return NULL
##' @export
##' @author Guangchuang Yu
generate_run_skleid_app <- function() {
    os <- Sys.info()[["sysname"]]
    if (os != "Windows") {
        stop("This function only works in Windows Platform...")
    }
    
    outfile <- paste0(file.path(Sys.getenv("USERPROFILE"), "Desktop"), "/run_skleid.bat")
    out <- file(outfile, "w")
    
    writeLines(sub("library/base", "bin/Rscript -e 'skleid::run_skleid()'", system.file(package="base")), out)
    close(out)
}
