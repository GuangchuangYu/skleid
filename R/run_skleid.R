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
    dir <- system.file(package="skleid")
    runApp(dir)
}
