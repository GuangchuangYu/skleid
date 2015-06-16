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
    runApp(dir)
}
