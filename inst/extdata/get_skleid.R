installPKG <- function(pkg) {
    pkgs <- installed.packages()[,1]
    if (length(which(pkg == pkgs))==0) {
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg)
    }
}

installPKG("Biostrings")
installPKG("devtools")

require("devtools")
install_github("GuangchuangYu/skleid")

suppressPackageStartupMessages(library("skleid"))
printInfo()
