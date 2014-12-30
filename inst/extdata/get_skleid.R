installPKG <- function(pkg) {
    pkgs <- installed.packages()[,1]
    if (length(which(pkg == pkgs))==0) {
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg)
    }
}

installPKG("Biostrings")
installPKG("devtools")

library("devtools")
install_github("GuangchuangYu/ggtree")
install_github("GuangchuangYu/skleid")

## suppressPackageStartupMessages(library("skleid"))
## printInfo()
