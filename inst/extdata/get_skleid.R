is.installed <- function(pkg) {
    pkgs <- installed.packages()[,1]
    if (length(which(pkg == pkgs))==0) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

installPKG <- function(pkg) {
    if (!is.installed(pkg)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkg)
    }
}

installPKG("Biostrings")
installPKG("devtools")

library("devtools")
if (!is.installed("ggtree")) {
    install_github("GuangchuangYu/ggtree")
}

install_github("GuangchuangYu/skleid")

## suppressPackageStartupMessages(library("skleid"))
## printInfo()
