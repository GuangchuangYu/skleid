##' generating gap file
##'
##' 
##' @title generateGapFile
##' @param out.folder output folder of autoReport, that contains consensus sequences,
##'                   or summary.csv from summarize output
##' @param ref.folder reference folder
##' @param read.fileName list of file names of read file
##' @param sample.size sample size
##' @param gap.file output gap file name 
##' @return NULL 
##' @author ygc
##' @importFrom Biostrings readBStringSet
##' @importFrom magrittr %>%
##' @export
generateGapFile <- function(out.folder="output", ref.folder="Ref", read.fileName, sample.size=200, gap.file="gaps.txt") {
    ## ref: AH1-P5_S5HA_Ref.fas
    ## read: "HA_CK_JX_22232_2014_S1HA.fa"
    ## consensus seq: AB53TL10F06_DK_JX_13730_2014_H4.fas , _mixed.fas for mixed 
    ##
    ## consensus seq name should contain virus strain that match read name
    ## read name should contain strain name and gene name, that can be matched by ref
    ##
    ##
    reads <- read.delim(read.fileName, header=F)
    reads <- as.character(unlist(reads))
    reads <- unique(reads)

    ref <- list.files(path=ref.folder, pattern=".fas$")

    if (file.info(out.folder)$isdir) {
        cs <- list.files(path=out.folder, pattern=".fas$")
        cs2 <- paste(out.folder, cs, sep="/")
    } else {
        sf <- read.csv(out.folder)
    }

    writeGap_ <- function(ss, reads, ref, pp, ii, gapfile, name) {
        
        if (pp == "NA.") {
            pp <- "NA"
        }
        rr <- reads[grep(ss, reads)]
        rr <- rr[grep(paste(".*[SRL]+\\d+(", pp, ")\\..*", sep=""), rr)]

        if (any(sapply(rr, function(a) grepl("RL\\d+", a)))) {
            ## in this case ss should be RLxx
            rr <- rr[grep(paste0(ss, pp), rr)]
        }
        
        if (length(rr) == 0) {
            warning("gene ", pp, " of ", ss, " not found...")
            return(NULL)
        } 
        
        sg <- get_sid_gn(rr)
        strain <- gsub("([SRL]+\\d+)[HNMP][APSB]*\\d*.*", "\\1", sg)
        strain <- paste(strain, pp, sep="") %>% unique
        
        ref2 <- ref[grep(strain, ref)]
        if (length(ref2) == 0) {
            warning("reference of ", strain, " not found...")
            return(NULL)
        }
        ff <- paste("File:", rr, sep="")
        fr <- paste("Ref:", ref2, sep="")

        idx <- grep("-", ii)
        pre <- numeric()
        post <- numeric()
        if (length(idx) > 0) {
            for (kk in idx) {
                xx <- strsplit(ii[kk], "-") %>% unlist %>% as.numeric
                pre <- c(pre, xx[1])
                post <- c(post, xx[2])
            }
            ii <- ii[-idx]
        }
        ii <- as.numeric(ii)
        pre <- c(pre, ii-10)
        post <- c(post, ii)
        jj <- paste(pre, post, sep=",")
        cat("-> ", length(jj), "gap(s) found in", name, "\n")
        if (length(jj) > 1) {
            fg <- paste(jj, collapse=" | ")
        } else {
            fg <- jj
        }
        fg <- paste("Gaps:", fg, sep="")
        fn <- paste("#:", sample.size, "\n", sep="")
        writeLines(ff, gapfile)
        writeLines(fr, gapfile)
        writeLines(fg, gapfile)
        writeLines(fn, gapfile)
    }

    gapfile <- file(gap.file, "w")
    if (file.info(out.folder)$isdir == FALSE) {
        for (i in 1:nrow(sf)) {
            for (j in c("HA", "NA.", "PB2", "PB1", "PA", "NP", "MP", "NS")) {
                entry <- sf[i, j]
                if (is.null(entry)) {
                    next
                }
                entry <- as.character(entry)
                if (entry != "" && entry != "full") {
                    x <- strsplit(entry, "\\|") %>% unlist %>% unique
                    if (x == "full") {
                        next
                    }
                    ss <- sub("^[a-zA-Z0-9]+_", "", sf[i, "name"])
                    ii <- strsplit(entry, ",") %>% unlist %>%
                        strsplit(split="\\|") %>% unlist
                    ii <- ii[ii!="full"]
                    if (length(ii) == 0) {
                        next
                    }
                    ii <- ii[ii!="mix"]
                    if (length(ii) == 0) {
                        next
                    }
                    ii <- sub("[A-Za-z]$", "", ii) 
                    nn <- paste0(as.character(sf[i, "name"]), "_", j)
                    notuse <- writeGap_(ss, reads, ref, j, ii, gapfile, nn)
                }
            }
        }
    } else {
        for (i in 1:length(cs)) {
            jj <- getAmbiguous.index.base(cs2[i])
            ii <- jj$index
            ## x <- toupper(toString(readBStringSet(cs2[i])))
            ## x <- substring(x, 1:nchar(x), 1:nchar(x))
            ## ii <- which(! x %in% c("A", "C", "G", "T"))
            ## if (length(ii) > 0) {
            ##     ii <- ii[ii > 10 && ii < (length(x)-10)]
            ## }
            if (length(ii) == 0) {
                next
            }

     
            ## remove prefix
            ss <- sub("^[a-zA-Z0-9]+_", "", cs[i])
            ## remove _mixed if any
            ss <- sub("_mixed", "", ss)
            ## remove gene name and .fas
            ss <- sub("_[A-Z0-9]+.fas$", "", ss)
            ## pp <- gsub(".*_(\\w+).fas", '\\1', cs[i])
            pp <- gsub(".*_([HNMP][APSB]*\\d*)[_mixed]*\\.fas", '\\1', cs[i])
            if (pp == "M") {
                pp <- "MP"
            } else if (length(grep("H", pp)) > 0) {
                pp <- "HA"
            } else if ( length(grep("N\\d", pp)) > 0 ) {
                ## pp != "NS" && pp != "NP" && length(grep("N", pp)) > 0) {
                pp <- "NA"
            }
            notuse <- writeGap_(ss, reads, ref, pp, ii, gapfile, cs[i])
        }
    }
      
    close(gapfile)
    cat("-> done... \n")
    cat("------------\n")
    says()
}

