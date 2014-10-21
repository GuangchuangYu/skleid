##' generating gap file
##'
##' 
##' @title generateGapFile
##' @param out.folder output folder of autoReport, that contains consensus sequences
##' @param ref.folder reference folder
##' @param read.fileName list of file names of read file
##' @param sample.size sample size
##' @param gap.file output gap file name 
##' @return NULL 
##' @author ygc
##' @export
generateGapFile <- function(out.folder="output", ref.folder="Ref", read.fileName, sample.size=200, gap.file="gaps.txt") {
    ## ref: AH1-P5_S5HA_Ref.fas
    ## read: "HA_CK_JX_22232_2014_S1HA.fa"
    ## consensus seq: AB53TL10F06_DK_JX_13730_2014_H4.fas
    ##
    ## consensus seq name should contain virus strain that match read name
    ## read name should contain strain name and gene name, that can be matched by ref
    ##
    ##
    reads <- read.delim(read.fileName, header=F)
    reads <- as.character(unlist(reads))

    ref <- list.files(path=ref.folder, pattern=".fas$")
    
    cs <- list.files(path=out.folder, pattern=".fas$")
    cs2 <- paste(out.folder, cs, sep="/") 

    gapfile <- file(gap.file, "w")
    for (i in 1:length(cs)) {
        x <- toupper(toString(readDNAStringSet(cs2[i])))
        x <- substring(x, 1:nchar(x), 1:nchar(x))
        ii <- which(! x %in% c("A", "C", "G", "T"))
        if (length(ii) > 0) {
            ii <- ii[ii > 10 && ii < (length(x)-10)]
        }
        if (length(ii) == 0) {
            next
        }
        ss <- sub("_[A-Z0-9]+.fas$", "", sub("^[a-zA-Z0-9]+_", "", cs[i]))
        rr <- reads[grep(ss, reads)]
        if (length(rr) == 0) {
            warning("file name, ", ss, " not match...")
            next
        }
        pp <- gsub(".*_(\\w+).fas", '\\1', cs[i])
        if (length(grep("H", pp)) > 0) {
            pp <- "HA"
        } else if (length(grep("N", pp)) > 0) {
            pp <- "NA"
        }
        ## rr <- rr[grep(pp, rr)]
        rr <- rr[grep(paste(".*_[SRL]+\\d+(", pp, ")\\..*", sep=""), rr)]
        
        if (length(rr) == 0) {
            warning("gene ", pp, " of ", ss, " not found...")
            next
        } 
        strain <- gsub(".*_([SRL]+\\d+)[HNMP][APSB]\\d*\\..*", "\\1", rr)
        strain <- paste(strain, pp, sep="")
        ref2 <- ref[grep(strain, ref)]
        if (length(ref2) == 0) {
            warning("reference of ", strain, " not found...")
            next
        } 
        ff <- paste("File:", rr, sep="")
        fr <- paste("Ref:", ref2, sep="")
        jj <- paste(ii - 10, ii, sep=",")
        if (length(jj) > 1) {
            fg <- paste(jj, collapse=" | ")
        } else {
            fg <- jj
        }
        fn <- paste("#:", sample.size, "\n", sep="")
        writeLines(ff, gapfile)
        writeLines(fr, gapfile)
        writeLines(fg, gapfile)
        writeLines(fn, gapfile)
    }
      
    close(gapfile)
}