##' auto report sequence alignment and consensus sequence
##'
##' 
##' @title autoReport
##' @param contig.folder contig folder
##' @param ref.folder reference folder
##' @param name.file name file
##' @param out.folder output folder 
##' @param filter logical
##' @param percentage filter percentage
##' @return NULL
##' @importFrom markdown markdownToHTML
##' @export
##' @author ygc
autoReport <- function(contig.folder, ref.folder, name.file, out.folder="output", filter=FALSE, percentage=1) {
    contig <- list.files(path=contig.folder)
    contig <- paste(contig.folder, contig, sep="/")

    nameMap <- read.delim(name.file, header=F, stringsAsFactor=FALSE)
    outfile <- "report.md"
    sink(outfile)
    cat("# Sequence Report by `skleid`\n")
    
    cat("## Summary \n")
    
    sink()

    sc <- gsub(pattern=".*/[^SRL]*([SRL]+\\d+).*", replacement="\\1", contig)
    
    ## move unknown files if any
    sink(outfile, append=TRUE)
    idx <- which(sc == contig)
    if (length(idx) >= 1) {
        ## print("found unknown files...")
        if (!file.exists("unknown"))
            dir.create("unknown")
        for (i in idx) {
            file.rename(contig[i], to=sub(contig.folder,"unknown", contig[i]))
        }
        contig <- contig[-idx]
        sc <- sc[-idx]
        cat(length(idx), " [unknown](unknown) file(s) found.\n")
    }
    
    sink()

    if (length(idx) >= 1)
        cat(">>", length(idx), " unknown file(s) found.\n")
    
    ## filter
    f454 <- contig[grep("454", contig)]
    if (filter == TRUE) {
        cat(">>", "filtering sequence read numbers below", percentage, "percentage", "\n") 
        readNumFilter(f454, percentage)
    }

    ## move mixed files
    idx <- getMixedFileIndex(f454)
    mix <- unique(gsub(pattern=".*/[^SRL]*([SRL]+\\d+).*", replacement="\\1", f454[idx]))
    
    
    if (length(mix) >= 1) {
        sink(outfile, append=TRUE)
        ## print("mixed files found...")
        if (!file.exists("mixed"))
            dir.create("mixed")
        ii <- which(sc %in% mix)
        for (i in ii) {
            file.rename(contig[i], to=sub(contig.folder, "mixed", contig[i]))
        }
        cat("\n", length(unique(sc[ii])), " [mixed](mixed) strain(s) found.\n")
        ms <- unique(sc[ii])
        cat("\n\t", paste(ms), "\n")
        sink()
        
        cat(">>", length(unique(sc[ii])), " mixed strain(s) found.\n")
        mixff <- nameMap[ nameMap[,1] %in% sub("[SRL]+", "", ms),]
        mixff[,2] <- paste(mixff[,2], "mixed", sep="_")
        write.table(mixff,
                    file="mixed_table.txt", sep="\t",
                    row.names=F, col.names=F, quote=F)        
        contig <- contig[-ii]
        sc <- sc[-ii]
        write.table(nameMap[! nameMap[,1] %in% sub("S", "", ms), ],
                    file="unmixed_table.txt", sep="\t",
                    row.names=F, col.names=F, quote=F)        
    }

    cat(">>", length(unique(sc)), " strains were processed in the following report.\n")
    ## processing 
    sink(outfile, append=TRUE)
    cat(length(unique(sc)), " strains were processed in the following report.\n")
    cat("\n\t", paste(unique(sc)), "\n\n")
    cat("## un-mixed Files\n")
    sink()

    
    ref <- list.files(path=ref.folder)
    ref <- paste(ref.folder, ref, sep="/")
    ref <- ref[grep("Ref.fas", ref)]

    
    processItems(outfile, contig, ref, nameMap, contig.folder, out.folder)

    mix.folder <- "mixed"
    if ( file.exists(mix.folder)) {
        if ( length(list.files(mix.folder)) > 0) {
            sink(outfile, append=TRUE)
            cat("\n\n")
            cat("## Mixed Files\n")
            sink()
    
    
            mix.contig <- list.files(path=mix.folder)
            mix.contig <- paste(mix.folder, mix.contig, sep="/")
            processItems(outfile, mix.contig, ref, nameMap, mix.folder, out.folder, outfile.suffix="mixed")
            
            addFooter(outfile)
            markdownToHTML(outfile, "report.html")
            file.remove(outfile)
            file.remove("tmp.log")
            cat(">> done...", "\t\t\t", format(Sys.time(), "%Y-%m-%d %X"), "\n")
        }
    }

}

processItems <- function(outfile, contig, ref, nameMap, contig.folder, out.folder, outfile.suffix="") {

    pc <- gsub(pattern=".*/[^SRL]*([SRL]+\\d+[a-zA-Z]+\\d*)_[4Cm].*", replacement="\\1", contig)
    ## if _Ref.fas in contig folder, it will move to missingFile folder in next step
    
    pr <- gsub(pattern=".*/[^SRL]*([SRL]+\\d+\\w+\\d*)_Ref.fas", replacement="\\1", ref)
    
    if (!file.exists(out.folder)) {
        dir.create(out.folder)
    }
    
    if (!file.exists("html")) {
        dir.create("html")
    }
    
    oldstrain <- ""
    for ( pp in unique(pc) ) {
        jj=contig[pc==pp]
        seqs <- c(jj[grep("Contigs", jj)],
                  ref[pr == pp ],
                  jj[grep("mira_delX", jj)]
                  )
        if (length(seqs) != 3) {
            if (!file.exists("missingFile")) {
                dir.create("missingFile")
            }
            for (i in jj) {
                file.rename(i, to=sub(contig.folder, "missingFile", i))
            }
            next
        }
        
        seqname <- sub(contig.folder, "", seqs)
        seqname <- sub("/", "", seqname)
        

        if (isMixed(jj[grep("454.fas", jj)]) == TRUE) {
            pn <- pp
            outhtml <- NULL
        } else {
            outinfo <- itemReport(seqs, seqname, pp, nameMap, out.folder, outfile.suffix)
            pn <- outinfo$pn
            outhtml <- outinfo$outhtml
        }

        
        strain <- gsub("([SRL+]\\d+).*", "\\1", pp)

        sink(outfile, append=TRUE)
        if ( strain != oldstrain) {
            cat("\n")
            if (is.null(outhtml)) {
                cat("### ", pn)
            } else {
                cat("### ", paste("[", pn, "]", "(", outhtml, ")", sep=""))
            }
        } else {
            if (is.null(outhtml)) {
                cat("\t", pn)
            } else {
                cat("\t", paste("[", pn, "]", "(", outhtml, ")", sep=""))
            }
        }
        sink()
        oldstrain <- strain
    }
}


itemReport <- function(seqs, seqname, pp, nameMap, out.folder, outfile.suffix) {
    ## seqs: sequence file name
    ## seqname: sequence name
    ## pp: protein name
    
    cat(">>", "processing", pp, "\t\t", format(Sys.time(), "%Y-%m-%d %X"), "\n",
        "\t\t", seqname[1], "\n",
        "\t\t", seqname[2], "\n",
        "\t\t", seqname[3], "\n")
    
    
    sink("tmp.log")
    aln <- doAlign(seqs)
    sink()
    
    if (length(grep("HA", pp)) > 0) {
        pn <- gsub(".*(H\\d+)N\\d+.*", "\\1", aln$seqs[1,1])
        pn <- sub("HA", pn, pp)
    } else if (length(grep("NA", pp)) > 0) {
        pn <- gsub(".*H\\d+(N\\d+).*", "\\1", aln$seqs[1,1])
        pn <- sub("NA", pn, pp)
    } else {
        pn <- pp
    }
    
    
    outmd <- paste("html/", pn, ".md", sep="")
    sink(outmd)
    cat("## ", pn, "\n")
    cat(paste("\n", "[", seqname, "]", "(", "../", seqs, ")", sep="", collapse="\n"), "\n")
        
    cat("\n\naligned sequences:\n```\n")
    printAlignedSeq(aln)        
    printConsensus(aln)
    cat("```\n")
    outfasta <- nameMap[ nameMap[,1] == gsub("\\D+(\\d+)\\w+", "\\1", pp), 2]
    outfasta <- paste(as.character(outfasta), sub("[SRL]+\\d+", "", pn), sep="_")
    if (outfile.suffix != "") {
        outfasta <- paste(outfasta, outfile.suffix, sep="_")
    }
    outfasta <- paste(outfasta, ".fas", sep="")
    cat("\n\nConsensus file", paste("[", outfasta, "]", sep=""))
    outfasta <- paste(out.folder, outfasta, sep="/")
    cat(paste("(", "../", outfasta, ")", sep=""), "generated.\n")
    sink()
    writeConsensus(aln, output=outfasta)
    outhtml <- sub(".md", ".html", outmd)
    addFooter(outmd)
    markdownToHTML(outmd, outhtml)
    file.remove(outmd)

    res <- list(pn=pn, outhtml=outhtml)
    return(res)
}



getMixedFileIndex <- function(files) {
    which(sapply(files, isMixed))
}


isMixed <- function(file) {
    prot <- gsub(pattern=".*/[^SRL]*[SRL]+\\d+(\\w+\\d*)_454.*", replacement="\\1", file)
    
    prot.cutoff <- c(rep(2000, 4), 1000, rep(3000, 3))
    names(prot.cutoff) <- c("HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2")
    
    res <- file.info(file)$size > prot.cutoff[prot]
    ## NDV that is not exists
    if (is.na(res))
        return(TRUE)
    return(res)
}
    
