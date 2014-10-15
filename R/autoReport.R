##' auto report sequence alignment and consensus sequence
##'
##' 
##' @title autoReport
##' @param contig.folder contig folder
##' files in contig.folder should have name contains _S\\d+[HNMP]+, or it will be move to unknown folder
##' should contains for example: XXX_S111NA_454.fas, XXX_S111NA_Contigs.fna and XXX_S111NA_mira_delX.fasta
##' @param ref.folder reference folder
##' should contains files like XXX__S111NA_Ref.fas
##' @param name.file name file
##' @param out.folder output folder 
##' @param filter logical
##' @param percentage filter percentage
##' @return NULL
##' @importFrom markdown markdownToHTML
##' @export
##' @author ygc
autoReport <- function(contig.folder, ref.folder, name.file, out.folder="output", filter=FALSE, percentage=1) {
    printInfo()
    
    nameMap <- read.delim(name.file, header=F, stringsAsFactor=FALSE)

    outfile <- "report.md"
    sink(outfile)
    cat("# Sequence Report by `skleid`\n")
    cat("## Summary \n")
    sink()

    unKnownFileReport(contig.folder, outfile)
    emptyFileReport(contig.folder, ref.folder, outfile)
    NDVFileReport(contig.folder, outfile)
    
    if(filter == TRUE)
        doFilter(contig.folder, percentage)

    
    mixedFileReport(contig.folder, outfile)

    generateStrainTable(contig.folder, nameMap)

    contig <- getFiles(contig.folder)
    ref <- getFiles(ref.folder)
    ref <- ref[grep("Ref.fas", ref)]
    sink(outfile, append=TRUE)
    cat("## un-mixed Files\n")
    sink()
    processItems(outfile, contig, ref, nameMap, contig.folder, out.folder)

    mix.folder <- "mixed"
    if ( file.exists(mix.folder)) {
        if ( length(list.files(mix.folder)) > 0) {
            sink(outfile, append=TRUE)
            cat("\n\n")
            cat("## Mixed Files\n")
            sink()
            mix.contig <- getFiles(mix.folder)
            processItems(outfile, mix.contig, ref, nameMap, mix.folder, out.folder, outfile.suffix="mixed")        }    
    }
    
    addFooter(outfile)
    markdownToHTML(outfile, "report.html")
    file.remove(outfile)
    cat(">> done...", "\t\t\t", format(Sys.time(), "%Y-%m-%d %X"), "\n")
    printInfo()
}

processItems <- function(outfile, contig, ref, nameMap, contig.folder, out.folder, outfile.suffix="") {

    pc <- gsub(pattern=".*/.*_([SRL]+\\d+[HNMP][APSB]\\d*)_[4Cm].*", replacement="\\1", contig)
    ## if _Ref.fas in contig folder, it will move to missingFile folder in next step
    
    pr <- gsub(pattern=".*/.*_([SRL]+\\d+[HNMP][APSB]\\d*)_Ref.fas", replacement="\\1", ref)
    
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
    
    aln <- doAlign(seqs)
    
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

generateStrainTable <- function(contig.folder, nameMap) {
    mix.sc <- getMixedStrain()
    if (!is.null(mix.sc)) {
        mixff <- nameMap[ nameMap[,1] %in% sub("[SRL]+", "", unique(mix.sc)),]    
        mixff[,2] <- paste(mixff[,2], "mixed", sep="_")
        write.table(mixff,
                    file="mixed_table.txt", sep="\t",
                    row.names=F, col.names=F, quote=F)        
    }
        
    contig <- getFiles(contig.folder)
    sc <- getSampleID(contig)
    write.table(nameMap[nameMap[,1] %in% sub("S", "", unique(sc)), ],
                file="unmixed_table.txt", sep="\t",
                row.names=F, col.names=F, quote=F)        
}

doFilter <- function(contig.folder, percentage) {
    contig <- getFiles(contig.folder)
     
    f454 <- contig[grep("_454.fa[sta]$", contig)] ## contig[grep("454", contig)]
    cat(">>", "filtering sequence read numbers below", percentage, "percentage", "\n") 
    readNumFilter(f454, percentage)
}

NDVFileReport <- function(contig.folder, report.file) {
    moveNDVFile(contig.folder)

    if (file.exists("NDV")) {
        ff <- list.files(path="NDV")
        if (length(ff) > 0) {
            ndv <- gsub(".*/*_([SRL]+\\d+NDV)_.*",'\\1', ff)
            ndv <- unique(ndv)
            sink(report.file, append=TRUE)
            cat("\n", length(ndv), " [NDV](NDV) gene(s) found.\n")
            cat("\n\t", paste(ndv), "\n")
            sink()
        }
    }
}


mixedFileReport <- function(contig.folder, report.file) {
    mix.sc <- getMixedStrain()  
    if ( !is.null(mix.sc)) {
        sink(report.file, append=TRUE)
        cat("\n", length(unique(mix.sc)), " [mixed](mixed) strain(s) found.\n")
        cat("\n\t", paste(unique(mix.sc)), "\n")
        sink()
        cat(">>", length(unique(mix.sc)), " mixed strain(s) found.\n")
    }
    
    contig <- getFiles(contig.folder)
    sc <- getSampleID(contig)
    if (length(sc) > 0) {
        cat(">>", length(unique(sc)), " strains were processed in the report.\n")
        ## processing 
        sink(report.file, append=TRUE)
        cat(length(unique(sc)), " strains were processed in the report.\n")
        cat("\n\t", paste(unique(sc)), "\n\n")
        sink()
    }
}

emptyFileReport <- function(contig.folder, ref.folder, report.file) {
    moveEmptyFile(contig.folder)
    moveEmptyFile(ref.folder)

    if (file.exists("empty")) {
        num <- length(list.files(path="empty"))
        if (num > 0) {
            sink(report.file, append=TRUE)
            cat(num, " [empty](empty) file(s) found.\n")
            sink()
            cat(">>", num, " empty file(s) found.\n")
        }
    }
}


unKnownFileReport <- function(contig.folder, report.file) {
    ## move unknown files if any
    moveUnknownFile(contig.folder)
    if (file.exists("unknown")) {
        num <- length(list.files(path="unknown"))
        if (num > 0) {
            sink(report.file, append=TRUE)
            cat(num, " [unknown](unknown) file(s) found.\n")
            sink()
            cat(">>", num, " unknown file(s) found.\n")
        }
    }
}


