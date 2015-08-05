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
##' @param markAmbiguous makr ambiguous with * or not
##' @return NULL
##' @importFrom markdown markdownToHTML
##' @export
##' @author ygc
autoReport <- function(contig.folder, ref.folder, name.file, out.folder="output", filter=FALSE, percentage=1, markAmbiguous=TRUE) {
    
    nameMap <- read.delim(name.file, header=FALSE, stringsAsFactor=FALSE)

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
    processItems(outfile, contig, ref, nameMap,
                 contig.folder, out.folder,
                 markAmbiguous=markAmbiguous)

    mix.folder <- "mixed"
    if ( file.exists(mix.folder)) {
        if ( length(list.files(mix.folder)) > 0) {
            sink(outfile, append=TRUE)
            cat("\n\n")
            cat("## Mixed Files\n")
            sink()
            mix.contig <- getFiles(mix.folder)
            processItems(outfile, mix.contig, ref, nameMap, mix.folder,
                         out.folder, outfile.suffix="mixed",
                         markAmbiguous=markAmbiguous)
        }    
    }
    
    addFooter(outfile)
    markdownToHTML(outfile, "report.html")
    file.remove(outfile)
    ## cat(">> done...", "\t\t\t", format(Sys.time(), "%Y-%m-%d %X"), "\n")

    ## 
    cat(">> done... \n")
    cat("------------\n")
    says()
}



processItems <- function(outfile, contig, ref, nameMap, contig.folder,
                         out.folder, outfile.suffix="",
                         markAmbiguous) {

    contig <- contig[grep("[SRL]+\\d+[HNMP][APSB]\\d*_[4Cm].*", contig)]
    ## pc <- gsub(pattern=".*/.*_([SRL]+\\d+[HNMP][APSB]\\d*)_[4Cm].*", replacement="\\1", contig)
    pc <- get_sid_gn(contig)

    ## sort contig
    pc.segment <- gsub("^[SRL]+\\d+([HNMP][APSB]\\d*)", '\\1', pc)

    segment_order <- c("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")
    idx <- c()
    for (i in segment_order) {
        idx <- c(idx, which(pc.segment == i))
    }
    contig <- contig[idx]
    pc <- pc[idx]

    pc.strain <- gsub("^[SRL]+(\\d+)[HNMP][APSB]\\d*", '\\1', pc)
    pc.strain <- as.numeric(pc.strain)
    idx <- order(pc.strain, decreasing = FALSE)
    contig <- contig[idx]
    pc <- pc[idx]

    
    ## if _Ref.fas in contig folder, it will move to missingFile folder in next step
    ## pr <- gsub(pattern=".*/.*_([SRL]+\\d+[HNMP][APSB]\\d*)_Ref.fas", replacement="\\1", ref)
    ref <- ref[grep("_Ref.fas", ref)]
    pr <- get_sid_gn(ref)
    
    if (!file.exists(out.folder)) {
        dir.create(out.folder)
    }
    
    if (!file.exists("html")) {
        dir.create("html")
    }
    
    oldstrain <- ""
    for ( pp in unique(pc) ) {
        jj=contig[pc==pp]
        seqs <- c(jj[grep("454\\.fas", jj)], ## jj[grep("Contigs", jj)],
                  ref[pr == pp ],
                  jj[grep("454M\\.", jj)]
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
        
        ## seqname <- sub(contig.folder, "", seqs)
        ## seqname <- sub("/", "", seqname)
        ## seqname <- sub("\\w+/", "", seqs)
        seqname <- basename(seqs)
        
        if (isMixed(jj[grep("454\\.fas", jj)]) == TRUE) {
            pn <- pp
            outhtml <- NULL
            outfasta <- NULL
        } else {
            outinfo <- itemReport(seqs, seqname, pp, nameMap, out.folder, outfile.suffix)
            pn <- outinfo$pn
            outhtml <- outinfo$outhtml
            outfasta <- outinfo$outfasta
        }

        
        strain <- gsub("([SRL+]\\d+).*", "\\1", pp)

        sink(outfile, append=TRUE)
        if ( strain != oldstrain) {
            cat("\n")
            if (is.null(outhtml)) {
                cat("### ", pn)
            } else {
                cat("### ")
            }
        } else {
            if (is.null(outhtml)) {
                cat("\t", pn)
            } else {
                cat("\t")
            }
        }

        if (!is.null(outfasta) && markAmbiguous == TRUE){
            ii <- getAmbiguous.index.base(outfasta)$index
            if (length(ii) > 0) {
                pn <- paste0(pn, "*")
            }
        }

        if (!is.null(outhtml)) {
            cat(paste("[", pn, "]", "(", outhtml, ")", sep=""))
        }
        
        sink()
        oldstrain <- strain
    }
    gc()
}


itemReport <- function(seqs, seqname, pp, nameMap, out.folder, outfile.suffix) {
    ## seqs: sequence file name
    ## seqname: sequence name
    ## pp: protein name
    
    cat(">>", "processing", pp, "\t\t", format(Sys.time(), "%Y-%m-%d %X"), "\n",
        "\t\t", seqname[1], "\n",
        "\t\t", seqname[2], "\n",
        "\t\t", seqname[3], "\n")

    ### sometimes in Windows, it may throw error for no permission to read temp.rafa file.
    ### but the alignment was no problem.
    ### when this happened, re-do the alignment is always (almost) works.
    aln <- tryCatch(doAlign(seqs), error=function(e) NULL)

    if (is.null(aln)) {
        Sys.sleep(3)
        aln <- tryCatch(doAlign(seqs), error=function(e) NULL)
        if (is.null(aln)) {
            Sys.sleep(3)
            aln <- doAlign(seqs)
        }
    }
    
    if (length(grep("HA", pp)) > 0) {
        ## pn <- gsub(".*(H\\d+)N\\d+.*", "\\1", aln$seqs[1,1])
        pn <- gsub(".*(H\\d+)N\\d+.*", "\\1", names(aln)[1])
        pn <- sub("HA", pn, pp)
    } else if (length(grep("NA", pp)) > 0) {
        ## pn <- gsub(".*H\\d+(N\\d+).*", "\\1", aln$seqs[1,1])
        pn <- gsub(".*H\\d+(N\\d+).*", "\\1", names(aln)[1])
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

    alnDir <- "alignment"
    alnFile <- paste(pn, "_alignment.fas", sep="")
    alnFile2 <- paste(alnDir, alnFile, sep="/")
    if(!file.exists(alnDir))
        dir.create(alnDir)
    writeAlignedSeq(aln, alnFile2)
    
    cat("```\n")
    outfasta <- nameMap[ nameMap[,1] == gsub("\\D+(\\d+)\\w+", "\\1", pp), 2]
    outfasta <- paste(as.character(outfasta), sub("[SRL]+\\d+", "", pn), sep="_")
    if (outfile.suffix != "") {
        outfasta <- paste(outfasta, outfile.suffix, sep="_")
    }
    cat("\n\nAlignment file", paste("[", alnFile, "]", sep = ""))
    cat(paste("(", "../", alnFile2, ")", sep = ""), "generated.\n")
    
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

    res <- list(pn=pn, outhtml=outhtml, outfasta=outfasta)
    return(res)
}

generateStrainTable <- function(contig.folder, nameMap) {
    mix.sc <- getMixedStrain()
    if (!is.null(mix.sc)) {
        ## mixff <- nameMap[ nameMap[,1] %in% sub("[SRL]+", "", unique(mix.sc)),]
        ii <- nameMap[,1] %in% gsub("[SRL]+(\\d+).*", "\\1", unique(mix.sc))
        if (sum(ii) == 0) {
            cat("--> mixed samples are not listed in name.txt...\n")
            cat("--> press ENTER to skip...\n")
            pause()
        } else {
            mixff <- nameMap[ii,]    
            mixff[,2] <- paste(mixff[,2], "mixed", sep="_")
            write.table(mixff,
                        file="mixed_table.txt", sep="\t",
                        row.names=F, col.names=F, quote=F)        
        }
    }
        
    contig <- getFiles(contig.folder)
    sc <- getSampleID(contig)
    write.table(nameMap[nameMap[,1] %in% gsub("[SRL]+(\\d+).*", "\\1", unique(sc)),],
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
            ndv <- gsub(".*/.*([SRL]+\\d+NDV)_.*",'\\1', ff)
            ndv <- unique(ndv)
            sink(report.file, append=TRUE)
            cat("\n", length(ndv), " [NDV](NDV) gene(s) found.\n")
            cat("\n\t", paste(ndv), "\n")
            sink()
        }
    }
}


mixedFileReport <- function(contig.folder, report.file) {
    moveMixedFile(contig.folder)
    
    mix.sc <- getMixedStrain()
    sampleID <- gsub("^[SRL](\\d+)$", '\\1', mix.sc)
    sampleID <- as.numeric(sampleID)
    idx <- order(sampleID, decreasing = FALSE)
    mix.sc <- mix.sc[idx]

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


