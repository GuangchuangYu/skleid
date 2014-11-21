##' merge files and generate annotation file for treeAnno
##'
##' 
##' currently support merging output of NetNGlyc server (http://www.cbs.dtu.dk/services/NetNGlyc/)
##' @title fileMerge
##' @param files input files
##' @param pattern input file pattern
##' @param output output file
##' @param removeColnames rows to remove
##' @param removeRownames columns to remove
##' @importFrom reshape2 melt
##' @export
##' @return data.frame
##' @author Guangchuang Yu
fileMerge <- function(files=NULL, pattern="txt$",
                      output,
                      removeColnames="Frequency",
                      removeRownames="Frequency") {
    if ( is.null(files) ) {
        files <- list.files(pattern=pattern)
    }
    
    data <- lapply(files, read.delim, row.names=1)
    if (!is.null(removeColnames) || !is.na(removeColnames)) {
        data <- lapply(data, function(d) {
            i <- sapply(removeColnames, grep, x=colnames(d))
            if (length(unlist(i)) != 0) {
                d[, -i]
            } else {
                d
            }
        })
    }
    
    if (!is.null(removeRownames) || !is.na(removeRownames)) {
        data <- lapply(data, function(d) {
            i <- sapply(removeRownames, grep, x=rownames(d))
            if (length(unlist(i)) != 0) {
                d[-i,]
            } else {
                d
            }
        })
    }

    allcn <- unique(unlist(lapply(data, colnames)))
    for (i in 1:length(data)) {
	cn = allcn[!allcn %in% colnames(data[[i]])] 	
	for (j in cn) {
            data[[i]][j] <- ""
	}
    }

    for (i in 1:length(data)) {
	idx=order(colnames(data[[i]]))
	data[[i]] <- data[[i]][,idx] 
    }
    
    res <- do.call("rbind", data)
    write.table(res, output, sep="\t")

    res$name <- rownames(res)
    res2 <- suppressWarnings(melt(res, id.vars=c("name")))
    res2 <- res2[res2$value != "",]
    res2 <- res2[res2$value != " ",]
    res2$variable <- as.character(res2$variable)
    res2$variable <- as.numeric(sub("X", "", res2$variable))
    colnames(res) <- c("name", "position", "seq")
    invisible(res2)
}


##' extract N-glyco sites from output of NetNGlyc 1.0
##'
##' 
##' @title extract_netNglyc
##' @param infile output of NetNGlyc
##' @param outfile output file
##' @param cutoff cutoff
##' @param simplify boolean to remove null column
##' @export
##' @return data.frame
##' @author Guangchuang Yu
extract_netNglyc <- function(infile, outfile, cutoff=0.5, simplify=TRUE) {
    x <- extractNglycLine(infile)
    info <- extractInfo(x, cutoff)
    res <- convertToLongFormat(info, simplify)
    write.table(res, file=outfile, sep="\t", quote=F, row.names=F)
    invisible(info)
}

extractNglycLine <- function(file) {
    ## input is the output of netNglyc
    x <- readLines(file)
    idx <- grep(pattern="[A-Za-z0-9[:punct:]]+\\s+\\d+\\s+\\w+\\s+[\\(\\/\\)0-9]+",x)
    return(x[idx])
}

extractLineInfo <- function(y) {
    name <- gsub("([A-Za-z0-9[:punct:]]+)\\s+.*","\\1", y, perl=T)
    pos <- gsub("[A-Za-z0-9[:punct:]]+\\s+(\\d+).*","\\1", y, perl=T)
    seq <- gsub("[A-Za-z0-9[:punct:]]+\\s+\\d+\\s+(\\w+).*","\\1", y, perl=T)
    score <- gsub("[A-Za-z0-9[:punct:]]+\\s+\\d+\\s+\\w+\\s+([\\d.]+).*","\\1", y, perl=T)
    
    return(list(name=name, position=as.numeric(pos),
                seq=seq, score=as.numeric(score)))
}

extractInfo <- function(x, cutoff) {
    yy <- lapply(x, extractLineInfo)

    info <- do.call(rbind, yy)
    info <- as.data.frame(info)
    info$name <- unlist(info$name)
    info$position <- as.numeric(info$position)
    info$score <- as.numeric(info$score)
    info <- info[info$score > cutoff, ]
    return(info)
}

convertToLongFormat <- function(info, simplify=TRUE) {
    res <- matrix(NA, nrow=length(unique(info$name)), ncol=max(info$pos)+2)
    
    res <- as.data.frame(res)
    colnames(res) <- c("Position", as.character(1:max(info$pos)), "Frequency")
    nn <- unique(info$name)
    
    for (i in 1:length(nn)) {
        res[i,1] = nn[i]
    }
    
    for (i in 1:nrow(info)) {
        idx <- which(info$name[i] == res[, 1])
        res[idx, info$pos[i]+1] = info$seq[i]
    }
    res$Frequency <- apply(res, 1, function(i) sum(!is.na(i))-1)

    
    if (simplify == TRUE) {
        j <- apply(res, 2, function(i) sum(!is.na(i)))
        res <- res[, j!=0]
    }
    
    for (i in 1:nrow(res)) {
        for (j in 1:ncol(res)) {
            if (is.na(res[i,j])) {
                res[i,j]=""
            }
        }
        
    }

    return (res)
}

##' plot functional site annotation
##'
##' 
##' @title plotAnno
##' @param anno1 annotation data.frame, require
##' @param anno2 annotation data.frame, optional
##' @param xlab xlab
##' @param ylab ylab
##' @param title title
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##' @export
##' @return ggplot2 object
##' @author Guangchuang
plotAnno <- function(anno1, anno2=NULL, xlab="Position", ylab="Potential", title="") {
    score_flag = TRUE
    if ("score" %in% names(anno1) == FALSE) {
        ## if anno was obtained from NetNGlyc 1.0 Server
        ## score is missing.
        anno1$score <- 0.5
        score_flag <- FALSE
    }
    position <- score <- name <- NULL
    p <- ggplot(anno1, aes(x=position, y=score, color = name))
    p <- p+geom_segment(aes(x=position, xend=position, y=0, yend=score))
    if (score_flag == TRUE) 
        p <- p+geom_hline(y=0.5, linetype="longdash", color="steelblue")

    
    if (! is.null(anno2)) {
        if (score_flag == FALSE)
            anno2$score <- 0.5
        pp <- p + geom_segment(data=anno2, aes(x=position, xend=position,
                                   y=0, yend = -1 * score))
        if (score_flag == TRUE)
            pp <- pp + geom_hline(y=-0.5, linetype="longdash", color="steelblue")
        pp <- pp + scale_y_continuous(limits=c(-1, 1),
                                      breaks=seq(-1, 1, by=.25),
                                      labels=c(1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1))
        pp <- pp+geom_hline(y=0, color="black")
    }

    
    ## pp <- pp + scale_x_continuous(limits=c(0, 550), breaks=seq(0, 550, 50))
    ## pp <- pp + annotate("text", x=260, y=0.85, label="Glyco_H1 (983taxa, 6697 sites)")
    ## pp <- pp + annotate("text", x=260, y=-0.85, label="Glyco_H3 (1585 taxa, 13289 sites)")

    pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(title)
    pp <- pp + theme(legend.position="none")
    return(pp)
}

             



