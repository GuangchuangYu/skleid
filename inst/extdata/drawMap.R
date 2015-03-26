
##' get map coordination
##'
##' 
##' @title get.map
##' @param map only 'china' is supported
##' @return data.frame
##' @import maps
##' @import mapdata
##' @export
##' @author Guangchuang Yu
get.map <- function(map="china") {
    if (map != "china") {
        stop("currently, only 'china' is supported...")
    }
    mapdata <- "mapdata"
    require(mapdata, character.only=TRUE)
    coord <- map.poly(map)
    df <- data.frame(lon=coord$x, lat=coord$y)
    return(df)
}

map.poly <- maps:::map.poly

##' plot map
##'
##' 
##' @title plotMap
##' @param mapdata map coordination 
##' @param xlim xlim
##' @param ylim ylim
##' @param location location matrix
##' @param lineWidth linewidth matrix
##' @param maxLineSize maxLineSize
##' @param color color matrix
##' @param colorIntensity color intensity matrix 
##' @param arrowHeadSize arrow head size
##' @param matrixFlag one of upper or lower
##' @param mapcol color of background map
##' @param ggmap whether use ggmap or not
##' @param maptype map type of ggmap
##' @return ggplot2 figure
##' @importFrom reshape2 melt
##' @importFrom geosphere gcIntermediate
##' @importFrom ggmap ggmap
##' @importFrom ggmap get_map
##' @importFrom grid arrow
##' @importFrom grid unit
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_path
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 ylim
##' @importFrom ggplot2 theme_classic
##' @import DOSE
##' @export
##' @author Guangchuang Yu
plotMap <- function(mapdata, xlim=c(110, 125), ylim=c(20,35),
                    location, lineWidth, maxLineSize = 3,
                    color, colorIntensity, arrowHeadSize = 1,
                    matrixFlag="upper", mapcol="darkgrey",
                    ggmap=FALSE, maptype="hybrid") {

    lon <- lat <- NULL
    if (ggmap == TRUE) {
        ## maptype: "terrain", "satellite", "roadmap", "hybrid"
        mm <- get_map(location=c(mean(xlim), mean(ylim)), zoom=5, maptype=maptype)
        p <- ggmap(mm)
        p <- p+xlim(xlim)+ylim(ylim)
    } else {
        p <- ggplot(mapdata, aes(lon, lat))
        p <- p+geom_path(size=0.5, color=mapcol)+xlim(xlim)+ylim(ylim)
    }
    
    ## loc2 <- location
    ## loc2$m <- rowMeans(mat, na.rm=T)
    ## p <- p+geom_point(data=loc2, aes(x=V3, y=V2, size=m), alpha=1/4, colour="#b5e521")+scale_size_continuous(range=c(10, 20), name=legendTitle)

    lwd.max <- max(lineWidth, na.rm=TRUE)
    colInt.max <- max(colorIntensity, na.rm=TRUE)
    if (matrixFlag == "upper") {
        for (i in 1:nrow(lineWidth)) {
            for (j in 1:i) {
                lineWidth[i,j] <- NA
            }
        }
    } else if (matrixFlag == "lower") {
        for (i in 1:nrow(lineWidth)) {
            for (j in i:ncol(lineWidth)) {
                lineWidth[i,j] <- NA
            }
        }
        
    }

    lineWidth$id <- rownames(lineWidth)
    lwd.df <- melt(lineWidth, id.vars="id")
    lwd.df <- lwd.df[!is.na(lwd.df$value),]
    lwd.df <- lwd.df[order(lwd.df$value),]
    lwd.df[,2] <- as.character(lwd.df[,2])
    lwd.df$value <- lwd.df$value/lwd.max * maxLineSize

    colorNN <- apply(color[,1:2], 1, paste0, collapse="")
    lwdNN <- apply(lwd.df[,1:2], 1, paste0, collapse="")
    ii <- match(colorNN, lwdNN)
    ii <- ii[!is.na(ii)]
    
    for (idx in ii) {
        lwd = lwd.df$value[idx]

        col <- as.character(color$V3[color$V1==lwd.df[idx,1] & color$V2==lwd.df[idx,2]])
        cols <- getCol(col)
        ci <- DOSE:::getIdx(colorIntensity[lwd.df[idx,1], lwd.df[idx,2]], 0, colInt.max)
        col2 = cols[ci]
                
        p1 <- rev(as.vector(location[lwd.df[idx,1],]))
        p2 <- rev(as.vector(location[lwd.df[idx,2],]))
        inter <- gcIntermediate(p1=p1, p2=p2, addStartEnd=TRUE)
        
        p <- p+geom_path(data=as.data.frame(inter), aes(lon, lat),
                         size=lwd, color=col2,
                         arrow=arrow(length=unit(lwd/lwd.max * arrowHeadSize, "cm")) ) 
                
    }

    p <- p+theme_classic()
    #p <- p+theme(legend.position="none")
    return(p)
}

getCol <- function(col) {
    pal <- colorRampPalette(c("white", col))
    colors <- pal(100)
    return(colors)
}



##' plot map using input files
##'
##' 
##' @title plotMap2
##' @param mapdata map coordination 
##' @param xlim xlim
##' @param ylim ylim
##' @param location_file location file
##' @param lineWidth_file line width file
##' @param maxLineSize maxLineSize
##' @param color_file color file
##' @param colorIntensity_file color intensity file
##' @param arrowHeadSize arrow head size
##' @param matrixFlag one of upper or lower
##' @param mapcol color of background map
##' @param ggmap whether use ggmap or not
##' @param maptype map type of ggmap
##' @return ggplot2 figure
##' @export
##' @author Guangchuang Yu
plotMap2 <- function(mapdata, xlim=c(110, 125), ylim=c(20,35),
                    location_file, lineWidth_file, maxLineSize = 3,
                    color_file, colorIntensity_file, arrowHeadSize = 1,
                    matrixFlag="upper", mapcol="darkgrey",
                    ggmap=FALSE, maptype="hybrid") {

    
    color <- read.delim(color_file, header=FALSE, stringsAsFactor=FALSE)
    location <- read.delim(location_file, header=FALSE, row.names=1, stringsAsFactor=FALSE)
    lineWidth <- read.delim(lineWidth_file, header=TRUE, row.names=1, stringsAsFactor=FALSE)
    colorIntensity <- read.delim(colorIntensity_file, header=TRUE, row.names=1, stringsAsFactor=FALSE)

    plotMap(mapdata, xlim, ylim, location, lineWidth, maxLineSize,
            color, colorIntensity, arrowHeadSize, matrixFlag,
            mapcol, ggmap, maptype)
}
