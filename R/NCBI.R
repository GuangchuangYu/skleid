##' download genbank (.gb) file by GI number
##'
##' 
##' @title download_genbank_gi 
##' @param gi gi number
##' @param db supported db, currently 'nuccore'
##' @param outfile output file, by default, gi_number.gb
##' @return NULL
##' @export
##' @author Guangchuang Yu
download_genbank_gi <- function(gi, db="nuccore", outfile=NULL) {
    if (db != "nuccore") {
        stop("currently, only nuccore is supported...")
    }
    url <- paste0("www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=",
                  db, '&dopt=genbank&sort=&val=', gi)

    if (is.null(outfile)) {
        outfile <- paste0(gi, ".gb")
    }
    
    download.file(url=url, destfile = outfile, method="curl")
}

##' download genbank (.gb) file by accession number
##'
##' 
##' @title download_genbank_acc
##' @param acc accession number
##' @param db supported db, currently 'nuccore'
##' @param outfile output file, by default, acc_number.gb
##' @return NULL
##' @export
##' @author Guangchuang Yu
download_genbank_acc <- function(acc, db="nuccore", outfile=NULL) {
    if (is.null(outfile)) {
        outfile <- paste0(acc, ".gb")
    }
    gi <- acc2gi(acc)
    download_genbank_gi(gi, db, outfile)
}

acc2gi <- function(acc) {
    url <- paste0("http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=", acc)
    x <- readLines(url)
    i <- grep('<Id>', x)
    if (length(i) == 0) {
        stop("GI number not found...")
    }
    id <- x[i]
    gi <- gsub('</*Id>', '', id)
    return(gi)
}

##' parse genbank file and extract infomation
##' (not a universal parse, just extract infomation of my own needs)
##'
##' 
##' @title read.genbank
##' @param gbfile genbank file
##' @return data.frame
##' @export
##' @author Guangchuang Yu
read.genbank <- function(gbfile) {
    gb <- readLines(gbfile)

    source <- get_gb_range(gb, "source")
    gene <- get_gb_range(gb, "CDS")

    gn <- get_gb_name(gb, "gene")
    org <- get_gb_organism(gb)
    
    year <- get_year_from_organism_name(org)

    seq <- get_gb_seq(gb)
    
    data.frame(acc       = get_gb_acc(gb),
               gi        = get_gb_gi(gb),
               name      =org,
               year      =year,
               start     =source[1],
               end       =source[2],
               geneName  = gn,
               geneStart =gene[1],
               geneEnd   =gene[2],
               seq       = seq)
}

get_gb_seq <- function(gb) {
    i <- grep("ORIGIN", gb)
    i <- i+1
    j <- grep("//", gb)
    j <- j-1
    seq <- gb[i:j]
    seq <- gsub("\\d+", "", seq) %>% gsub("\\s+", "", .)
    seq <- paste0(seq, collapse="")
    return(seq)
}

get_gb_range <- function(gb, field) {
    gb[grep(paste0(field, "\\s"), gb)] %>%
        strsplit(., "\\s+") %>% unlist %>% `[`(3) %>%
            strsplit(., "\\.\\.>*") %>% unlist %>% as.numeric
}

get_gb_name <- function(gb, field) {
    gb[grep(paste0(field, "="), gb)][1] %>%
        strsplit(., "=") %>% unlist %>% `[`(2) %>%
            gsub('\"', "", .)
}

get_year_from_organism_name <- function(org) {
    sub("\\(H\\d+N\\d+\\)", "", org) %>%
        sub("\\)", "", .) %>%
            gsub(".*/(\\d+)$", "\\1", .) %>%
                as.numeric
}

get_gb_organism <- function(gb) {
    i <- grep("SOURCE", gb)
    j <- grep("ORGANISM", gb)
    gb[i:(j-1)] %>% paste0(., collapse = " ") %>%
        sub("SOURCE\\s+", "", .) %>% gsub("\\s+", " ", .)
}

get_gb_acc <- function(gb) {
    gb[grep("ACCESSION", gb)] %>%
        strsplit("\\s+") %>% unlist %>%
            `[`(2)
}

get_gb_gi <- function(gb) {
    gb[grep("VERSION", gb)] %>%
        strsplit("\\s+") %>% unlist %>%
            `[`(3)
}
