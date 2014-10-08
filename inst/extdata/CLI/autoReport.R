#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-c", "--contig_folder"), help="Contig folder"),
    make_option(c("-r", "--ref_folder"), help="Reference folder"),
    make_option(c("-n", "--name_file"), help="consensus fasta name scheme file"),
    make_option(c("-o", "--out_folder"), help="output folder"),
    make_option(c("-f", "--filter"), help="TRUE or FALSE, perform read number filter or not"),
    make_option(c("-p", "--percentage"), help="below the specific percentage of read number  will be filter, only effective when --filter=TRUE")
    )


opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("skleid"))
autoReport(opt$contig_folder, opt$ref_folder, opt$name_file, opt$out_folder, opt$filter, opt$percentage)


    
