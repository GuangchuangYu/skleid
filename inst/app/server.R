library("skleid")
library("shiny")
library("shinyFiles")

shinyServer(function(input, output, session) {
    roots <- c(root=getVolumes()())
    shinyDirChoose(input, 'wd_folder', session=session, roots=roots, filetypes=c('', '.txt'))

    output$set_wd <- renderUI({
        if (!is.null(input$wd_folder)) {
            textInput("wd", "", parseDirPath(roots=roots, input$wd_folder))
        }
    })

    get_files <- function() {
        if (!is.null(input$wd)) {
            files <- list.files(input$wd)
        } else {
            files <- ""
        }
        return(files)
    }
    
    output$choose_contig_folder <- renderUI({
        files <- get_files()
        if (input$action == "autoReport") {
            selectInput("contig_folder", "Contig Folder", as.list(files))
        }        
    })

    output$choose_output_folder <- renderUI({
        files <- get_files()
        if (input$action == "summarize") { 
            selectInput("out_folder", "Output Folder", as.list(files))
        } else if (input$action == "generateGapFile") {
            selectInput("out_folder", "Output Folder/Summary File", as.list(files))
        }
    })
    
    output$choose_ref_folder <- renderUI({
        files <- get_files()
        if (input$action != "summarize") {
            selectInput("ref_folder", "Ref Folder", as.list(files))
        } 
    })

    output$set_read_fileName <- renderUI({
        if (input$action == "generateGapFile") {
            files <- get_files()
            selectInput("read_fileName", "Filenames of Reads", as.list(files))
        }
    })

    output$set_sample_size <- renderUI({
        if (input$action == "generateGapFile") {
            textInput("sample_size", "Sample Size", value=200)
        }
    })
    
    output$choose_name_file <- renderUI({
        if (input$action != "generateGapFile") {
            files <- get_files()
            selectInput("name_file", "Name File", as.list(files))
        }
    })
    
    ## output$choose_summary_name_file <- renderUI({
    ##     if (input$action == "summarize") {
    ##         files <- get_files()
    ##         selectInput("sname_file", "Name File", as.list(files))
    ##     }
    ## })

    output$set_filter_percentage <- renderUI({
        if (input$action == "autoReport") {
            files <- get_files()
            textInput("percentage", "Filter Percentage", value = 1)
        }
    })

    is_empty <- function(x) {
        if (is.null(x)) {
            return(TRUE)
        }
        if (x == "") {
            return(TRUE)
        }
        return(FALSE)
    }

    run_skleid <- eventReactive(input$goButton, {
        setwd(input$wd)

        if (input$action == "autoReport") {
            autoReport(input$contig_folder, input$ref_folder, input$name_file, "output", filter=TRUE, percentage=input$percentage)
        } else if (input$action == "summarize") {
            summarize(input$out_folder, input$name_file)
        } else if (input$action == "generateGapFile") {
            generateGapFile(input$out_folder, input$ref_folder, input$read_fileName, input$sample_size)
        }
                               
    })

    output$skleid_info <- renderPrint({
        cat(paste("  skleid package, version = ", as.character(packageVersion("skleid")),
        sep = ""), "\n")
        cat("  Author: Guangchuang Yu (gcyu@connect.hku.hk)")
        ## printInfo()
    })
    
    output$skleid_output <- renderPrint({
        run_skleid()
    })

    
    output$file_output <- renderUI({
        if (input$action == "summarize") {
            dataTableOutput("summary")
        }
        
        file_url <- paste0("file://", input$wd, "/", "report.html")
        tags$link(file_url)
        tags$html(
            tags$body(a(href=file_url))
            )

        tags$iframe(src=file_url, seamless=NA)
        ##HTML(paste0("<p><a href='", file_url, "'>report.html</a></p>"))
    })

    
    
    output$summary <- renderDataTable({
        if (input$action == "summarize") {
            setwd(input$wd)
            read.csv("summary.csv")
        }
    })

    output$gaps_lines <- renderUI({
        if (!is.null(input$wd) && input$action == "generateGapFile") {
            setwd(input$wd)
            gaps <- readLines('gaps.txt')
            sliderInput("nlines", "Number of lines to view:", min=0, max=length(gaps),value=20)
        }
    })

    
    output$gaps <- renderPrint({
        if (!is.null(input$wd) && !is.null(input$nlines) && input$action == "generateGapFile") {
            setwd(input$wd)
            gaps <- readLines('gaps.txt')
            for (i in 1:input$nlines) {
                cat("[", i, "]", " ", gaps[i], "\n")
            }
            ## head(gaps, 50)
        }
    })

    output$browser <- renderUI({
        if (!is.null(input$wd) && input$action == "autoReport") {
            actionButton("browser_button", "view report.html in separate browser")
        }
    })

    run_browser <- eventReactive(input$browser_button, {
        url <- paste0("file://", input$wd, "/", "report.html")
        browseURL(url)
        br()
    })

    output$view_browser <- renderUI({
        run_browser()
    })
    
    output$html <- renderUI({
        if (!is.null(input$wd) && input$action == "autoReport") {
            setwd(input$wd)
            includeHTML("report.html")           
        }
    })
    
})


