shinyUI(fluidPage(
    titlePanel("SKLEID"),
    
    sidebarLayout(
        sidebarPanel(
            radioButtons("action", "Action:",
                         list("autoReport" = "autoReport",
                              "summarize" = "summarize",
                              "generateGapFile" = "generateGapFile"
                              ),
                         ),
            ## selectInput("action", "Action:",
            ##             list("autoReport" = "autoReport",
            ##                  "summarize" = "summarize",
            ##                  "generateGapFile" = "generateGapFile"
            ##                  ),
            ##             ),   
            shinyDirButton("wd_folder", label="Working Directory", title = "Please select a folder"),
            uiOutput("set_wd"),
            uiOutput("choose_contig_folder"),
            uiOutput("choose_ref_folder"),
            uiOutput("choose_name_file"),
            uiOutput("set_filter_percentage"),
            uiOutput("choose_output_folder"),
            uiOutput("set_read_fileName"),
            uiOutput("set_sample_size"),
            br(),
            actionButton("goButton", "Go!")
            ),
        
        mainPanel(
            verbatimTextOutput("skleid_info"),
            tabsetPanel(
                tabPanel("Console output",
                         verbatimTextOutput("skleid_output")
                         ),
                tabPanel("Final output",
                         uiOutput("browser"),
                         uiOutput("view_browser"),
                         uiOutput("html"),
                         dataTableOutput('summary'),
                         uiOutput("gaps_lines"),
                         verbatimTextOutput("gaps")
                         )
                
                )
            )
        )
    ))

