#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# devtools::source_gist('7f63547158ecdbacf31b54a58af0d1cc', filename = 'util.R')
# load needed packages. Note that for the shiny app, these need to be preinstalled with root priviliges (sudo su - -c "R -e \"install.packages('pacman', repos='http://cran.rstudio.com/')\"")
# install.deps("pacman") 

CRAN_packages <- c("shiny", "dplyr", "purrr", "sqldf", "dbplyr", "DT", "RSQLite", "shinydashboard", 
                   "shinyWidgets", "Biostrings", "R.utils", "ggplot2", "plotly", "Cairo") # tidyverse
pacman::p_load(char=CRAN_packages)


# Find all available files
db_files <- list.files("/mnt/Shinotate_data/Trinotate_dbs", pattern = "[tT]rinotate.*\\.sqlite", 
                       full.names = TRUE)
dbs_list <- setNames(as.list(db_files), basename(db_files))
##### Specify default parameters #####
def_query <- "SELECT * FROM blast_annotation" # ORF_annotation / BlastDbase / blast_annotation
def_db <- "P_maxima_mantle_Trinotate_local.sqlite"
def_con <- DBI::dbConnect(RSQLite::SQLite(), dbs_list[[def_db]]) # /srv/shiny-server/bioinfo03/data/Paspaley/Annotation_dbs
num_results <- 1
# annot_query <- def_query
get_annotation <- function(db_src, query){
  dplyr::tbl(db_src, dplyr::sql(query)) %>% 
    dplyr::collect(n=Inf)
}

# Calculate ExN50
calc_N50 <- function(tx_table){
  N50 <- tx_table %>% mutate(trans_length=nchar(sequence)) %>% arrange(desc(trans_length)) %>% 
    # collect(n=Inf) %>%
    mutate(cum_sum=cumsum(trans_length)) %>% 
    mutate(N50=ifelse(cum_sum <= .$cum_sum[nrow(.)]/2, TRUE, FALSE)) %>% filter(N50==TRUE) %>%
    slice(n())
  N50$trans_length
} 
calc_ExN50 <- function(ex, exp_data, tx_table){
  tx_ids <- exp_data %>% filter(percentile<=ex) %>% .$feature_name
  tx_subset <- tx_table %>% filter(TrinityID %in% tx_ids) #%>% collect(n=Inf)
  dplyr::tibble(Ex=ex, ExN50=calc_N50(tx_subset), 
         num_transcripts=length(tx_ids))
} 

# ExN50 <- dplyr::tbl(def_con, "ExN50_stats") %>% dplyr::rename(E.N50=`E-N50`) %>%
#   collect() %>% mutate_at(-1, as.numeric) %>% mutate(Ex=as.numeric(gsub("E", "", .[["#E"]])))



# Define UI for application that draws a histogram
ui <- dashboardPage(skin = "purple",
  dashboardHeader(title = "Shinotate - Transcriptome annotation app", titleWidth = 450),

  dashboardSidebar(sidebarMenu(id = "sidebar",
                     menuItem("Sample Information", tabName = "samples",
                                  icon = icon("leaf", lib="glyphicon")),#, badgeLabel = "advanced", badgeColor = "green"),
                     menuItem("Transcript Annotation", tabName = "annotation",
                              icon = icon("list-alt", lib="glyphicon")), # icons: icon("dna", lib="font-awesome)
                     menuItem("Transcript Expression", tabName = "expression", 
                              icon = icon("bar-chart-o")),
                     menuItem("Transcriptome Stats", tabName = "stats",
                              icon = icon("info-circle", lib="font-awesome")),
                     menuItem("Experimental Metadata", tabName = "metadata",
                              icon = icon("book", lib="glyphicon")), #, badgeLabel = "advanced", badgeColor = "green"),
                     br(),
                     imageOutput("image"),
                     dropdownButton(
                       uiOutput('resetable_input'), size="default",
                       circle = TRUE, status = "success", icon = icon("cog", lib="font-awesome"),  # database 
                       width = "600px",
                       tooltip = tooltipOptions(title = "Trinotate db settings")  ),
                     
                     menuItem("Shinotate code and wiki", icon = icon("file-code-o"), 
                              href = "https://github.com/IdoBar/shiny-server/wiki")#, 
                              #badgeLabel = "info", badgeColor = "green")
                    )
                   ),
  dashboardBody(
    tags$head(tags$style(
      type="text/css",
      "#image img {max-width: 100%; width: 90%; height: auto}"
    )),
    tabItems(
      # Sample info tab content
      tabItem(tabName = "samples",
              fluidRow(
                box(
                  h2("Sample information"),
                  p("Mantle tissues from pearl oysters were collected from two pearl families, differing in the average quality of their produced pearls. The mantle of each individual was dissected to 3 distinct regions: proximal, distal and central to the foot base."), width=10),
                box(title = "Samples table",
                    div(DT::dataTableOutput("samples_table")),
                    br(),
                    div(
                      actionButton("excl_selected", label = "Exclude selected samples")# ,
                      # actionButton("select_all", label = "Select all")# ,
                    ) , width=5),
                box(, width=5))),
      # Annotation tab content
      tabItem(tabName = "annotation",
              fluidRow(
                box(
                  h2("Pearl Oyster (", em("Pinctada maxima"), ") mantle transcriptome database"), 
                  p("Mantle tissues from pearl oysters were collected, RNA was extracted and was subjected to RNA-Sequencing."),
                  p("The resulting data was ", em("de novo"), 
                    " assembled into a reference transcriptome, annotated and stored in a ",
                    a("Trinotate", href="http://trinotate.github.io/") ," database"),
                  p("Use the search bar on the right of the table to retrieve transcripts matching the keyword (at any field). The table can then be copied/printed/exported using the buttons on the left."),
                  p("Finally, selected transcripts can be exported in FASTA format (use the ", 
                    strong("Selection and Download"), " buttons below the table)."), width=11
                ),
                br(), br(),
                box(title = "Transcript annotation table",
                  div(DT::dataTableOutput("output_table")),
                  br(),
                  div(
                    downloadButton("download_fasta", 
                                   label = "Download sequences of selected transcripts (FASTA)"),
                    actionButton("clear_selection", label = "Clear selection"),
                    actionButton("select_all", label = "Select all")# ,
                  ) , width=11)
              )
      ),
      # Transcriptome stats tab content
      tabItem(tabName = "stats",
              fluidRow(
                box(
                  h2("Transcriptome summary statistics"),
                  h3("ExN50 plot:"),
                  p("An alternative to the Contig Nx statistic that could be considered more appropriate for transcriptome assembly data is the ExN50 statistic. Here, the N50 statistic is computed as above but limited to the top most highly expressed transcripts that represent x% of the total normalized expression data (Ex).",
                    "Additional information on the ExN50 statistic can be found at the ", a("Trinity documentation wiki", href="https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats"), "."),
                  actionButton("plot_ExN50", label = "Calculate and Plot ExN50"),
                  br(),
                  uiOutput('ExN50'), width=11 ))
              ),
    # Expression tab content
    tabItem(tabName = "expression",
            h2("Transcript expression")),
    # Metadata (analysis pipeline) tab content
    tabItem(tabName = "metadata",
            fluidRow(
              box(
                h2("Metadata (experimental design and analysis pipeline)"),
                h3("Pearl Oyster (", em("Pinctada maxima"), ") mantle transcriptome database"), 
                p("Mantle tissues from pearl oysters were collected, RNA was extracted and was subjected to RNA-Sequencing to identify differentially expressed genes between different regions of the mantle muscle, and high/low quality pearl producing genetic families."),
                p("The resulting data was ", em("de novo"), 
                  " assembled into a reference transcriptome, which was annotated and stored in a ",
                a("Trinotate", href="http://trinotate.github.io/") ," database. Transcript and gene counts were estimated (using RSEM) and added to the database." ), width=11)))# ,
    
    # Trinotate setup tab content (currently settings are in a dropdownButton in the sidebar)
    # tabItem(tabName = "settings") #,
            # uiOutput('resetable_input')
            
            # h2("Trinotate db setup"),
            # box(div(selectInput("db", "Trinotate db file", choices = dbs_list,
            #                     selected = dbs_list$P_maxima_mantle_Trinotate_local.sqlite)),
            #     h4("Trinotate db information:"),
            #     div(textOutput("db_info")), width = 10),
            #   box(div( textInput("query", "Annotation query", def_query),
            #          actionButton("update_query", "Update table"))),
            # tags$hr(),
            # br(),
            # box(div(actionButton("restore_defaults", "Restore defaults")), width = 2)
            # )
    )
  )
)
    

  
      
   
# )

# Define server logic required to show a results table of transcripts matching the search
server <- function(input, output, session) {
  # Render image
  output$image <- renderImage({
    list(
      src = "/srv/shiny-server/shinotate/www/1200px-Pinctada_margaritifera_MHNT.CON.2002.893.jpg",
      # # contentType = "image/jpg", 
      height = 120,
      alt = "Pinctada maxima"
    )}, deleteFile = FALSE)
  
  
  # resettable dropdownButton
  output$resetable_input <- renderUI({
    times <- input$restore_defaults
    div(id=letters[(times %% length(letters)) + 1], style = "color: black !important;",
        h2("Trinotate db setup"),    
        
    box(
        div(selectInput("db", "Trinotate db file", choices = dbs_list,
                        selected = dbs_list[[def_db]])),
        h4("Trinotate db information:"),
        div(verbatimTextOutput("db_info")), width = 11), # , title = "Trinotate db setup"
        
        # textInput("query", "Annotation query", def_query),
        # actionButton("update_query", "Update table", width = '150px'),
        # actionButton("restore_defaults", "Restore defaults",  width = '150px'))
        
    box(textInput("query", "Annotation query", def_query),
        actionButton("update_query", "Update table", width = '150px', 
                     style = "display: inline-block !important;"),
             actionButton("restore_defaults", "Restore defaults",  width = '150px', 
                          style = "display: inline-block !important;"), 
        width = 11)
    )
  }) 
  
  
 #### Pull annotation table from Trinotate ####
  # make a reactive connection to the db 
  C1  <- reactive({
    dbplyr::src_dbi(DBI::dbConnect(RSQLite::SQLite(), input$db))
  })
  # print db information
  output$db_info <- renderText({
    db_info <- R.utils::captureOutput(print(C1()))
    sprintf("%s\n%s", db_info[1], paste(db_info[-1], collapse = "\n"))
  }
    )
  # initialise the table with default settings
  values <- reactiveValues(annot_table = get_annotation(dbplyr::src_dbi(def_con), def_query),
                           tx_table = dplyr::tbl(def_con, "Transcript"),
                           orf_table = dplyr::tbl(def_con, "ORF"),
                           exp_table = dplyr::tbl(def_con, "Expression"),
                           rep_table = dplyr::tbl(def_con, "Replicates"),
                           sample_table = dplyr::tbl(def_con, "samples_table"),
                           show_exn50 = FALSE
   )  
  # update the table after db and/or query change
  observeEvent(c(input$update_query), { # , input$restore_defaults
    values$annot_table <- get_annotation(C1(), input$query)
    values$tx_table <- dplyr::tbl(C1(), "Transcript")
    values$orf_table <- dplyr::tbl(C1(), "ORF")
    values$exp_table <- dplyr::tbl(C1(), "Expression")
    values$rep_table = dplyr::tbl(def_con, "Replicates")
    values$sample_table = dplyr::tbl(def_con, "samples_table")
    
  })
  observeEvent(input$plot_ExN50,{
    values$show_exn50 <- TRUE
  })
  
  
  #### Samples table ####
  output$samples_table <- DT::renderDataTable(values$sample_table %>% collect(n=Inf),  
                                              rownames= F, # filter = "top",
                             extensions = c("Buttons", "Scroller",'FixedColumns',"FixedHeader"),
                             options = list(
                               dom = 'Bfrtip',
                               responsive = TRUE,
                               # autoWidth = TRUE,
                               # columnDefs = list(list(width = '100px', targets = c(1, 5))),
                               pageLength = 12,
                               buttons =list('copy','print',
                                             list( extend = 'collection',
                                                   buttons = c('csv', 'excel'), # , 'pdf'
                                                   text = 'Download')) #,
                                             # list(extend = 'colvis', columns = 1:5)),
                               # scrollY = 350,# deferRender = TRUE,
                               # scroller = TRUE,
                               # scrollX = FALSE,
                               # "search" = list("search" = "Mantle"),
                               # # paging=FALSE,
                               # fixedHeader=TRUE#,
                               # fixedColumns = list(leftColumns = 1, rightColumns = 0)#,  lengthMenu=20
                            ))
  # # Row selection in table 
  # set the table as proxy to be able to update it
  samples_proxy = dataTableProxy('samples_table')
  # Identify button click 
  observeEvent(input$excl_selected, {
    samples_proxy %>% selectRows(list())
  })
  
  # observeEvent(input$select_all, {
  #   # print(as.numeric(input$output_table_rows_all))
  #   proxy %>% selectRows(as.numeric(input$output_table_rows_all))
  # })
  
  output$ExN50 <- renderUI({
    if (values$show_exn50==TRUE){
      withProgress(message = 'Calculating ExN50', value = 0,{
        incProgress(1/4, detail = "Retrieving transcripts")
        # retrieve transcripts
        tx_data <- values$tx_table %>% dplyr::select(TrinityID=transcript_id, sequence) %>%
          dplyr::collect(n=Inf)
        incProgress(2/4, detail = "Retrieving expression data")
        # Retrieve Expression data
        exp_data <-  values$exp_table %>% dplyr::filter(feature_type=="T") %>% 
          dplyr::group_by(feature_name) %>% dplyr::summarise(feature_count=sum(fpkm)) %>% 
          dplyr::filter(feature_count>0) %>% dplyr::collect(n=Inf) %>% 
          dplyr::arrange(desc(feature_count)) %>% 
          dplyr::mutate(cum_sum=cumsum(feature_count)) %>% 
          dplyr::mutate(percentile=cum_sum/sum(feature_count)*100)
        # Calculate ExN50 (see https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)
        incProgress(3/4, detail = "Calculating ExN50")
        ExN50_table <<- purrr::map_df(1:100, ~ dplyr::bind_rows(calc_ExN50(., exp_data, 
                                                                           tx_data)))
        E100N50 <<- ExN50_table %>% dplyr::filter(Ex==100)
        EmaxN50 <<- ExN50_table %>% dplyr::filter(ExN50==max(ExN50))
      })
      div(plotlyOutput('ExN50Plot', height = "350px"),
          br(),
      p("The ExN50 per each Ex plot shows that the effective N50 'improved' from ", 
            strong(sprintf(" %s bp to %s bp (at E%sN50)", 
                           prettyNum(E100N50$ExN50, big.mark = ","), 
                           prettyNum(EmaxN50$ExN50, big.mark = ","),EmaxN50$Ex )),
        "compared with the N50 caluclated with all transcripts (E100)."))
    } else {
      p("Click on the button above to calculate and show ExN50 plot")
    }
    
  })
  
  # Render plotly 
  output$ExN50Plot <- renderPlotly({
    
    g <- ggplot(ExN50_table, aes(x=Ex, y=ExN50)) +
      geom_line(col = "blue", size=1) + ylab("ExN50 (bp)")+ xlab("Ex (%)") +
      theme_bw(base_size = 16) +
      scale_y_continuous(expand = c(0,0), limits = c(0,ceiling(EmaxN50$ExN50/1000)*1000)) + 
      scale_x_continuous(breaks=c(0,25,50,75,EmaxN50$Ex,100)) +
      geom_segment(data=EmaxN50, xend=EmaxN50$Ex, yend=0,colour="pink", linetype=2, lwd=0.8) +
      geom_point(data = EmaxN50, size=2.5, shape=21, fill="pink", colour="darkred") +
      # annotate("text", x=EmaxN50$Ex-2,y=EmaxN50$E.N50+50, label=sprintf("E%sN50 = %sbp",EmaxN50$Ex,
      #                                                  prettyNum(EmaxN50$E.N50, big.mark = ",")), hjust=1) +
      geom_segment(data=E100N50, xend=E100N50$Ex, yend=0,colour="pink", linetype=2, lwd=0.8) +
      geom_point(data = E100N50, size=2.5, shape=21, fill="pink", colour="darkred")# +
    # annotate("text", x=E100N50$Ex-2,y=E100N50$E.N50+50, label=paste("N50 =",
    #                                                 prettyNum(E100N50$E.N50, big.mark = ","), "bp"), hjust=1)
    NULL    
    ggplotly(g)
  })
  
  #### Annotation table ####
  output$output_table <- DT::renderDataTable(values$annot_table,  rownames= F, # filter = "top",
                               extensions = c("Buttons", "Scroller",'FixedColumns',"FixedHeader"),
                               options = list(
                                 dom = 'Bfrtip',
                                 # autoWidth = TRUE,
                                 columnDefs = list(list(width = '150px', targets = c(2, 5, 7))),
                                 # pageLength = 15,
                                 buttons =list('copy','print', 
                                               list( extend = 'collection',
                                       buttons = c('csv', 'excel'), # , 'pdf'
                                       text = 'Download'),
                                       list(extend = 'colvis', columns = 1:8)),
                                 scrollY = 350,# deferRender = TRUE,
                                 scroller = TRUE,
                                 scrollX = TRUE,
                                 # paging=FALSE,
                                 fixedHeader=TRUE#,
                                 # fixedColumns = list(leftColumns = 1, rightColumns = 0)#,  lengthMenu=20
                                 
                                 ))
  # Row selection in table 
  # set the table as proxy to be able to update it
  proxy = dataTableProxy('output_table')
  # Identify button click 
  observeEvent(input$clear_selection, {
    proxy %>% selectRows(list())
  })
  
  observeEvent(input$select_all, {
    # print(as.numeric(input$output_table_rows_all))
    proxy %>% selectRows(as.numeric(input$output_table_rows_all))
  })
  
  #### Download sequences as fasta  ####
  # curr
  output$download_fasta <- downloadHandler(
    filename = "selected_sequences.fasta",
    content = function(file) {
      # input_output_table_rows_selected = 1:10
      ids <- values$annot_table$TrinityID[input$output_table_rows_selected]
      # ids <- mantle_annot$TrinityId[input_output_table_rows_selected]
      tx_ids <- unique(sub("\\|m\\..+", "", ids))
      cds_ids <- ids[grepl("\\|m\\..+", ids)]
      # extract transcript sequences
      tx <- values$tx_table %>% dplyr::filter(transcript_id %in% tx_ids) %>% 
        dplyr::select(transcript_id, sequence) %>% dplyr::collect(n=Inf)
      # extract just the coding sequences (reverse complement if on the minus strand)
      cds <- values$orf_table %>% dplyr::filter(orf_id %in% cds_ids) %>%  
        dplyr::collect(n=Inf) %>% dplyr::mutate(cds=substr(tx$sequence[match(transcript_id, tx$transcript_id)], 
                                 lend, rend), cds_2=ifelse(strand=="+", cds, chartr("ATGC","TACG",reverse(cds)))) %>% 
        dplyr::select(transcript_id=orf_id, sequence=cds_2) 
      # combine both whole transcripts and cds
      seqs_out <- tx %>% dplyr::filter(transcript_id  %in% ids[!grepl("\\|m\\..+", ids)]) %>% 
        dplyr::bind_rows(cds)
      # save as fasta file (nucleotide)
      seqinr::write.fasta(as.list(seqs_out$sequence), names = seqs_out$transcript_id, 
                          file.out = file)
      
      
    }
  )
  
}


# Run the application 
shinyApp(ui = ui, server = server)

