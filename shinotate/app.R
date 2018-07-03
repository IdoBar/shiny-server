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
# Install and load bioconductor packages
# bioc_packages <- c("Biostrings")
# install.deps(bioc_packages, repo="bioc")
CRAN_packages <- c("shiny", "dplyr", "sqldf", "dbplyr", "DT", "RSQLite", "shinydashboard",                              "Biostrings", "R.utils") # tidyverse
pacman::p_load(char=CRAN_packages)


# Find all available files
db_files <- list.files("/mnt/Shinotate_data/Trinotate_dbs", pattern = "[tT]rinotate.*\\.sqlite", 
                       full.names = TRUE)
dbs_list <- setNames(as.list(db_files), basename(db_files))

def_query <- "select * from blast_annotation" # ORF_annotation / BlastDbase / blast_annotation
num_results <- 1
annot_query <- def_query

# Define UI for application that draws a histogram
ui <- dashboardPage(skin = "purple",
  dashboardHeader(title = "Shinotate - Transcriptome annotation app", titleWidth = 450),

  dashboardSidebar(sidebarMenu(id = "sidebar",
                     menuItem("Transcript Annotation", tabName = "annotation",
                              icon = icon("list-alt", lib="glyphicon")),
                     menuItem("Transcript Expression", tabName = "de", icon = icon("bar-chart-o")),
                     
                     br(),
                     imageOutput("image"),
                     menuItem("Trinotate db setup", tabName = "settings",
                              icon = icon("gear"), badgeLabel = "advanced", badgeColor = "green"),
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
      # First tab content
      tabItem(tabName = "annotation",
    fluidRow(
      box(
    h2("Pearl Oyster (", em("Pinctada maxima"), ") mantle transcriptome database"), 
                p("Mantle tissues from pearl oysters were collected, RNA was extracted and was subjected to RNA-Sequencing."),
                p("The resulting data was ", em("de novo"), " assembled into a reference transcriptome, annotated and stored in a ",
                  a("Trinotate", href="http://trinotate.github.io/") ," database"),
                p("Use the search bar on the right of the table to retrieve transcripts matching the keyword (at any field). The table can then be copied/printed/exported using the buttons on the left."),
            p("Finally, selected transcripts can be exported in FASTA format (use the ", strong("Selection and Download"), " buttons below the table)."), width=11),
                br(), br(),
    
                box(div(DT::dataTableOutput("output_table")),
                br(),
                div(downloadButton("download_fasta", 
                        label = "Download sequences of selected transcripts (FASTA)"),
                    # div(class="horizontalgap", style="width:5px"),# hr(),
                    actionButton("clear_selection", label = "Clear selection"),
                    actionButton("select_all", label = "Select all")) , width=11)
                # box(textOutput("selected_var"))
        )
      )
    ,
    
    # Second tab content
    tabItem(tabName = "de",
            h2("Differential expression tab content")),
    tabItem(tabName = "settings",
            uiOutput('resetable_input'),
            
            # h2("Trinotate db setup"),
            # box(div(selectInput("db", "Trinotate db file", choices = dbs_list,
            #                     selected = dbs_list$P_maxima_mantle_Trinotate_local.sqlite)),
            #     h4("Trinotate db information:"),
            #     div(textOutput("db_info")), width = 10),
            #   box(div( textInput("query", "Annotation query", def_query),
            #          actionButton("update_query", "Update table"))),
            tags$hr(),
            br(),
            box(div(actionButton("restore_defaults", "Restore defaults")), width = 3)
            )
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
  
  # res <- reactive({mantle_annot %>% dplyr::filter(grepl(input$keyword, Description, ignore.case = TRUE))})
  output$resetable_input <- renderUI({
    times <- input$restore_defaults
        div(id=letters[(times %% length(letters)) + 1],
            h2("Trinotate db setup"),
        box(div(selectInput("db", "Trinotate db file", choices = dbs_list,
                            selected = dbs_list$P_maxima_mantle_Trinotate_local.sqlite)),
            h4("Trinotate db information:"),
            div(verbatimTextOutput("db_info")), width = 8),
        box(div( textInput("query", "Annotation query", def_query),
                 actionButton("update_query", "Update table")), width = 8))
  })
  
  con <- DBI::dbConnect(RSQLite::SQLite(), "/mnt/Shinotate_data/Trinotate_dbs/P_maxima_mantle_Trinotate_local.sqlite") # /srv/shiny-server/bioinfo03/data/Paspaley/Annotation_dbs
  
  # annot_query <- eventReactive(input$update_query, {
  #   input$query
  # })
  # observeEvent(D1(),{
  #   updateSelectInput(session, "selectinputid", "Language to Select:",  choices = unique(D1()$Language),selected = unique(D1()$Language)[1])
  # })
  
  C1  <- reactive({
    dbplyr::src_dbi(DBI::dbConnect(RSQLite::SQLite(), input$db))
  })
  output$db_info <- renderText({
    db_info <- R.utils::captureOutput(print(C1()))
    sprintf("%s\n%s", db_info[1], paste(db_info[-1], collapse = "\n"))
  }
    )
  # con_values <- reactiveValues(con = DBI::dbConnect(RSQLite::SQLite(), input$db))
  values <- reactiveValues(annot_table = {
    dplyr::tbl(dbplyr::src_dbi(con), dplyr::sql(annot_query)) %>% 
      dplyr::collect(n=Inf) # %>% dplyr::group_by(1) %>%
      # dplyr::arrange(dplyr::desc(BitScore)) %>%
      # dplyr::slice(1:num_results)
  } )  # %>% dplyr::group_by(TrinityId) %>%
    # dplyr::arrange(dplyr::desc(BitScore)) %>%
    # dplyr::slice(1:num_results)
  # dplyr::filter(between(row_number(), 1, num_results)) %>%
  observeEvent(input$update_query, {
    values$annot_table <- dplyr::tbl(dbplyr::src_dbi(con), dplyr::sql(input$query)) %>% 
      dplyr::collect(n=Inf)
    
  })
  # return settings to default
  D2 <- eventReactive(input$goButton1,{
    D1()[D1()$Language %in% input$selectinputid,]
  })
  
  output$output_table <- DT::renderDataTable(values$annot_table,  rownames= F, # filter = "top",
                               extensions = c("Buttons", "Scroller",'FixedColumns',"FixedHeader"),
                               options = list(
                                 dom = 'Bfrtip',
                                 # autoWidth = TRUE,
                                 columnDefs = list(list(width = '150px', targets = c(2, 5, 7))),
                                 # pageLength = 15,
                                 buttons =list('copy','print', 
                                               list( extend = 'collection',
                                       buttons = c('csv', 'excel', 'pdf'),
                                       text = 'Download'),
                                       list(extend = 'colvis', columns = 1:8)),
                                 scrollY = 350,# deferRender = TRUE,
                                 scroller = TRUE,
                                 scrollX = TRUE,
                                 # paging=FALSE,
                                 fixedHeader=TRUE#,
                                 # fixedColumns = list(leftColumns = 1, rightColumns = 0)#,  lengthMenu=20
                                 
                                 ))
  # set the table as proxy to be able to update it
  proxy = dataTableProxy('output_table')
  # Identify button click 
  observeEvent(input$clear_selection, {
    proxy %>% selectRows(list())
  })
  
  observeEvent(input$select_all, {
    print(as.numeric(input$output_table_rows_all))
    proxy %>% selectRows(as.numeric(input$output_table_rows_all))
  })
  
  # output$selected_var <- renderText({ 
  #   paste("Your filtered table contains the following rows: ", 
  #         paste(as.numeric(input$output_table_rows_all), collapse = ","))
  # })
  
  output$download_fasta <- downloadHandler(
    filename = "selected_sequences.fasta",
    content = function(file) {
      # input_output_table_rows_selected = 1:10
      ids <- annot_table$TrinityID[input$output_table_rows_selected]
      # ids <- mantle_annot$TrinityId[input_output_table_rows_selected]
      tx_ids <- unique(sub("\\|m\\..+", "", ids))
      cds_ids <- ids[grepl("\\|m\\..+", ids)]
      # extract transcript sequences
      tx <- dplyr::tbl(src_dbi(con), "Transcript") %>% dplyr::filter(transcript_id %in% tx_ids) %>% 
        dplyr::select(transcript_id, sequence) %>% dplyr::collect(n=Inf)
      # extract just the coding sequences
      cds <- dplyr::tbl(src_dbi(con), "ORF") %>% dplyr::filter(orf_id %in% cds_ids) %>%  
        dplyr::collect(n=Inf) %>% dplyr::mutate(cds=substr(tx$sequence[match(transcript_id, tx$transcript_id)], lend, rend), cds_2=ifelse(strand=="+", cds, chartr("ATGC","TACG",reverse(cds)))) %>% 
        dplyr::select(transcript_id=orf_id, sequence=cds_2) 
      # combine both whole transcripts and cds
      seqs_out <- tx %>% dplyr::filter(transcript_id  %in% ids[!grepl("\\|m\\..+", ids)]) %>% 
        dplyr::bind_rows(cds)
      # file="test.fasta"
      # tempFile <- file.path(tempdir(), "selected_sequences.fasta")
      
      seqinr::write.fasta(as.list(seqs_out$sequence), names = seqs_out$transcript_id, 
                          file.out = file)
      # file.copy(tempFile, file)
      # file.copy(filename(),file )
      
    }
  )
  
}

# the indices of the selected rows are available through input$tableId_rows_selected
# Run the application 
shinyApp(ui = ui, server = server)

