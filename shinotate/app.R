#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
install.deps("pacman")
CRAN_packages <- c("shiny",  "tidyverse", "sqldf", "dbplyr", "DT", "RSQLite", "shinydashboard")
pacman::p_load(char=CRAN_packages)

# Define UI for application that draws a histogram
ui <- dashboardPage(skin = "purple",
  dashboardHeader(title = "Shinotate - Transcriptome annotation app", titleWidth = 450),

  dashboardSidebar(sidebarMenu(
                     menuItem("Annotation", tabName = "annotation",
                              icon = icon("education", lib="glyphicon")),
                     menuItem("Expression", tabName = "de", icon = icon("bar-chart-o"))),
                   imageOutput("image")
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
                p("Use the search bar on the right of the table to retrieve transcripts matching the keyword (at any field). The table can then be copied/printed/exported using the buttons on the left. Finally, selected transcripts can be exported in FASTA format (", strong("Selection and Download"), " buttons below)."), width=12),
                br(), br(),
    
                box(div(DT::dataTableOutput("output_table"), style = "font-size:75%"),
                br(),
                div(downloadButton("download_fasta", 
                                   label = "Download sequences of selected transcripts (FASTA)"),
                    actionButton("clear_selection", label = "Clear selection"),
                    actionButton("select_all", label = "Select all")) , width=12)
        )
      )
    ),
    
    # Second tab content
    tabItem(tabName = "de",
            h2("Differential expression tab content")
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
  output$output_table <- DT::renderDataTable(mantle_annot,  rownames= F, # filter = "top",
                               extensions = c("Buttons", "Scroller"),
                               options = list(
                                 dom = 'Bfrtip',
                                 autoWidth=FALSE,
                                 buttons =list('copy','print', 
                                               list( extend = 'collection',
                                       buttons = c('csv', 'excel', 'pdf'),
                                       text = 'Download'),
                                       list(extend = 'colvis', columns = 2:8)), 
                                 deferRender = TRUE,scrollY = 200,
                                      scroller = TRUE,  lengthMenu=20
                                 
                                 ))
  # set the table as proxy to be able to update it
  proxy = dataTableProxy('output_table')
  # Identify button click 
  observeEvent(input$clear_selection, {
    proxy %>% selectRows(list())
  })
  
  observeEvent(input$select_all, {
    proxy %>% selectRows(as.numeric(input$output_table_rows_all))
  })
  
  output$download_fasta <- downloadHandler(
    filename = "selected_sequences.fasta",
    content = function(file) {
      # input_output_table_rows_selected = 1:10
      ids <- mantle_annot$TrinityId[input$output_table_rows_selected]
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
      seqinr::write.fasta(as.list(seqs_out$sequence), names = seqs_out$transcript_id, file.out = file)
    }
  )
  
}

# the indices of the selected rows are available through input$tableId_rows_selected
# Run the application 
shinyApp(ui = ui, server = server)

