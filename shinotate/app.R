#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
CRAN_packages <- c("shiny", "data.table", "tidyverse", "sqldf", "dbplyr", "DT")
install.deps(CRAN_packages)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Trinotate Transcriptome Database Query App" ),
   
   # Sidebar with a search bar 
   sidebarLayout(
        sidebarPanel(width=2,
          img(src = "1200px-Pinctada_margaritifera_MHNT.CON.2002.893.jpg", height = 151, width = 200),
          br(),
          br(),
          textInput("keyword", h3("Keyword search"), 
                    value = "...")),
      
      # Show a results table of transcripts matching the search
      mainPanel(
        
        h2("Pearl Oyster (", em("P. maxima"), ") Trinotate transcriptome database"), 
        p("Mantle tissues from pearl oysters were collected, RNA was extracted and was subjected to RNA-Sequencing."),
        p("The resulting data was ", em("de novo"), " assembled into a reference transcriptome, annotated and stored in a ",
          a("Trinotate", href="http://trinotate.github.io/") ," database"),
        p("Use the search bar on the left panel to retrieve transcripts whose BLAST annotation matches the keyword."),
        br(), br(),
        DT::dataTableOutput("output_table")
      )
   )
)

# Define server logic required to show a results table of transcripts matching the search
server <- function(input, output) {
  mantle_db <- src_sqlite("/mnt/bioinfo03/data/Paspaley/Mantle/Trinity_assembly/Annotation/Trinotate_db/P_maxima_mantle_Trinotate_local.sqlite")
  # retrieve blast and pfam results (remove P_maxima_Mantle... to match the accessions in the databases)
  
  mantle_annot <- tbl(mantle_db, sql("SELECT TrinityId, Description, BitScore, Pfam.HMMERDomains, Pfam.DomainDescriptions, H.GO_ids, H.GO_names, K.KO_ID AS KEGG_id, K.KO_name AS KEGG_name  FROM Blast_tax B LEFT OUTER JOIN (SELECT QueryProtID,GROUP_CONCAT(HMMERDomain, ';')  HMMERDomains ,GROUP_CONCAT(HMMERTDomainDescription, ';')  DomainDescriptions FROM HMMERDbase WHERE FullSeqScore>20 AND FullDomainScore>20 GROUP BY QueryProtID) Pfam ON B.TrinityId=Pfam.QueryProtID LEFT OUTER JOIN (SELECT * FROM HMMERDbase JOIN (SELECT pfam_acc, group_concat(name) AS GO_names, group_concat(go_id) AS GO_IDs FROM (SELECT * FROM pfam2go JOIN go ON pfam2go.go_id=go.id) GROUP BY pfam_acc) G ON substr(HMMERDbase.pfam_id,1,7)=G.pfam_acc WHERE FullSeqScore>20 AND FullDomainScore>20 GROUP BY QueryProtID HAVING MAX(FullSeqScore) AND MAX(FullDomainScore) ) H ON B.TrinityId=H.QueryProtID LEFT OUTER JOIN (SELECT * FROM ORF_KO WHERE KO_id != 'None') K ON B.TrinityId=K.orf_id GROUP BY TrinityId HAVING MAX(BitScore)")) %>% 
    collect(n=Inf) # %>% dplyr::filter(Description %likeci% "Leucine") 
  # 
  res <- reactive({mantle_annot %>% filter(grepl(input$keyword, Description, ignore.case = TRUE))})
  output$output_table <- DT::renderDataTable(res(), filter = "top", rownames= F,
                               extensions = list("Buttons" = NULL,
                                                 "FixedColumns" = list(leftColumns=1)),
                               options = list(
                                 dom = 'Brltpi',
                                 autoWidth=TRUE,
                                 buttons =list('copy','print', 
                                               list( extend = 'collection',
                                       buttons = c('csv', 'excel', 'pdf'),
                                       text = 'Download')
                                     ),  deferRender = TRUE,scrollY = 200,scroller = TRUE
                                 ))
  
}

# Run the application 
shinyApp(ui = ui, server = server)

