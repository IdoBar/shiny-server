#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

devtools::source_gist('7f63547158ecdbacf31b54a58af0d1cc', filename = 'util.R')
# load needed packages. Note that for the shiny app, these need to be preinstalled with root priviliges (sudo su - -c "R -e \"install.packages('pacman', repos='http://cran.rstudio.com/')\"")
# install.deps("pacman") 

CRAN_packages <- c("shiny", "dplyr", "purrr", "sqldf", "dbplyr", "DT", "RSQLite",
                   "shinydashboard", "shinyWidgets", "Biostrings", "R.utils", 
                   "ggplot2", "plotly", "Cairo", "DESeq2") # tidyverse
pacman::p_load(char=CRAN_packages, install = FALSE)


# Find all available files
db_files <- list.files("/mnt/Shinotate_data/Trinotate_dbs", 
                       pattern = "[tT]rinotate.*\\.sqlite", 
                       full.names = TRUE)
dbs_list <- setNames(as.list(db_files), basename(db_files))
##### Specify default parameters #####
def_query <- "SELECT * FROM blast_annotation" # ORF_annotation / BlastDbase / blast_annotation
def_db <- "P_maxima_mantle_Trinotate_local.sqlite"
# def_con <- DBI::dbConnect(RSQLite::SQLite(), dbs_list[[def_db]]) # /srv/shiny-server/bioinfo03/data/Paspaley/Annotation_dbs
num_results <- 1
max_intgroups <- 2
brew_pals <- row.names(RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category=="qual",])
brew_pals <- setNames(brew_pals, brew_pals)
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

# Extract annotation summary
get_number_features <- function(con, table, cols){
  table_data <- dplyr::tbl(con, table) %>% 
    dplyr::select(one_of(cols)) %>% collect(n=Inf) 
  setNames(list(c("Total" = nrow(table_data), "Unique" = length(unique(table_data[[1]])))), cols)
  
}
# Features and the tables that holds them (for summarising)
queries <- c("transcript_id" = "Transcript", "gene_id" = "Transcript", "orf_id" = "ORF")
# ExN50 <- dplyr::tbl(def_con, "ExN50_stats") %>% dplyr::rename(E.N50=`E-N50`) %>%
#   collect() %>% mutate_at(-1, as.numeric) %>% mutate(Ex=as.numeric(gsub("E", "", .[["#E"]])))



# Define UI for the application 
ui <- dashboardPage(skin = "purple", 
  dashboardHeader(title = "Shinotate - Transcriptome annotation app", titleWidth = 450),

  dashboardSidebar(sidebarMenu(id = "sidebar",
                     menuItem("Transcript Annotation", tabName = "annotation",
                              icon = icon("list-alt", lib="glyphicon")), # icons: icon("dna", lib="font-awesome)
                     menuItem("Transcriptome Stats", tabName = "stats",
                              icon = icon("info-circle", lib="font-awesome")),
                     menuItem("Sample Information", tabName = "samples",
                              icon = icon("leaf", lib="glyphicon")),#, badgeLabel = "advanced", badgeColor = "green"),
                     menuItem("Transcript Expression", tabName = "expression", 
                              icon = icon("bar-chart-o")),
                     menuItem("Experimental Metadata", tabName = "metadata",
                              icon = icon("book", lib="glyphicon")), #, badgeLabel = "advanced", badgeColor = "green"),
                     br(),
                     # imageOutput("image"),
                     dropdownButton(
                       uiOutput('resetable_input'), size="default",
                       circle = TRUE, status = "success", icon = icon("database", lib="font-awesome"),  # cog, database 
                       width = "600px",
                       tooltip = tooltipOptions(title = "Trinotate db settings")  ),
                     
                     menuItem("Shinotate code and wiki", icon = icon("file-code-o"), 
                              href = "https://github.com/IdoBar/shiny-server/wiki")#, 
                              #badgeLabel = "info", badgeColor = "green")
                    )
                   ),
  dashboardBody(
    tags$head(
      tags$style(
      type="text/css",
      "#image img {max-width: 100%; width: 90%; height: auto}"
      ),
      # tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"
      # ),
      
      tags$style(HTML("
          /* Smaller font for preformatted text */
                      pre, table.table {
                      font-size: medium;
                      }
                      
                      body {
                      min-height: 2000px;
                      font-size: large;
                      }
                      
                      .option-group {
                      border: 1px solid #ccc;
                      border-radius: 6px;
                      padding: 0px 5px;
                      margin: 5px -10px;
                      background-color: #f5f5f5;
                      }
                      
                      .option-header {
                      color: #79d;
                      text-transform: uppercase;
                      margin-bottom: 5px;
                      }
                      "))),
    tabItems(
      # Sample info tab content
      tabItem(tabName = "samples",
              fluidRow(
                box(
                  h2("Sample information"),
                  p("Mantle tissues from pearl oysters were collected from two pearl families, differing in the average quality of their produced pearls. The mantle of each individual was dissected from 3 distinct regions: proximal, distal and central to the foot base."), width=11)) ,
              fluidRow(
                column(width=5,
                  box(#title = "Samples table",
                      div(
                        h3("Samples table"),
                        DT::dataTableOutput("samples_table")),
                      br(),
                      div(
                        actionButton("excl_selected", label = "Exclude selected samples",
                                     style = "display: inline-block !important;")
                        # actionButton("restore_excluded", "Restore sample table",  width = '150px', 
                        #              style = "display: inline-block !important;")# ,
                        # actionButton("select_all", label = "Select all")# ,
                      ) , width=NULL)),
                column(width=6,
                  box(#title = "PCA plot", 
                      h3("PCA plot"),
                      uiOutput("PCA_ui"), width=12),
                  column(width=5,
                  div(class = "option-group",
                      h3("PCA parameters"),
                  # box(title = "PCA parameters", width=NULL,
                      
                      div(style="display:inline-block",
                          awesomeRadio("feature", "Feature type",
                                       choices = c("Transcript", "Gene"), inline = TRUE, 
                                       status="primary")),
                      div(style="display: inline-block; width: 50px;",
                          HTML("<br>")),
                      div(style="display:inline-block",
                          actionBttn("recalc_pca", "Calculate PCA", icon = icon("refresh"),
                                     style = "jelly", size = "sm", color="primary")), # settings for shinywidgets::actionBttn
                          sliderInput("min_counts", "Min count to keep gene/transcript",
                                      min=1, max=100, value=10, step=1),
                      div(style="display:inline-block",
                          awesomeRadio("trans_fun", "Count transformation function",
                                             c("vst", "rlog"), selected = "vst", inline = TRUE, status="primary"
                                             # thick = TRUE, animation = "smooth", bigger=FALSE,outline=TRUE
                                             )),
                  
                      div(style="display: inline-block; width: 10px;",
                          HTML("<br>")),
                      div(style="display:inline-block",
                          materialSwitch(inputId = "blind_trans",label = "Blind", 
                                         status = "primary", inline = TRUE, value=FALSE)))),
                      # div()#)
                        # ),
                column(width=7,
                       # box(title = "Modify plot appearance", width=NULL,
                           div(class = "option-group", 
                               h3("Modify plot"),
                               uiOutput("design_selection"),
                               awesomeRadio("prin_comp", "Principal components to plot",
                                                  c("1:2", "3:4"), selected = "1:2", inline = TRUE, 
                                            status = "primary"),
                                                  # thick = TRUE, animation = "smooth", bigger=FALSE,
                                                  # outline=FALSE),
                               # radioButtons("prin_comp", "Principal components to plot",
                               #              c("1:2", "3:4"), inline = TRUE),
                               div(style="display:inline-block; width: 35%;",
                               selectInput("brew_pal", "Brewer palette", choices = brew_pals, 
                                           selected = "Set1")),
                               div(style="display: inline-block; width: 10%;",
                                   HTML("<br>")),
                               div(style="display:inline-block; width: 35%;",
                               selectInput("pca_theme", "ggplot2 theme", 
                                           choices = list("Black-White" = "bw", "Grey" = "grey", 
                                                          "Classic" = "classic", "Dark" = "dark", "Light" = "light"), 
                                           selected = "grey")), # , width = '45%'
                               sliderInput("plot_font_size", "Font base size",
                                           min=5, max=30, value=15, step=1)
                               )) # verbatimTextOutput("brush_selection")
                      # )
                  ))),
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
                box(# title = "Transcript annotation table",
                  div(h3("Transcript annotation table"),
                    DT::dataTableOutput("annot_table")),
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
                # box(
                #   h2("Transcriptome summary"),
                #   # h3("Annotation summary"),
                #   
                #    width=11),
                #   br(),
                box(width=4, #title = "Annotation summary",
                    div(h3("Annotation summary"),
                        p("Each annotation table was summarised to find the unique number of annotations for each feature across the transcriptome. The information is summarised in the table below."),
                        DT::dataTableOutput("annot_sum"))),
                box(width=6,
                  h3("ExN50 plot"),
                  p("An alternative to the Contig Nx statistic that could be considered more appropriate for transcriptome assembly data is the ExN50 statistic. It is computed as the traditional N50, but limited to the top most highly expressed transcripts that represent x% of the total normalized expression data (Ex).",
                    "Additional information on the ExN50 statistic can be found at the ", a("Trinity documentation wiki", href="https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats"), "."),
                  actionBttn("plot_ExN50", "Calculate and Plot ExN50",# icon = icon("refresh"),
                             style = "jelly", size = "sm", color="primary"),
                  # actionButton("plot_ExN50", label = "Calculate and Plot ExN50"),
                  uiOutput('ExN50_plot'))
                  )
              ),
    # Expression tab content
    tabItem(tabName = "expression",
            h2("Transcript expression")),
            # imageOutput("exp_image")),
    # Metadata (analysis pipeline) tab content
    tabItem(tabName = "metadata",
            fluidRow(
              box(
                h2("Metadata (experimental design and analysis pipeline)"),
            # includeHTML("/mnt/Shinotate_data/Assembly_stats/P_maxima_tissues_assembly_stats.html"),
            width=12)))# ,
    
    # Trinotate setup tab content (currently settings are in a dropdownButton in the sidebar)
    # tabItem(tabName = "settings") #,
            # uiOutput('resetable_input')
            
            # h2("Trinotate db setup"),
            # box(div(selectInput("db", "Trinotate db file", choices = dbs_list,
            #                     selected = dbs_list$P_maxima_mantle_Trinotate_local.sqlite)),
            #     h3("Trinotate db information:"),
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
  
  
  output$exp_image <- renderImage({
    list(
      src = "/srv/shiny-server/shinotate/www/Buccalin_KiSS-1r_genes_panel_05_06_2017.png",
      # # contentType = "image/jpg", 
      # height = 120,
      alt = "Transcript expression"
    )}, deleteFile = FALSE)
    # When input$n is 3, filename is ./images/image3.jpeg
    # filename <- normalizePath(file.path('/mnt/Shinotate_data/plots',
    #                                     'Buccalin_KiSS-1r_genes_panel_05_06_2017.png'))

  # resettable dropdownButton
  output$resetable_input <- renderUI({
    times <- input$restore_defaults
    div(id=letters[(times %% length(letters)) + 1], style = "color: black !important;",
        h2("Trinotate db setup"),    
    box(
        div(selectInput("db", "Trinotate db file", choices = dbs_list,
                        selected = dbs_list[[def_db]])),
        h3("Trinotate db information:"),
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
  
  
 ##### Define reactive and observed events/values ###### 
  # make a reactive connection to the db 
  C1  <- reactive({
    dbplyr::src_dbi(DBI::dbConnect(RSQLite::SQLite(), input$db))
  })
  
  # initialise the connection with default settings
  values <- reactiveValues(annot_conn = DBI::dbConnect(RSQLite::SQLite(), dbs_list[[def_db]]),
                    annot_table = get_annotation(dbplyr::src_dbi(DBI::dbConnect(RSQLite::SQLite(), 
                                                  dbs_list[[def_db]])), def_query),
                    samples_table = dplyr::tbl(DBI::dbConnect(RSQLite::SQLite(), dbs_list[[def_db]]), 
                                                      "samples_table") %>%  collect(n=Inf),
                           show_exn50 = FALSE)
  # Update value after user selection
  observeEvent(c(input$update_query), {
    values$annot_conn <- C1() # dbplyr::src_dbi(DBI::dbConnect(RSQLite::SQLite(), input$db))
    values$annot_table <- get_annotation(C1(), input$query)
    values$samples_table <- dplyr::tbl(C1(), "samples_table") %>% collect(n=Inf)
                                
  })
  
  
  
  # Update value
  # observeEvent(input$plot_ExN50,{
  #   values$show_exn50 <- TRUE
  # })
  
  # print db information
  
  output$db_info <- renderText({
    db_info <- R.utils::captureOutput(print(C1()))
    sprintf("%s\n%s", db_info[1], paste(db_info[-1], collapse = "\n"))
  }
  )
 
  
  
  #### Samples table ####
  output$design_selection <- renderUI({
    
    awesomeCheckboxGroup("pca_intgroups", "Select up to 2 variables to map to Colour and Shape",
                       colnames(values$samples_table), selected = c("Family", "Mantle"), inline = TRUE, 
                       status = "primary" #,thick = TRUE, animation = "smooth", bigger=FALSE, outline=TRUE
                       )
  })
  output$samples_table <- DT::renderDataTable(values$samples_table, server = FALSE,  
                                              rownames= T, # filter = "top",
                             extensions = c("Buttons", "Scroller"), # ,'FixedColumns',"FixedHeader"
                             options = list(
                               dom = 'Bfrtip',
                               responsive = TRUE,
                               # autoWidth = TRUE,
                               # columnDefs = list(list(width = '100px', targets = c(1, 5))),
                               pageLength = 12,
                               paging=TRUE,scrollY = 420,
                               scroller = TRUE,
                               buttons =list('copy','print',
                                             list( extend = 'collection',
                                                   buttons = c('csv', 'excel'), # , 'pdf'
                                                   text = 'Download')) #,
                            ) ) 
  # # Row selection in table 
  # set the table as proxy to be able to update it
  samples_proxy = dataTableProxy('samples_table', session)
  
  
  output$PCA_ui <- renderUI({ # 
    plotOutput("PCA_plot", height=400,
               click = "plot_click",
               dblclick = dblclickOpts(id = "plot_dblclick"),
               hover = hoverOpts(id = "plot_hover",nullOutside = TRUE),
               brush = brushOpts(id = "plot_brush")#,
                 # delay = input$brush_delay,
                 # delayType = input$brush_policy,
                 # direction = input$brush_dir,
                 # resetOnNew = input$brush_reset
               )
    # )
            })
  
  # Identify button click
  observeEvent(input$plot_click, {
    res <- nearPoints(PCA_data(), input$plot_click,
                      threshold = 5)
    samples_proxy %>% selectRows(which(values$samples_table[[1]] %in% res$name))
  })
  
  observeEvent(input$plot_brush, {
    res <- brushedPoints(PCA_data(), input$plot_brush)
    samples_proxy %>% selectRows(which(values$samples_table[[1]] %in% res$name))
  })
  # Remove samples from table
  observeEvent(input$excl_selected,{
    if (!is.null(input$samples_table_rows_selected)) {
      values$samples_table <- values$samples_table[-input$samples_table_rows_selected,]
    }
    # ids <- values$annot_table$TrinityID[input$output_table_rows_selected]
  })
  
  # Restore excluded samples
  # observeEvent(input$restore_excluded,{
  #   values$samples_table <- dplyr::tbl(C1(), "samples_table") %>% collect(n=Inf)
  # })
  
  
  # output$brush_selection <- renderPrint({
  #   cat("Brush selection:\n")
  #   # str(input$plot_brush)
  #   print(brushedPoints(PCA_data(), input$plot_brush))
  # })
  
  #### Calculate and plot PCA ####
  calc_PCA <- eventReactive(c(input$recalc_pca, input$excl_selected), {
                      
    feature_selection <- substr(input$feature, 1, 1) #"T"
    filter_trans_counts <- input$min_counts # 10
    blind_flag <- input$blind_trans # FALSE
    design_vars <- paste(input$pca_intgroups, collapse="+") # "Group"
    design_f <- paste("~", design_vars)
    trans_func <- get(input$trans_fun) # rlog "vst"
    withProgress(message = 'Calculating PCA: ', value = 0,{
      replicate_info <- dplyr::tbl(values$annot_conn, "Replicates") %>% dplyr::collect() 
      rep_dict <- setNames(object = replicate_info$replicate_name, nm=replicate_info$replicate_id)
      # samples_table <- dplyr::tbl(values$annot_conn, "samples_table") %>% dplyr::collect() 
      incProgress(1/4, detail = "Retrieving count matrix")
      count_table <- dplyr::tbl(values$annot_conn, "Expression") %>% 
        dplyr::filter(feature_type==feature_selection) %>% 
        dplyr::select(TrinityID=feature_name, replicate_id, frag_count) %>% 
        dplyr::collect(n=Inf) %>%
        dplyr::mutate(replicate_id=rep_dict[replicate_id], frag_count=round(frag_count)) %>% 
        tidyr::spread("replicate_id", "frag_count", fill=0) %>% # include in progress bar 
        .[c("TrinityID", values$samples_table[[1]])] # order columns as the same order of samples in samples_table
      # create the DESeq2 object using the data frames
      incProgress(2/4, detail = "Constructing dds object")
      dds <- DESeqDataSetFromMatrix(countData = count_table %>% as.data.frame() %>% 
                                      tibble::column_to_rownames("TrinityID") ,
                                    colData = values$samples_table %>% as.data.frame() %>% 
                                      tibble::column_to_rownames("Sample"),
                                    design = as.formula(design_f))
      filtered_dds <- dds[ rowSums(counts(dds)) > filter_trans_counts, ]
      incProgress(3/4, detail = "Performing count transformation")
      # Variance stabilizing transformation
      trans_dds <- trans_func(filtered_dds, blind=blind_flag) # include in progress bar 
      # Principal component
      incProgress(4/4, detail = "Calculating and plotting PCA")
      # plotPCA4(trans_dds, intgroup = input$pca_intgroups, returnData=TRUE, 
      #                      PCs = prin_comps, ntop = 100000)
      })
    return(trans_dds)
  })
  
  # Get all the user input
  PCA_data <- reactive(plotPCA4(calc_PCA(), intgroup = input$pca_intgroups, returnData=TRUE, 
                       PCs = prin_comps, ntop = 100000))

  # Make sure that a maximum of 2 intgroups can be selected
  observe({
    if(length(input$pca_intgroups) > max_intgroups){
      updateCheckboxGroupInput(session, "pca_intgroups", selected= tail(input$pca_intgroups,max_intgroups))
    }
    })
  
  output$PCA_plot <- renderPlot({
    
      # Reduce margins
      # par(mar = c(4, 4, 0.5, 0.5))
    
      # Get all the user input
      # PCA_data <- plotPCA4(calc_PCA(), intgroup = input$pca_intgroups, returnData=TRUE, 
      #                      PCs = prin_comps, ntop = 100000)
      user_theme <- get(paste0("theme_", input$pca_theme))
      percentVar <- 100 * attr(PCA_data(), "percentVar")
    
      
      
      switch (length(input$pca_intgroups)+1,
        {PCA_intgroups2 <- NULL
         PCA_intgroups1 <- NULL },
        {PCA_intgroups2 <- NULL
            PCA_intgroups1 <-  input$pca_intgroups[1] },
        {PCA_intgroups1 <-  input$pca_intgroups[1]
             PCA_intgroups2 <-  input$pca_intgroups[2]}

      )
      # 
      # PCA_intgroups1 <- case_when(length(input$pca_intgroups)==0 ~ NULL,
      #                             length(input$pca_intgroups)>0 ~ input$pca_intgroups[1])
      # 
      # PCA_intgroups2 <- case_when(length(input$pca_intgroups)==0 ~ NULL,
      #                             length(input$pca_intgroups)==1 ~ NULL,
      #                             length(input$pca_intgroups)==2 ~ input$pca_intgroups[2])
      
      
      
      brew_pal <- input$brew_pal
      prin_comps <- eval(parse(text=input$prin_comp)) # 1:2  
      # Create text for hovering (not implemented currently)
      mytext=paste0("Sample Name: ", PCA_data()$name, "\n")    
      
      # Make the plot in ggplot2
      g <- ggplot(PCA_data(), aes_string(sprintf("PC%s", prin_comps[1]), sprintf("PC%s", prin_comps[2]), 
                                       color=PCA_intgroups1, shape=PCA_intgroups2)) +
        xlab(sprintf("PC%s: %.2f%% variance", prin_comps[1], percentVar[prin_comps[1]])) + 
        geom_point(size=4) +
        ylab(sprintf("PC%s: %.2f%% variance", prin_comps[2], percentVar[prin_comps[2]])) + 
        scale_color_brewer(palette = brew_pal) +
        user_theme(input$plot_font_size) +
        coord_fixed()
      g
      # })
  })

  ##### Annotation Summary rendering  #####
  
  ##### ExN50 data and plot ####
  exn50_data <- eventReactive(input$plot_ExN50, {
    withProgress(message = 'Calculating ExN50: ', value = 0,{
      incProgress(1/4, detail = "Retrieving transcripts")
      # retrieve transcripts
      tx_data <- dplyr::tbl(values$annot_conn, "Transcript") %>% 
        dplyr::select(TrinityID=transcript_id, sequence) %>%
        dplyr::collect(n=Inf)
      incProgress(2/4, detail = "Retrieving expression data")
      # Retrieve Expression data
      exp_data <-  dplyr::tbl(values$annot_conn, "Expression") %>% 
        dplyr::filter(feature_type=="T") %>% 
        dplyr::group_by(feature_name) %>% dplyr::summarise(feature_count=sum(fpkm)) %>% 
        dplyr::filter(feature_count>0) %>% dplyr::collect(n=Inf) %>% 
        dplyr::arrange(desc(feature_count)) %>% 
        dplyr::mutate(cum_sum=cumsum(feature_count)) %>% 
        dplyr::mutate(percentile=cum_sum/sum(feature_count)*100)
      # Calculate ExN50 (see https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)
      incProgress(3/4, detail = "Calculating ExN50")
      purrr::map_df(1:100, ~ dplyr::bind_rows(calc_ExN50(., exp_data, tx_data)))
    })
  })
  
  output$ExN50_plot <- renderUI({
    # if (values$show_exn50==TRUE){
      
        E100N50 <<- exn50_data() %>% dplyr::filter(Ex==100)
        EmaxN50 <<- exn50_data() %>% dplyr::filter(ExN50==max(ExN50))
      # })
      div(plotlyOutput('ExN50Plot', height = "350px"),
          br(),
      p("The ExN50 per Ex plot above shows that the effective N50 'improved' from ", 
            strong(sprintf(" %s bp to %s bp (at E%sN50)", 
                           prettyNum(E100N50$ExN50, big.mark = ","), 
                           prettyNum(EmaxN50$ExN50, big.mark = ","),EmaxN50$Ex )),
        "compared with the N50 calculated with all transcripts (E100)."))
    # } else {
    #   p("Click on the button to calculate and show ExN50 plot")
    # }
    
  })
  
  # Render plotly 
  output$ExN50Plot <- renderPlotly({
    # Reduce margins
    old_mar <- par(mar = c(1, 2, 1, 2))
    g <- ggplot(exn50_data(), aes(x=Ex, y=ExN50)) +
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
    
    if (names(dev.cur()) != "null device") dev.off()
    pdf(NULL)
    g
    # ggplotly(g)
  })
  
  #### Annotation summary table ####
  
    #### Summarise stats #####
  output$annot_sum <- DT::renderDataTable({
    withProgress(message = 'Summarising annotations', value = 0,{
      incProgress(1/8)
    # summarise number of transcripts, genes and orfs
    sum_annot <- sapply(names(queries), 
                        function(nq) get_number_features(values$annot_conn, queries[nq], nq ), 
                        USE.NAMES = FALSE) %>% map_dbl("Unique") 
    incProgress(2/8)
    blast_tax_sum <- dplyr::tbl(values$annot_conn, "Blast_tax") %>% 
      dplyr::select(TrinityId, DatabaseSource) %>% collect(n=Inf) %>% 
      group_by(DatabaseSource) %>% summarise(Unique=n()) %>% 
      mutate(Feature=paste(c( "orf", "transcript"), DatabaseSource, "annotations", sep = "_")) %>%
      dplyr::select(Feature, Unique)
    # update progress
    incProgress(3/8)
    blastp_sum <- dplyr::tbl(values$annot_conn, "BlastDbase") %>% 
      dplyr::select(TrinityID, Feature=DatabaseSource) %>% 
      collect(n=Inf) %>%  distinct() %>% group_by(Feature) %>% summarise(Unique=n())
    # update progress
    incProgress(4/8)
    pfam_sum <- dplyr::tbl(values$annot_conn, sql("SELECT H.QueryProtID, substr(H.pfam_id,1,7) pfam_acc, G.GO_terms FROM (SELECT * FROM HMMERDbase WHERE FullSeqScore>20 AND FullDomainScore>20 GROUP BY QueryProtID HAVING MAX(FullSeqScore) AND MAX(FullDomainScore)) H JOIN (SELECT pfam_acc, group_concat(go_id) AS GO_terms FROM pfam2go GROUP BY pfam_acc) G ON substr(H.pfam_id,1,7)=G.pfam_acc")) %>%
      dplyr::select(pfam_acc) %>% distinct() %>% collect() %>%
      nrow() %>% tibble(Feature="Pfam", Unique=.)
    # update progress
    incProgress(5/8)
    go_sum <- dplyr::tbl(values$annot_conn, sql("SELECT H.QueryProtID, substr(H.pfam_id,1,7) pfam_acc, G.GO_terms FROM (SELECT * FROM HMMERDbase WHERE FullSeqScore>20 AND FullDomainScore>20 GROUP BY QueryProtID HAVING MAX(FullSeqScore) AND MAX(FullDomainScore)) H JOIN (SELECT pfam_acc, group_concat(go_id) AS GO_terms FROM pfam2go GROUP BY pfam_acc) G ON substr(H.pfam_id,1,7)=G.pfam_acc")) %>%
      dplyr::select(GO_terms) %>% collect() 
    # update progress
    incProgress(6/8)
    go_df <- tibble(Feature="GO", Unique=length(unique(unlist(go_sum$GO_terms))))
    ko_sum <- dplyr::tbl(values$annot_conn, "ORF_KO") %>% filter(KO_ID!="None") %>% 
      dplyr::select(Feature=KO_ID) %>% collect() %>% distinct() %>%
      nrow() %>% tibble(Feature="KEGG", Unique=.)
    incProgress(3/8)
    signalP_sum <- dplyr::tbl(values$annot_conn, "SignalP") %>% filter(prediction=="YES") %>% 
      dplyr::select(Feature=query_prot_id) %>% collect() %>% distinct() %>%
      nrow() %>% tibble(Feature="SignalP", Unique=.)
    # update progress
    incProgress(7/8)
    tmhmm_sum <- dplyr::tbl(values$annot_conn, "tmhmm") %>% filter(PredHel=="PredHel=1") %>% 
      dplyr::select(Feature=queryprotid) %>% collect() %>% distinct() %>%
      nrow() %>% tibble(Feature="TmHMM", Unique=.)
    })
    # Summarise all together
    tibble(Feature=names(sum_annot), Unique=sum_annot) %>% 
      bind_rows(blast_tax_sum,blastp_sum, pfam_sum,                                                             go_df, ko_sum, signalP_sum, tmhmm_sum) }, 
    rownames= FALSE, server = FALSE, 
    extensions = c("Buttons", "Scroller"), # ,'FixedColumns',"FixedHeader"
    options = list(
      dom = '<"bottom"B>rtip',
      responsive = TRUE,
      # autoWidth = TRUE,
      # columnDefs = list(list(width = '100px', targets = c(1, 5))),
      pageLength = 12,
      paging=TRUE,scrollY = 420,
      scroller = TRUE,
      buttons =list('copy','print',
                    list( extend = 'collection',
                          buttons = c('csv', 'excel'), # , 'pdf'
                          text = 'Download'))
    
  ))
  
  
  
  #### Annotation table ####
  output$annot_table <- DT::renderDataTable(values$annot_table,  server = TRUE,
                               rownames= F, # filter = "top",
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
  annot_proxy = dataTableProxy('annot_table')
  # Identify button click 
  observeEvent(input$clear_selection, {
    annot_proxy %>% selectRows(list())
  })
  
  observeEvent(input$select_all, {
    # print(as.numeric(input$output_table_rows_all))
    annot_proxy %>% selectRows(as.numeric(input$annot_table_rows_all))
  })
  
  #### Download sequences as fasta  ####
  # curr
  output$download_fasta <- downloadHandler(
    filename = "selected_sequences.fasta",
    content = function(file) {
      # input_output_table_rows_selected = 1:10
      ids <- values$annot_table$TrinityID[input$annot_table_rows_selected]
      # ids <- mantle_annot$TrinityId[input_output_table_rows_selected]
      tx_ids <- unique(sub("\\|m\\..+", "", ids))
      cds_ids <- ids[grepl("\\|m\\..+", ids)]
      # extract transcript sequences
      tx <- dplyr::tbl(values$annot_conn, "Transcript") %>% 
        dplyr::filter(transcript_id %in% tx_ids) %>% 
        dplyr::select(transcript_id, sequence) %>% dplyr::collect(n=Inf)
      # extract just the coding sequences (reverse complement if on the minus strand)
      cds <- dplyr::tbl(values$annot_conn, "ORF") %>% dplyr::filter(orf_id %in% cds_ids) %>%  
        dplyr::collect(n=Inf) %>% 
        dplyr::mutate(cds=substr(tx$sequence[match(transcript_id, tx$transcript_id)], 
                                 lend, rend), 
                      cds_2=ifelse(strand=="+", cds, chartr("ATGC","TACG",reverse(cds)))) %>% 
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

