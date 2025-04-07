
if (!requireNamespace("SGSEA", quietly = TRUE) || utils::packageVersion("SGSEA") < "2.0.0") {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ShellsheDeng/SGSEA", force = TRUE)
}


library(shiny)
library(waiter)

library(limma)
library(fgsea)
library(DT)



##################### Define UI  ###################

ui <- fluidPage(
  # APP TITLE
  titlePanel(HTML("<center>Survival-based Gene Set Enrichment Analysis App (SGSEA)</center>")),

  # UI layout
  sidebarLayout(
    # Sidebar for guiding steps
    sidebarPanel(
      # NEW: Download Example Dataset Button
      downloadButton("download_example_data", "Download Example Dataset (KIRC)"),
      helpText("Click to download the example dataset (KIRC) with gene expression and survival data."),

      # NEW: Download Example R Script Button
      downloadButton("download_example_script", "Download Example Script"),
      helpText("Click to download an R script that demonstrates how to use the SGSEA package."),

      # step 1 upload gene file
      fileInput(inputId = "upload_file_gene", label = h5("Step 1: Upload data containing Gene expression and survival information")),
      helpText("The first column must be 'ID', the second 'survtime', and the third 'status'. All other columns should contain numeric gene expression values with gene symbols as column names. Click 'Input Data' (top-right panel) to verify the upload before proceeding."),
      radioButtons(inputId = "data_type_gene", label = h5("Select File Type"),
                   choices = c(".csv", ".txt", ".xlsx"), selected = ".csv"),
      helpText("The uploaded file should match the format of the KIRC dataset, including 'ID', 'survtime', 'status', and gene expression columns.
      The status indicator, normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). 
               For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored."),

      # step 2 filtering
      radioButtons("filtering", "Step 2: Apply Filtering?", choices = c("Yes" = 1, "No" = 2)),

      # Step 3 Normalization
      radioButtons("normalization", "Step 3: Apply Normalization?", choices = c("Yes" = 1, "No" = 2)),

      # Action button to start analysis
      waiter::use_waiter(),
      actionButton("gobutton", "Go!"),

      # Help & Support Section
      h5("Results will appear in the tabs on the right."),
      h5("Need Help?"),
      helpText("Contact: shellshe_deng@hotmail.com"),


    ),

    # main panel for output display
    mainPanel(
      tabsetPanel(
        tabPanel("Input Data", DT::dataTableOutput("output_kirc_data")),
        tabPanel("Filtering", DT::dataTableOutput("output_filtering1"), textOutput("output_filtering2")),
        tabPanel("Normalization", DT::dataTableOutput("output_normalization1"), textOutput("output_normalization2")),
        tabPanel("SGSEA Results", DT::dataTableOutput("output_file_sgsea")),
        tabPanel("Top10 Significant Pathways", plotOutput("output_sgsea_top10", width = "80%"))
      )
    )
  ),
  h4("License"),
  helpText("This software is released under the GPLv3 license."),
  tags$a("View License", href="LICENSE.md", target="_blank")

)


################ Define server  ##################


server <- function(input, output) {
  options(shiny.maxRequestSize = 100*1024^2)


  # NEW: Function to Download Example Dataset (KIRC)
  data("KIRC", package = "SGSEA")

  # Fallback to example if no upload
  kirc_data_reac <- reactive({
    req(input$upload_file_gene)
    file <- input$upload_file_gene$datapath
    type <- input$data_type_gene

    if (type == ".xlsx") {
      readxl::read_xlsx(file)
    } else {
      data.table::fread(file, header = TRUE, data.table = FALSE)
    }
  })

  output$download_example_data <- downloadHandler(
    filename = function() { "KIRC_example.csv" },
    content = function(file) {
      load("data/KIRC.rda")  # Loads KIRC object
      write.csv(KIRC, file, row.names = FALSE)}
  )

  output$download_example_script <- downloadHandler(
    filename = function() { "SGSEA_example_script.R" },
    content = function(file) {
      script_path <- system.file("scripts", "SGSEA_example_script.R", package = "SGSEA")
      file.copy(script_path, file)
    }
  )


  # step1 gene data output
  # if not choosing any file, show blank; otherwise, show data table
  output$output_kirc_data <- DT::renderDataTable({
    req(kirc_data_reac())
    DT::datatable(kirc_data_reac()[1:50, 1:10], options = list(lengthMenu = c(5, 30), pageLength = 8))
  })


  # step2 filtering
  # drop genes having count < (# of samples)/10, column is gene
  filtering_reac <- reactive({
    dat <- kirc_data_reac()
    if (input$filtering == 1) {
      gene_data <- dat[, -(1:3)]
      keep_genes <- colSums(gene_data) > (nrow(dat)/10)
      filtered <- cbind(dat[, 1:3], gene_data[, keep_genes])
    } else {
      filtered <- dat
    }
    filtered
  })
  # if needs filtering, show data
  output$output_filtering1 <- DT::renderDataTable({
    req(input$filtering)
    if (input$filtering == 1) {
      DT::datatable(filtering_reac()[1:50, 1:10], options = list(lengthMenu = c(5, 30), pageLength = 8))
    }
  })
  # if dont need filtering, show guidance text
  output$output_filtering2 <- renderText({
    if (input$filtering == 2) "Please proceed to normalization."
  })




  # step 2 normalization
  normalization_reac <- reactive({
    w <- waiter::Waiter$new(html = div(style = "color:green;", waiter::spin_3(), h3("Normalizing data...")))
    w$show(); on.exit(w$hide())
    Sys.sleep(1)

    dat <- filtering_reac()
    if (input$normalization == 1) {
      gene_data <- dat[, -(1:3)]
      gene_data <- gene_data[, sapply(gene_data, is.numeric)]  # keep only numeric columns

      validate(
        need(ncol(gene_data) > 0, "No numeric gene expression columns found for normalization.")
      )

      voom_data <- limma::voom(gene_data)
      normed <- cbind(dat[, 1:3], as.data.frame(voom_data$E))
    } else {
      normed <- dat
    }
    normed
  })
  # if needs normalization, show data
  output$output_normalization1 <- DT::renderDataTable({
    if (input$normalization == 1) {
      DT::datatable(round(normalization_reac()[1:50, 4:13], 4), options = list(lengthMenu = c(5, 30), pageLength = 8))
    }
  })
  # if dont need normalization, show guidance text

  output$output_normalization2 <- renderText({
    if (input$normalization == 2) "Click 'GO!' to start analysis."
  })




  # step 3 get log hazard ratio
  # computing lhr with progress bar once the user hits go button

  lhr_reac <- eventReactive(input$gobutton, {
    w <- waiter::Waiter$new(html = div(style = "color:red;", waiter::spin_3(), h3("Calculating LHR...")))
    w$show(); on.exit(w$hide())
    Sys.sleep(1)
    dat <- normalization_reac()
    survTime <- as.numeric(dat[, 2])
    survStatus <- as.numeric(dat[, 3])
    gene_expr <- dat[, 4:ncol(dat)]
    
    # Normalize before LHR
    lhr <- SGSEA::getLHR(normalizedData = gene_expr, survTime = survTime, survStatus = survStatus, covariates = NULL)
    #na.omit(lhr)
  })



  #step 4 get pathways from reactome db with progress bar once the user hits go button

  getReactome_reac <- eventReactive(input$gobutton, {
    w <- waiter::Waiter$new(html = div(style = "color:orange;", waiter::spin_3(), h3("Loading Reactome...")))
    w$show(); on.exit(w$hide())
    Sys.sleep(1)
    #readRDS(system.file("extdata", "spathways.rds", package = "SGSEA"))
   
    readRDS("www/spathways.rds")
    


  })

  

  
  

  # step 5 enrichment tabl3
  sgsea_result_reac <- eventReactive(input$gobutton, {
    w <- waiter::Waiter$new(html = div(style = "color:green;", waiter::spin_3(), h3("Running SGSEA...")))
    w$show(); on.exit(w$hide())
    Sys.sleep(1) #allows R to temporarily be given very low priority and not to interfere with more important foreground tasks
    fgsea::fgsea(pathways = getReactome_reac(), stats = lhr_reac(), minSize = 5, maxSize = 500)
  })


  # show the results data or re-run if user clicks go button
  output$output_file_sgsea <- DT::renderDataTable({
    result <- sgsea_result_reac()
    result <- result[, c("pathway", "pval", "padj", "NES", "size", "ES", "leadingEdge")]
    result$NES <- round(result$NES, 4)
    DT::datatable(result, filter = 'top', options = list(lengthMenu = c(5, 30), pageLength = 8))
  })

  sgsea_top10_reac <- eventReactive(input$gobutton, {
    SGSEA::getTop10(sgsea_result_reac(), getReactome_reac(), lhr_reac(), 0.15)
  })

  output$output_sgsea_top10 <- renderPlot({
    sgsea_top10_reac()
  })
}

shinyApp(ui = ui, server = server)


# deploy
# the R PROJECT  that opens app.R should be within the same file(same getwd())
#library(rsconnect)
#rsconnect::deployApp("/Users/shellshe/Desktop/Research/SGSEAPaper/sgseaShiny/SGSEAapp")


