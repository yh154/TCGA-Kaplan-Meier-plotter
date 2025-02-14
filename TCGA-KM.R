library(shiny)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(DT)
library(ggplot2)

# Define UI
ui <- fluidPage(
  titlePanel("Survival Analysis Based on Gene Expression"),
  sidebarLayout(
    sidebarPanel(
      textInput("cancer_type", "Enter Cancer Type (e.g., TCGA-UCS):", value = "TCGA-UCS"),
      textInput("gene_name", "ENSEMBL Gene ID (e.g., ENSG00000141510.18):", value = "ENSG00000141510.18"),
      actionButton("run", "Run Analysis"),
      br(), br(),
      downloadButton("download_plot", "Download Survival Plot")
    ),
    mainPanel(
      plotOutput("survival_plot"),
      DTOutput("survival_data_table")
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Reactive value to store the survival plot
  survival_plot <- reactiveVal(NULL)

  observeEvent(input$run, {
    # Get user inputs
    cancer_type <- input$cancer_type
    gene_of_interest <- input$gene_name

    # Load clinical data
    clinic <- read.table("clinical_TCGA-GDC_cleaned.txt", sep = "\t", header = TRUE, quote = '', check.names = FALSE) %>%
      filter(`Project Identifier` == cancer_type) %>%
      select(`Patient ID`, `Overall Survival (Months)`, `Overall Survival Status`)

    # Download gene expression data (RNA-Seq)
    query <- GDCquery(project = cancer_type,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

    GDCdownload(query)
    data <- GDCprepare(query)

    # Extract the expression data for the gene of interest
    gene_expression <- assay(data, "fpkm_uq_unstrand")[gene_of_interest, ]
    gene_expression <- gene_expression %>%
      as.data.frame() %>%
      mutate(`Patient ID` = data$patient,
             sample_type = data$sample_type,
             stage = data$ajcc_pathologic_stage) %>%
      dplyr::rename(expression = ".")

    # Merge clinical and gene expression data
    survival_data <- clinic %>%
      inner_join(y = gene_expression %>% filter(sample_type != "Solid Tissue Normal"), by = "Patient ID") %>%
      group_by(`Patient ID`) %>%
      mutate(mean_expression = mean(expression, na.rm = TRUE)) %>%
      ungroup() %>% unique() %>%
      mutate(median_expression = median(.[['mean_expression']])) %>%
      mutate(group = ifelse(mean_expression > median_expression, "High", "Low")) %>%
      mutate(time = `Overall Survival (Months)`, status = substr(`Overall Survival Status`, 1, 1)) %>%
      mutate_at("status", as.numeric)

    # Fit the survival model
    fit <- survfit(Surv(time, status) ~ group, data = survival_data)

    # Generate the survival plot
    plot <- ggsurvplot(fit, data = survival_data,
                       pval = TRUE,
                       conf.int = FALSE,
                       legend.labs = c("High", "Low"),
                       risk.table = TRUE, # Add risk table
                       risk.table.col = "strata", # Change risk table color by groups
                       surv.median.line = "hv", # Specify median survival
                       palette = c("#F8766D", "#2E9FDF"),
                       legend.title = cancer_type)

    # Store the plot in the reactive value
    survival_plot(plot)

    # Render the survival plot
    output$survival_plot <- renderPlot({
      print(plot)
    })

    # Render the survival data table
    output$survival_data_table <- renderDT({
      datatable(survival_data %>% select(-mean_expression,-time,-status), options = list(pageLength = 10))
    })
  })

  # Download handler for saving the survival plot
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("survival_plot_", input$cancer_type, "_", input$gene_name, ".png", sep = "")
    },
    content = function(file) {
      # Save the plot as a PNG file
      ggsave(file, plot = survival_plot()$plot, device = "png", width = 10, height = 6)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
