######--------------- UI-----------------------#####

ui <- dashboardPage(
  
  dashboardHeader(
    title = "Microbiome Studio"
  ),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      
      menuItem("File Upload", tabName = "file_upload", icon = icon("upload")),
      menuItem("ANCOMBC Analysis", tabName = "ancombc", icon = icon("chart-line")),
      menuItem("Compositional Analysis", tabName = "composition", icon = icon("chart-pie")),
      menuItem("Alpha Diversity", tabName = "alpha", icon = icon("chart-area")),
      menuItem("Beta Diversity & PERMANOVA", tabName = "beta", icon = icon("project-diagram"))
    )
  ),
  
  dashboardBody(
    theme = bs_theme(preset = "cerulean"),
    
    tabItems(
      
      ## ================= FILE UPLOAD =================
      tabItem(
        tabName = "file_upload",
        
        fluidRow(
          box(
            width = 3,
            h5("Upload META, OTU and TAXA Files (.csv Format)"),
            fileInput("meta_file", "Choose the META File", accept = c("text/csv", ".csv")),
            fileInput("otu_file", "Choose the OTU File", accept = c("text/csv", ".csv")),
            fileInput("tax_file", "Choose the TAXA File", accept = c("text/csv", ".csv")),
            actionButton("process", "Upload Files"),
            hr(),
            selectInput(
              "preview_choice",
              "Select which Table to Display",
              choices = c("Meta Table", "Otu Table", "Taxonomy Table")
            )
          ),
          
          box(
            width = 9,
            tabsetPanel(
              tabPanel("File Preview",
                       DTOutput("preview_table")
              )
            )
          )
        )
      ),
      
      ## ================= ANCOMBC =================
      tabItem(
        tabName = "ancombc",
        
        fluidRow(
          box(
            width = 3,
            h5("AncomBC Analysis for Differential Abundance"),
            selectInput("tax_level", "Select Taxonomy Level", choices = NULL),
            selectInput(
              "formula",
              label = tooltip(
                trigger = list("Select Formula (Group)", bs_icon("info-circle")),
                "Select which variable you want to compare in the log fold change plot and dataframe",
                placement = "bottom"
              ),
              choices = NULL
            ),
            actionButton("analyze", "Run ANCOMBC"),
            hr(),
            selectInput(
              "result_choice",
              "Select ANCOMBC Result to Display",
              choices = c(
                "Log Fold Changes",
                "Standard Errors",
                "Test Statistics",
                "P-values",
                "Adjusted P-values",
                "Differentially Abundant Taxa",
                "Bias-corrected Abundances"
              )
            )
          ),
          
          box(
            width = 9,
            tabsetPanel(
              tabPanel("ANCOMBC Results",
                       DTOutput("ancombc_table")),
              tabPanel("Log Fold Change Data",
                       DTOutput("lfc_table")),
              tabPanel("Log Fold Change Plot",
                       plotOutput("FoldChange"))
            )
          )
        )
      ),
      
      ## ================= COMPOSITIONAL =================
      tabItem(
        tabName = "composition",
        
        fluidRow(
          box(
            width = 3,
            h5("Venn Diagram"),
            selectInput("meta_column", "Select Grouping Variable", choices = NULL),
            uiOutput("group_select_ui"),
            actionButton("generate_venn", "Generate Diagram"),
            hr(),
            h5("Select Parameters for Compositional Analysis"),
            selectInput(
              "plot_choice",
              "Select Plot:",
              choices = c(
                "Core Composition Heatmap",
                "Order Bar Plot",
                "Pie Charts",
                "Relative Abundance Plot by Patient",
                "Relative Abundance Divided by Groups"
              )
            ),
            selectInput("comp_group", "Select Grouping Variable", choices = NULL),
            selectInput(
              "tax_level",
              "Taxonomic Level:",
              choices = c("Phylum", "Class", "Order", "Family", "Genus")
            ),
            sliderInput("detection", "Detection Threshold (%):", 1, 0, 100, 5),
            sliderInput("prevalence", "Prevalence (%):", 50, 0, 100, 5)
          ),
          
          box(
            width = 9,
            tabsetPanel(
              tabPanel("Venn Diagram",
                       plotOutput("venn_diagram")),
              tabPanel("Compositional Analysis",
                       plotOutput("selected_plot"))
            )
          )
        )
      ),
      
      ## ================= ALPHA =================
      tabItem(
        tabName = "alpha",
        
        fluidRow(
          box(
            width = 3,
            h5("Alpha Diversity Analysis"),
            selectInput("group_selection_alpha", "Select Grouping Variable", choices = NULL),
            hr(),
            actionButton("alpha_process", "Run")
          ),
          
          box(
            width = 9,
            tabsetPanel(
              tabPanel(
                "Alpha Diversity Plots",
                plotOutput("alpha_plots", brush = "plot_brush"),
                DTOutput("alpha_statistics")
              ),
              tabPanel(
                "Alpha Diversity Indexes",
                DTOutput("alpha_index")
              )
            )
          )
        )
      ),
      
      ## ================= BETA =================
      tabItem(
        tabName = "beta",
        
        fluidRow(
          box(
            width = 3,
            h5("PCoA"),
            selectInput("group_selection_beta", "Select Grouping Variable", choices = NULL),
            hr(),
            radioButtons(
              "beta_index",
              "Select The Beta Diversity Metrics",
              choices = c(
                "Bray-Curtis Dissimilarity" = "bray",
                "Jaccard Distance" = "jaccard"
              )
            ),
            hr(),
            h5("Dispersion Analysis"),
            selectInput(
              "disp_plot_sel",
              "Select which Plot to Display",
              choices = c(
                "Distance to Centroid",
                "Ordination Centroids and Dispersion Labeled"
              )
            )
          ),
          
          box(
            width = 9,
            tabsetPanel(
              tabPanel("PCoA Plot",
                       plotOutput("pcoa_plot")),
              tabPanel("PERMANOVA",
                       verbatimTextOutput("permanova_out")),
              tabPanel(
                "Dispersion Analysis",
                plotOutput("disp_plot"),
                verbatimTextOutput("dispr"),
                verbatimTextOutput("permutest")
              )
            )
          )
        )
      )
    )
  )
)