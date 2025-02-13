library(shiny)
library(ape)
library(Biostrings)
library(ggtree)

# Define UI for the app
ui <- fluidPage(
  titlePanel("Phylogenetic Tree Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("fasta_choice", 
                  "Choose a FASTA File:", 
                  choices = c("Spike Protein" = "www/AlignedSpikeProteins.fa", 
                              "Example File" = "www/example.fasta",
                              "Upload Your Own File")),
      fileInput("fasta_file", "Upload FASTA File (if selected above)", 
                accept = c(".fasta", ".fa", ".faa")),
      actionButton("generate_tree", "Generate Phylogenetic Tree"),
      hr(),
      helpText("Select a preloaded FASTA file or upload your own to calculate and visualize the phylogenetic tree.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("APE Tree Visualization",
                 plotOutput("ape_tree_plot", height = "600px")  # Increase height
        ),
        tabPanel("GGTREE Visualization",
                 plotOutput("ggtree_plot", height = "600px")  # Consistent height
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Add resource paths for embedded FASTA files
  addResourcePath("files", "www")
  
  # Reactive to determine the FASTA file path
  fasta_path <- reactive({
    if (input$fasta_choice == "Upload Your Own File") {
      req(input$fasta_file)  # Ensure a file is uploaded
      return(input$fasta_file$datapath)
    } else {
      # Use the selected preloaded FASTA file
      return(input$fasta_choice)
    }
  })
  
  # Reactive value to store the tree
  tree_data <- reactiveVal(NULL)
  
  # Generate the phylogenetic tree
  observeEvent(input$generate_tree, {
    req(fasta_path())  # Ensure a valid file path
    
    # Read the DNA sequences from the FASTA file
    sequences <- readDNAStringSet(fasta_path())
    
    # Pairwise distance matrix (using pairwise alignment with Hamming distance)
    dist_matrix <- dist.dna(as.DNAbin(sequences), model = "raw", pairwise.deletion = TRUE)
    
    # Build a phylogenetic tree using the Neighbor-Joining method
    tree <- nj(dist_matrix)
    tree_data(tree)
  })
  
  # APE Tree Plot
  output$ape_tree_plot <- renderPlot({
    req(tree_data())
    
    # Ladderize the tree for cleaner visualization
    ladderized_tree <- ladderize(tree_data())
    
    # Calculate vertical limits
    num_tips <- length(ladderized_tree$tip.label)
    vertical_space <- num_tips * 1.2
    
    # Plot the tree
    plot.phylo(ladderized_tree, 
               main = paste("Phylogenetic Tree (", input$tree_method, ")", sep = ""),
               cex = 0.8,
               y.lim = c(0, vertical_space),
               no.margin = TRUE)
  })
  
  
  # GGTREE Visualization
  output$ggtree_plot <- renderPlot({
    req(tree_data())
    
    ggtree(tree_data()) + 
      geom_tiplab() + 
      theme_tree2() + 
      ggplot2::labs(title = paste("Phylogenetic Tree (", input$tree_method, ")", sep = ""))
  })
}

# Run the application
shinyApp(ui = ui, server = server)
