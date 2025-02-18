library(shiny)
library(ape)

# Define the UI
ui <- fluidPage(
  titlePanel("DNA Trees"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "fasta_file",
        label = "Choose an aligned FASTA file:",
        choices = list(
          "FiveSpecies" = "FiveSpecies_Simple.fasta",
          "SpikeProteins" = "AlignedSpikeProteins.fa"
        )
      ),
      actionButton("generate", "Generate Tree"),
      helpText(
        "Start with the Five Species, then choose the Spike Protein"
      )
    ),
    
    mainPanel(
      strong("Identify the Relationships"),
      p("The DNA sequencea data allows us to visualise relationships among samples"),
      p("There are two types of changes among species or individuals.  First, the letter can change at the same position"),
      p("Second, there can be changes at more than one position"),
      hr(),
      code("Five Species"),
      p("Note the single change of letter between AACAA, AATAA and AAAAA; then note how AACAA and AACGA are different.  Finally, note how AACAA and AATTA are the most different"),
      hr(),
      code("Spike Protein"),
      p("Note how all the Covid-19 sequences are identical. Then note how the SARS and MERS are different from Covid; Note also that the Bats are different too"),
      plotOutput("tree_plot", height = "600px"), # Adjusted height for larger trees
      verbatimTextOutput("status")
    )
  )
)

# Define the server
server <- function(input, output, session) {
  observeEvent(input$generate, {
    # Ensure a valid file is selected
    req(input$fasta_file)
    
    # File path in the 'www' directory
    file_path <- file.path("www", input$fasta_file)
    
    # Read the aligned sequences and compute the tree
    tryCatch({
      aligned_sequences <- read.dna(file_path, format = "fasta")
      dist_matrix <- dist.dna(aligned_sequences, model = "JC69")
      phylo_tree <- nj(dist_matrix)
      
      # Render the tree plot
      output$tree_plot <- renderPlot({
        plot(phylo_tree, main = paste("Phylogenetic Tree for", input$fasta_file), cex = 1.2)
      })
      
      # Display status
      output$status <- renderText("Tree generated successfully!")
    }, error = function(e) {
      output$status <- renderText(paste("Error:", e$message))
      output$tree_plot <- renderPlot(NULL)  # Clear the plot in case of error
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)