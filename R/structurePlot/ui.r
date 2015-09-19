library(shiny)
# Read input file
load("~/Dropbox/GitHub/singleCell-method/project/rdas/batch_effect_cell_cycle_genes.rda")
load("~/Dropbox/GitHub/singleCell-method/project/rdas/meta-data.rda")
nclust <- topics_list$clust2$K
omega <- topics_list$clust2$omega
phenoType <- data.frame(cellCycle = cell_phase_vector,
                        individual = 1)

choices <- setNames(1:ncol(phenoType), names(phenoType))


# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("Structure plot"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(position = "left",

                # Sidebar panel widgets
                sidebarPanel(
                    # Select input file
                    fileInput("file",
                              label = h4("Admixture output object") ),

                 # Replace the following with selectInput
                 # where the range of the clusters is automatically
                 # generated once the data is imported
                 # The choice of this selectInput then goes to
                 # inform the number of clusters to be plotted

                    # Number of clusters
                    h4("Range of the number of clusters"),
                    sliderInput("Minimum",
                                label = h5("Minimum:"),
                                min = 2,
                                max = 15,
                                value = 2 ),
                    sliderInput("Maximum",
                                label = h5("Maximum:"),
                                min = 2,
                                max = 15,
                                value = 7 ),

                    # Phenotypes of interest
                    selectInput("select",
                                label = h4("Select phenotype"),
                                choices = choices,
                                selected = 1)
                ),

                # Show a plot of the generated distribution
                mainPanel(
                  plotOutput("structurePlot")
                )
  )

))
