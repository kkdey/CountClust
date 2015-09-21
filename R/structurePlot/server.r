library(shiny)

load("~/Dropbox/GitHub/singleCell-method/project/rdas/batch_effect_cell_cycle_genes.rda")
load("~/Dropbox/GitHub/singleCell-method/project/rdas/meta-data.rda")
nclust <- topics_list$clust2$K
omega <- topics_list$clust2$omega
phenoType <- data.frame(cellCycle = cell_phase_vector,
                        individual = 1)

choices <- setNames(1:ncol(phenoType), names(phenoType))


# Define server logic required to draw a Structure plot
shinyServer(function(input, output) {

  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot

  output$structurePlot <- renderPlot({
     # Order according to the first phenotype entry
     whichPheno <- as.numeric(input$select)
     phenoOrdered <- phenoType[ , whichPheno]

     # TBA: order samples within each phenotype group
#     phenoOrdered <- order(phenoOrdered, omega[,1])
     phenoOrdered <- order(phenoOrdered)
     omegaOrdered <- omega[phenoOrdered, ]

     # TBA: order samples with more than one phenotype

     # TBA: beautify the plot; use broman colors
     barplot(t(omegaOrdered), axisnames = FALSE,
             space = 0, border = NA, col = c(2:3))

  })
})



