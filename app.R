#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(dplyr)
library(DT)
library(ggplot2)
library(purrr)
library(reshape)

all_data <- as.data.frame(load(file = "Data/21600-0002-Data.rda"))

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Old Faithful Geyser Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("bins",
                     "Number of bins:",
                     min = 1,
                     max = 50,
                     value = 30)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  ### CURE algorithm
  CURE <- function(dataset,  # data
                   k,     # number clusters (>2)
                   alpha, # factor (0 - 1)
                   p,     # number of partitions (>1)
                   f,     # sampling fraction (0 - 1)
                   delta, # 1-delta is the probability of sampling at least f*100% points of each cluster (0-1)
                   q      # number of clusters that should be found in a partition -> number of clusters equals 1/q of the original partition size (>1)
                   )
  {
    N <- length(dataset[[1]]) # TODO adjust to data
    print(N)
    Ni <- N/k #TODO evaluate if this function could fit (probably Ni has to be smaller)
    inv_delta <- 1/delta
    sample_size <- f*N + N/Ni * log(inv_delta) + N/Ni * sqrt(log(inv_delta)*log(inv_delta) + 2*f*Ni*log(inv_delta))
    print(as.integer(sample_size))
    
    # take a random sample of size n from a dataset
    sampleSet <- sample(1:nrow(dataset), as.integer(sample_size), replace=FALSE)
    mysample <- dataset[sampleSet,] 
    print(length(mysample[[1]]))
    
    # split sample into p equally sized partitions
    partitionsSize <- as.integer(length(sampleSet) / p)
    print(partitionsSize)
    partitions <- list()
    for (i in 1:p) {
      if (i == p ){
        partitions[[i]] <- mysample
        mysample <- NULL
        break()
      }
      sampleSet <- sample(1:nrow(mysample), as.integer(partitionsSize), replace=FALSE)
      partitions[[i]] <- mysample[sampleSet,]
      mysample <- mysample[-sampleSet,]
    }
    #print(length(mysample[[1]]))
    #print(length(partitions[[1]][[1]]))
    #print(length(partitions[[4]][[1]]))
    
    # Cluster points in each partition into N/(pq) clusters
    
    # Cluster previously found clusters until k clusters remain
    
    # assign remaining points that were not sampled to nearest cluster
    
  }
  
  print(typeof(all_data))
  
  Clustering <- CURE(da21600.0002, 3, 0.3, 4, 0.2, 0.3, 7)
  
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

