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
library(ggfortify)
library(purrr)
library(reshape)
library(prettyR)
library(parallelDist)
library(RANN)



all_data <- as.data.frame(load(file = "Data/21600-0002-Data.rda"))


# preprocess factors into numerics
num_attrib <- length(da21600.0002)

for (i in 1:num_attrib) {
  if(is.factor(da21600.0002[[i]])){
    lbls <- sort(levels(da21600.0002[[i]]))
    lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
    da21600.0002[[i]] <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", da21600.0002[[i]]))
    da21600.0002[[i]] <- add.value.labels(da21600.0002[[i]], lbls)
  }
  # TODO some kind of normalization (if enough time for this)
}

# Define UI for Add health Cure clustering
ui <- dashboardPage(
  dashboardHeader(title = "Add Health - CURE"),
  dashboardSidebar(
    sidebarMenu(
      sliderInput(inputId = "k", label = "Number Cluster", min = 2, max = 20, value = 3, step = 1, round = TRUE),
      sliderInput(inputId = "alpha", label = "alpha", min = 0, max = 1, value = 0.5),
      sliderInput(inputId = "p", label = "Number of partitions", min = 0, max = 50, value = 5),
      sliderInput(inputId = "f", label = "f", min = 0, max = 1, value = 0.5),
      sliderInput(inputId = "delta", label = "delta", min = 0, max = 1, value = 0.5),
      sliderInput(inputId = "q", label = "Cluster per partition", min = 0, max = 20, value = 7)
    )
  ),
  dashboardBody(
    fluidRow(
      # output table
      plotOutput(outputId = "clusterPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #Clustering part of the CURE algorithm (without sampling or partitioning)
  CURE_cluster <- function(dataset, #data
                           k) #number clusters
  {
    dataset$cluster <- 1:nrow(dataset) # each point is it's own cluster int he beginning
    dataset$rep <- TRUE #every cluster it it's own representive point
    dataset$closest <- NA # closest cluster to the cluster each point is in
    dataset$dist <- NA # distance to closest cluster
    addedCols <- c("cluster", "rep", "closest", "dist") # list of columns we added, and need to remove before we compute distances etc.
    
    dataset <- dataset[1:20 ,] #reduce size for debug
  
    nClusters <- nrow(table(dataset$cluster)) #number of current clusters
    while(nClusters > k)
    {
      #compute for each cluster it's closest cluster
      #for each representative point of a cluster, find nearest point of another cluster
      #we need to do this first for each cluster before we can merge the two closest clusters
      for(i in 1:length(unique(dataset$cluster)))
      {
        currentCluster <- unique(dataset$cluster)[i]
        currentRep <- dataset[dataset$cluster == currentCluster & dataset$rep == TRUE,] # representative points of the current cluster
        
        #compute closest cluster only if no closest cluster known. For example because the closest cluster was merged recently
        if(!is.na(currentRep[1, "closest"]))
        {
          next
        }
        
        pointsNotInCurrentCluster <- dataset[!(dataset$cluster == currentCluster) & dataset$rep == TRUE,]
        
        #use kd-tree to find nearest point for each representative point
        #TODO not sure if Kd-Tree really is a good idea here. It seems to me that the tree is build anew every time and I didn't find a way to keep it
        nearest <- nn2(data = pointsNotInCurrentCluster[, !names(dataset) %in% addedCols],
                       query = currentRep[, !names(dataset) %in% addedCols],
                       k = 1)
        
        #extract closest point
        idxClosestPoint <- nearest$nn.idx[ which(nearest$nn.dist %in% c(min(nearest$nn.dist))) ]
        closestPoint <- pointsNotInCurrentCluster[idxClosestPoint,]
        
        #set closest cluster and distance
        dataset[dataset$cluster == currentCluster, "closest"] <- closestPoint$cluster #closest cluster to this cluster
        dataset[dataset$cluster == currentCluster, "dist"] <- min(nearest$nn.dist) # distance to closest cluster
      }
      
      #merge the two closest clusters
      clustersToMerge <- dataset[dataset$dist == min(dataset$dist),][1,][names(dataset) %in% c("cluster", "closest")]
      dataset[dataset$cluster == clustersToMerge$closest, "cluster"] <- clustersToMerge$cluster #add points of closest cluster to cluster
      dataset[dataset$closest == clustersToMerge$closest, "closest"] <- NA #invalidate closest cluster of clusters that had the cluster that was merged as closest cluster 
      dataset[dataset$cluster == clustersToMerge$cluster, "closest"] <- NA #invalidate closest cluster of the new merged cluster
    
      
      #TODO: compute new representative points
      
    
      
      nClusters <- nrow(table(dataset$cluster)) #update number of current clusters
     # nClusters <- k #remove this later, just for debugging
    }
    
  
    
    View(dataset)
  }
  
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
    
    # Cluster points in each partition into N/(pq) clusters
    ### TODO: do clustering for all partitions, not only for the first one.
    part_length <- length(partitions[[1]][[1]])
    df_part <- data.frame(reps = list(part_length), dist = numeric(part_length), closest_cluster = numeric(part_length), stringsAsFactors = FALSE)
    #dist_mat <- matrix(data = NA, nrow = part_length, ncol = length(partitions[[1]]))
    for (i in 1:part_length) {
      df_part$reps[[i]] <- as.vector(partitions[[1]][i,])
      
    }
    dist_mat <- sapply(partitions[[1]], function(x){return(x)})
    
    #View(dist_mat)
    d <- parallelDist(dist_mat)
    # todo check structure of d, to get distances between points
    
    
    #   start with each point as an own cluster 
    #   compute closest cluster to each cluster
    #   store in ascending order
    
    # Cluster previously found clusters until k clusters remain
    
    # assign remaining points that were not sampled to nearest cluster
    
  }
  
  # convert to matrix
  data_mat  <-  matrix(unlist(da21600.0002), nrow=length(unlist(da21600.0002[1])))
  # remove first col, which is just an ID
  data_mat <- data_mat[,-1]
  df <- as.data.frame(data_mat)
  df[is.na(df)] <- 0 # TODO check if NA <- 0 makes sense
  
  #Clustering <- CURE(df, 3, 0.3, 4, 0.2, 0.3, 7)
  CURE_cluster(df, 5)
   
   output$clusterPlot <- renderPlot({
     # dimension reduction
     pca <- prcomp(df)
     
     #plotting
     ggplot(data = pca, mapping = aes(x = pca$x[,1], y = pca$x[,2])) +
       geom_point(data = df, color = 1, size = 1) # TODO color based on clustering result, and size is changed for representatives
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

