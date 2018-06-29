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
library(cluster)



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
      sliderInput(inputId = "k", label = "Number Cluster", min = 2, max = 20, value = 6, step = 1, round = TRUE),
      sliderInput(inputId = "alpha", label = "alpha", min = 0, max = 1, value = 0.6),
      sliderInput(inputId = "p", label = "Number of partitions", min = 1, max = 50, value = 20),
      sliderInput(inputId = "f", label = "f", min = 0, max = 1, value = 0.1),
      sliderInput(inputId = "delta", label = "delta", min = 0, max = 1, value = 0.7),
      sliderInput(inputId = "q", label = "q", min = 0, max = 20, value = 5)
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
                           k, #number of clusters
                           numberRep, # number of representative points per cluster,
                           alpha, #factor (0 - 1)
                           centers
  ) 
  {
    if(k < 1)
    {
      k <- 1
    }
    #other variables we will use
    addedCols <- c("cluster", "rep", "closest", "dist") # list of columns we added, and need to remove before we compute distances etc.
    
    #if centers don't exist yet, create it. In the beginning each point is it's own center
    if(!is.data.frame(centers))
    {
      centers <- dataset[, !names(dataset) %in% c("rep", "closest", "dist")] #centroids of each cluster
    }
    
    handle_outliers <- TRUE
    outlier_fraction <- 1/3 #when only (outlier_fraction*initial clusters) clusters are left, handle outliers
    
    #start clustering
    nClusters <- nrow(table(dataset$cluster)) #number of current clusters
    while(nClusters > k)
    {
      
      #outlier handling #1
      if(handle_outliers & nClusters < outlier_fraction * nrow(dataset))
      {
        handle_outliers <- FALSE
        count <- as.data.frame(table(dataset$cluster))
        clustersToRemove <- count[count[,2] < 2, 1]
        
        dataset <- dataset[!dataset$cluster %in% clustersToRemove,] #remove found outliers
        centers <- centers[!centers$cluster %in% clustersToRemove,] #remove centroids of found outliers
        
        nClusters <- nrow(table(dataset$cluster)) #number of current clusters
        next #we might have < k clusters now
      }
      
      
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
        nearest <- nn2(data = pointsNotInCurrentCluster[, !names(dataset) %in% addedCols],
                       query = currentRep[, !names(dataset) %in% addedCols],
                       k = 1)
        
        #extract closest point
        idxClosestPoint <- nearest$nn.idx[ which(nearest$nn.dist %in% c(min(nearest$nn.dist))) ]
        closestPoint <- pointsNotInCurrentCluster[idxClosestPoint,][1,]
        
        #set closest cluster and distance
        dataset[dataset$cluster == currentCluster, "closest"] <- closestPoint$cluster #closest cluster to this cluster
        dataset[dataset$cluster == currentCluster, "dist"] <- min(nearest$nn.dist) # distance to closest cluster
      }
      
      #clusters that should be merged
      clustersToMerge <- dataset[dataset$dist == min(dataset$dist),][1,][names(dataset) %in% c("cluster", "closest")]
      
      #before we merge the two clusters we have to move the representative points of the two clusters away from their respective centroid
      #this is necessary because after the merge we only use the old representative points of both clusters to compute new reresentative points as proposed in the paper
      #the new representative points of the merged cluster will then be moved towards the new centroid
      idx <- which(dataset$cluster == clustersToMerge$cluster & dataset$rep == TRUE) #for the first cluster
      for(i in idx)
      {
        p <- dataset[i, !names(dataset) %in% addedCols]
        centroid <- centers[centers$cluster == clustersToMerge$cluster, !names(centers) %in% c("cluster")]
        newP <- p - alpha * (centroid - p)
        dataset[i, !names(dataset) %in% addedCols] <- newP
      }
      idx <- which(dataset$cluster == clustersToMerge$closest & dataset$rep == TRUE) #for the second cluster
      for(i in idx)
      {
        p <- dataset[i, !names(dataset) %in% addedCols]
        centroid <- centers[centers$cluster == clustersToMerge$closest, !names(centers) %in% c("cluster")]
        newP <- p - alpha * (centroid - p)
        dataset[i, !names(dataset) %in% addedCols] <- newP
      }
      
      #merge the two closest clusters
      dataset[dataset$cluster == clustersToMerge$closest, "cluster"] <- clustersToMerge$cluster #add points of closest cluster to cluster
      dataset[dataset$closest == clustersToMerge$closest, "closest"] <- NA #invalidate closest cluster of clusters that had the cluster that was merged as closest cluster 
      dataset[dataset$cluster == clustersToMerge$cluster, "closest"] <- NA #invalidate closest cluster of the new merged cluster
      centers <- centers[!(centers$cluster == clustersToMerge$closest),] #remove centroids of clusters that don't exist anymore
      
      #compute new representative points for new merged cluster
      oldRepPoints <- dataset[dataset$cluster == clustersToMerge$cluster & dataset$rep == TRUE,]
      rowIdxOfOldRepPoints <- which(dataset$cluster == clustersToMerge$cluster & dataset$rep == TRUE)
      
      #compute centroid
      centroid <- as.data.frame(colMeans(dataset[dataset$cluster == clustersToMerge$cluster, !names(dataset) %in% addedCols]))
      #we have to reshape the centroid data. The centroid should be representated by a single row, as all the other points
      colnames(centroid) <- c("value")
      centroid$attribute <- rownames(centroid)
      newCentroid <- tidyr::spread(data = centroid, attribute, value)
      centroid <- newCentroid[,centroid$attribute] #reorder columns, because spread orders them alphabetically
      centers[centers$cluster == clustersToMerge$cluster, !names(centers)  %in% c("cluster")] <- centroid #update centroid in centers
      
      #don't compute new representative points if all points in the cluster will end up as representative points
      if(nrow(oldRepPoints) > numberRep) 
      {
        dataset[dataset$cluster == clustersToMerge$cluster, "rep"] <- FALSE #invalidate old representative points
        
        #compute distances
        oldRepPoints <- rbind(oldRepPoints[, !names(dataset) %in% addedCols], centroid) #add centroid because we need distance of all points to centroid 
        d <- parallelDist(as.matrix(oldRepPoints), diag = TRUE, upper = TRUE)
        distMat <- as.data.frame(as.matrix(d))
        oldRepPoints <- oldRepPoints[-nrow(oldRepPoints),] #remove centroid
        
        #find new representative points. Mostly done as in the pseudocode of the paper
        newRepPointsIdx <- NULL
        for(i in 1:numberRep)
        {
          maxDist <- 0
          maxPointIdx <- NULL
          for(j in 1:nrow(oldRepPoints))
          {
            if(i == 1)
            {
              #first new representative point is farthest from centroid
              distancesToCentroid <- distMat[nrow(distMat), -ncol(distMat)] #distances from all old representative points to centroid
              newRepPointsIdx <- c(which(max(distancesToCentroid) == distancesToCentroid)[1])
              distMat <- distMat[-nrow(distMat),-ncol(distMat)] #remove distances to centroid
              break
            }
            else
            {
              distToRepPoints <- distMat[j, newRepPointsIdx] #only distances to the other newly chosen representative points
              minDist <- min(distToRepPoints)
            }
            if(minDist >= maxDist)
            {
              maxDist <- minDist
              maxPointIdx <- j
            }
          }
          newRepPointsIdx <- c(newRepPointsIdx, maxPointIdx) #add newly found representative point
        }
        
        #set all newly found representative points as represetative points in the dataset
        dataset[rowIdxOfOldRepPoints[newRepPointsIdx], "rep"] <- TRUE 
      
      }
      
      #shrink representative points towards centroid
      idx <- which(dataset$cluster == clustersToMerge$cluster & dataset$rep == TRUE)
      for(i in idx)
      {
        p <- dataset[i, !names(dataset) %in% addedCols]
        newP <- p + alpha * (centroid - p)
        dataset[i, !names(dataset) %in% addedCols] <- newP
      }
      
      nClusters <- nrow(table(dataset$cluster)) #update number of current clusters
    }
    returnVal <- list()
    returnVal$dataset <- dataset
    returnVal$centers <- centers
    return(returnVal)
  }
  
  ### CURE algorithm
  CURE <- function(dataset,  # data
                   k,     # number clusters (>2)
                   num_reps, # number of representative points per cluster
                   alpha, # factor (0 - 1)
                   p,     # number of partitions (>1)
                   f,     # sampling fraction (0 - 1)
                   delta, # 1-delta is the probability of sampling at least f*100% points of each cluster (0-1)
                   q      # number of clusters that should be found in a partition -> number of clusters equals 1/q of the original partition size (>1)
                   )
  {
    N <- length(dataset[[1]]) # TODO adjust to data
    #print(N)
    Ni <- N/k #TODO evaluate if this function could fit (probably Ni has to be smaller)
    inv_delta <- 1/delta
    sample_size <- f*N + N/Ni * log(inv_delta) + N/Ni * sqrt(log(inv_delta)*log(inv_delta) + 2*f*Ni*log(inv_delta))
    
    # take a random sample of size n from a dataset
    sampleSet <- sample(1:nrow(dataset), as.integer(sample_size), replace=FALSE)
    initial_sampleSet <- sampleSet
    mysample <- dataset[sampleSet,]
    
    #adding some columns to sample
    mysample$cluster <- 1:nrow(mysample) # each point is it's own cluster int he beginning
    mysample$rep <- TRUE #every cluster it it's own representive point
    mysample$closest <- NA # closest cluster to the cluster each point is in
    mysample$dist <- NA # distance to closest cluster
    addedCols <- c("cluster", "rep", "closest", "dist") # list of columns we added, and need to remove before we compute distances etc.
    
    
    # split sample into p equally sized partitions
    partitionsSize <- as.integer(length(sampleSet) / p)
    
    #print(partitionsSize)
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
    
    #cluster each partition
    clusters_per_partition <- partitionsSize / q
    partitions_combined <- NULL
    centers_combined <- NULL
    for(i in 1:length(partitions))
    {
      print(sprintf("processing partition: %d", i))
      result <- CURE_cluster(dataset = partitions[[i]],
                                k = clusters_per_partition,
                                numberRep = num_reps,
                                alpha = alpha,
                                centers = NA)
      
      centers_combined <- rbind(centers_combined, result[["centers"]])
      partitions_combined <- rbind(partitions_combined, result[["dataset"]])
    }
    
    #outlier handling #2
    #remove small clusters
    threshold <- 6 #not sure how to choose this value
    count <- as.data.frame(table(partitions_combined$cluster))
    clustersToRemove <- count[count[,2] < threshold, 1]
    
    partitions_combined <- partitions_combined[!partitions_combined$cluster %in% clustersToRemove,] #remove found outliers
    centers_combined <- centers_combined[!centers_combined$cluster %in% clustersToRemove,] #remove centroids of found outliers
    
    #combine partitions
    partitions_combined[,"closest"] <- NA #forget every closest cluster
    print("clustering partitions ...")
    result <- CURE_cluster(partitions_combined,
                                    k,
                                    num_reps,
                                    alpha,
                                    centers_combined)
    partitions_combined <- result[["dataset"]]
    centers_combined <- result[["centers"]]
    
    #we still have to assign each of the remaining points that were not in the initial sample to the clusters
    print("labeling unsampled data ...")
    dataset <- dataset[-initial_sampleSet,] #points that were not assigned to any cluster yet
  
    #for each point that was not in the initial sample search closest representative point
    haystack <- partitions_combined[partitions_combined$rep == TRUE,]
    nearest <- nn2(data = haystack[, !names(partitions_combined) %in% addedCols],
                     query = dataset,
                     k = 1)
    
    #extract closest representative points and assign the remaining points to the cluster the closest point belongs to
    closestPoints <- haystack[nearest$nn.idx,]
    dataset$cluster <- closestPoints$cluster
    
    #all points are assigned to a cluster now. 
    #Since We reduced the dataset to points that were not in the inital sample earlier, 
    #we now have to re-add the points that were in the initial sample to the dataset
    #for that we first need to add the remaining columns to the dataset
    dataset$rep <- FALSE
    dataset$closest <- NA
    dataset$dist <- NA
    
    #then we can add the points that were in the initial sample to the dataset
    dataset <- rbind(dataset, partitions_combined)
    
    #reshape for final output
    dataset$closest <- NULL
    dataset$dist <- NULL
    
    returnVal <- list()
    returnVal$data <- dataset
    returnVal$cluster <- dataset$cluster
    returnVal$centers <- centers_combined
    
    returnVal$silinfo <- silhouette(returnVal$cluster, dist(dataset[, !names(dataset) %in% addedCols]))[, 3]
    
    return(returnVal)
    
  }
  
  # convert to matrix
  data_mat  <-  matrix(unlist(da21600.0002), nrow=length(unlist(da21600.0002[1])))
  # remove first col, which is just an ID
  data_mat <- data_mat[,-1]
  df <- as.data.frame(data_mat)
  df[is.na(df)] <- 0 # TODO check if NA <- 0 makes sense
  
  #            data, num_cluster, num_reps, alpha, num_partition, sampling fraction, delta, q 
   output$clusterPlot <- renderPlot({
     clusters <- CURE(dataset = df,
                      k = input$k,
                      num_reps = 10,
                      alpha = input$alpha,
                      p = input$p,
                      f = input$f,
                      delta = input$delta,
                      q = input$q)
     
     
     # dimension reduction
     cluster_data <- clusters$data
     drops <- c("rep","cluster")
     cluster_data <- cluster_data[, !(names(cluster_data) %in% drops)]
     pca <- prcomp(cluster_data)
     
     #plotting
     ggplot(data = pca, mapping = aes(x = pca$x[,1], y = pca$x[,2])) +
       geom_point(aes(color = as.factor(clusters$cluster),
                      size = as.factor(map_dbl(
                            clusters$data$rep, 
                            function(x){ifelse(x == TRUE, return(5), return(2))}))
                      )
                 )
                  
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

