# This script will store functions for the SNF lung capacity project

# function to combine number of hospital visis with number of ER visits into one column then remove them
combine_with_sum <- function(column_1, column_2){
  # create a list to store results 
  new_column_results <- list()
  # loop through each row and add columns together based on a few conditions
  for(i in 1:length(column_1)){
    # get the sub columns (just one element from the row)
    sub_col_1 <- column_1[i]
    sub_col_2 <- column_2[i]
    # condtion: if that both are NA, outcome is NA, but if one is NA then treat it as zero
    if(is.na(sub_col_1) & is.na(sub_col_2)){
      temp_new <- NA
    } else {
      if(is.na(sub_col_1)){
        sub_col_1 <- 0
      } 
      if(is.na(sub_col_2)){
        sub_col_2 <- 0
      }
      # add the two columns together, at this point should be no NAs
      temp_new <- sub_col_1 + sub_col_2
    }
    new_column_results[[i]] <- temp_new
  }
  new_column <- do.call('rbind', new_column_results)
  return(new_column)
}

# Function to implement SNF with options for number of clusters and distance metric type
SNFClustering <- function(data, 
                          numClus = 4, 
                          sampleRows = TRUE,
                          dist_type) {
  # Transpose data for the dist function
  # The dist function takes distances between rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- lapply(data, t)
  }
  
  # Calculate the distance between samples
  if(dist_type == 'normal') {
    distances <- lapply(data, function(x) as.matrix(dist(x)))
  } 
  if (dist_type == 'daisy') {
    distances <- lapply(data, function(x) as.matrix(daisy(x, metric = 'gower')))
  } 
  
  # Convert the distances to affinities
  affinities <- lapply(distances, affinityMatrix)
  
  # Fuse the affinity matrices
  fusedMatrix <- SNF(affinities)
  
  # Cluster the fused matrix
  labels <- spectralClustering(fusedMatrix, numClus)
  
  # cbind data to create a final data frame with corresponding cluster labels
  final_data <- do.call('cbind', data)
  
  # add a column to data for the labels
  final_data$clustering_label <- labels
  
  return(final_data)
}

# function togroup by cluster label and summarise mean results, do same thing again for factors 
summarise_cluster_groups <- function(temp_dat){
  
  # summarise numeric and factor variables and then join by clustering_label
  temp_num <- 
    temp_dat %>% group_by(clustering_label) %>% summarise_all(funs(mean))
  
  temp_num$childgender_group_1_binary <- temp_num$lbf_3m_status_group_2_factor <- 
    temp_num$lbf_6m_status_group_2_factor <- NULL
  
  temp_fac <- 
    temp_dat %>% 
    group_by(clustering_label) %>% 
    summarise(counts = n(),
              percent_male <-  
                round((sum(childgender_group_1_binary == 
                             'M')/counts)*100, 2),
              percent_lbf_3m_exclusive <-  
                round((sum(lbf_3m_status_group_2_factor == 
                             'Exclusive')/counts)*100, 2),
              percentlbf_3m_exclusive_after_hospital <-  
                round((sum(lbf_3m_status_group_2_factor == 
                             'Exclusive After Hospital')/counts)*100, 2),
              percentlbf_3m_partial <-  
                round((sum(lbf_3m_status_group_2_factor == 
                             'Partial')/counts)*100, 2),
              percentlbf_3m_zero <-  
                round((sum(lbf_3m_status_group_2_factor == 
                             'Zero')/counts)*100, 2),
              percentlbf_6m_exclusive <- 
                round((sum(lbf_6m_status_group_2_factor == 
                             'Exclusive')/counts)*100, 2),
              percentlbf_6m_exclusive_after_hospital <-  
                round((sum(lbf_6m_status_group_2_factor == 
                             'Exclusive After Hospital')/counts)*100, 2),
              percentlbf_6m_partial <- 
                round((sum(lbf_6m_status_group_2_factor == 
                             'Partial')/counts)*100, 2),
              percentlbf_6m_zero <- 
                round((sum(lbf_6m_status_group_2_factor == 
                             'Zero')/counts)*100, 2))
  
  # join by clustering label
  temp_final <- inner_join(temp_fac, temp_num, by = 'clustering_label')
  
  return(temp_final)
}

