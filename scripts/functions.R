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
  new_column <- as.character(do.call('rbind', new_column_results))
  return(as.character(new_column))
}



# Function to implement SNF with options for number of clusters and distance metric type
snf_clustering <- function(data, 
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

# create function that takes a vector of same type trajectories columns and imputes spline or avg on that vector of trajectories
impute_trajectory <- function(column_vector){
  # subset data with column vector
  sub_dat <- temp_dat[,names(temp_dat) %in% column_vector]
  
  # for age, weight, and height just impute the mean 
  if(any(grepl('age|weight|height', column_vector))){
    # loop through columns and get the mean to fill NA
    # sub_dat <- sub_dat[rowSums(is.na(sub_dat)) != ncol(sub_dat), ]
  
    for(j in 1:ncol(sub_dat)){
      sub_dat[is.na(sub_dat[,j]), j] <- mean(sub_dat[,j], na.rm = TRUE)
    }
  } else {
    # determine if variables are factors or numeric 
    data_type_indicator <- unique(sub_dat[,1])
    data_type_indicator <- data_type_indicator[!is.na(data_type_indicator)]
    if(length(data_type_indicator) > 5){
      # use imputation function from klm3d
      sub_dat <- longitudinalData::imputation(as.matrix(sub_dat), method="copyMean",lowerBound="globalMin",upperBound="globalMax")
    } else {
      # categorical data: majority votes, impute on similarity matrix, linear regression?
      
    }
    
  }
  
  return(sub_dat)
}


# function togroup by cluster label and summarise mean results, do same thing again for factors 
summarise_cluster_groups <- function(temp_dat){
  
  temp_dat$clustering_label <- as.numeric(temp_dat$clustering_label)
  
  # get all numeric data 
  nums <- unlist(lapply(temp_dat, is.numeric))  
  temp_num <- temp_dat[, nums]
  temp_num$clustering_label <- as.factor(temp_num$clustering_label)
  temp_dat$clustering_label <- as.factor(temp_dat$clustering_label)
  
  # get all facctor data
  facs <- unlist(lapply(temp_dat, is.factor))  
  temp_facs <- temp_dat[, facs]
  
  
  # summarise numeric and factor variables and then join by clustering_label
  # get numeric varoiables 
  temp_num <- 
    temp_num %>% group_by(clustering_label) %>% summarise_all(funs(mean))
  
  # get factor variables 
  temp_facs <- 
    temp_facs %>% 
    group_by(clustering_label) %>% 
    summarise(counts = n(),
              percent_male =  
                round((sum(childgender_g1 == 
                             'M')/counts)*100, 2),
              percent_lbf_3m_exclusive =  
                round((sum(lbf_3m_status_g4 == 
                             'Exclusive')/counts)*100, 2),
              percent_lbf_3m_exclusive_after_hospital =  
                round((sum(lbf_3m_status_g4 == 
                             'Exclusive After Hospital')/counts)*100, 2),
              percent_lbf_3m_partial =  
                round((sum(lbf_3m_status_g4 == 
                             'Partial')/counts)*100, 2),
              percent_lbf_3m_zero =  
                round((sum(lbf_3m_status_g4 == 
                             'Zero')/counts)*100, 2),
              percent_mom_anemia_g4 =  
                round((sum(mom_anemia_g4 != 
                             '0')/counts)*100, 2),
              percent_prenatal_smoke_g5 =  
                round((sum(prenatal_smoke_g5 == 
                             '1')/counts)*100, 2),
              percent_smk_inside_36m_g5=  
                round((sum(smk_inside_36m_g5 == 
                             '1')/counts)*100, 2),
              percent_smk_inside_3m_g5=  
                round((sum(smk_inside_3m_g5 == 
                             '1')/counts)*100, 2),
              percent_bronchodilator_2_3y_g6 =  
                round((sum(bronchodilator_2_3y_g6 == 
                             '1')/counts)*100, 2),
              percent_wheeze3m_g6 =  
                round((sum(wheeze3m_g6 == 
                             '1')/counts)*100, 2),
              percent_wheeze3y_g6 =  
                round((sum(wheeze3y_g6 == 
                             '1')/counts)*100, 2))

  
  # join by clustering label
  temp_final <- inner_join(temp_facs, temp_num, by = 'clustering_label')
  
  return(temp_final)
}

# create a function that regresses each variable on our clustering labels. 
get_group_stats <- function(temp_dat, method){
  

  if(method == 'lm') {
    # get variable names
    var_names <- names(temp_dat)
    # remove outcom from var_names
    var_names <- var_names[var_names != 'clustering_label']
    result_list <-list()
    outcome <- temp_dat[, ncol(temp_dat)]
    
    for(v in 1:length(var_names)){
      temp_var <- temp_dat[,v]
      var_name <- names(temp_dat)[v]
      mod_dat <- as.data.frame(cbind(cluster_group = outcome, var_name = temp_var))
      results_dat <- tidy(summary(lm(cluster_group ~ var_name, data = mod_dat)))
      results_dat$var_name <- var_name
      result_dat <- results_dat[results_dat$term != '(Intercept)',]
      result_list[[v]] <- results_dat
    }
  }
  
  if(method == 'kruskal'){
    # get outcome variable
    outcome <- temp_dat[, ncol(temp_dat)]
    
    # get all numeric data 
    nums <- unlist(lapply(temp_dat, is.numeric))  
    temp_num <- temp_dat[, nums]
    
    # groups 1, 2, and 3 are trajectories, subset and get the mean HERE
    temp_traj <- temp_num[,grepl('g1|g2|g3', names(temp_num))]
    
    # change names to match group exactly for group by
    names(temp_traj) <- c('age', 'age', 'age', 'age', 'height', 'height', 'height', 'height', 'weight', 'weight', 'weight', 'weight',
                          'g2', 'g2', 'g2', 'g2', 'g2', 'g2', 'g2', 'g3', 'g3', 'g3', 'g3')
    
    unique_names <- unique(names(temp_traj))
    mean_list <- list()
    # create mean variables
    for(i in 1:length(unique_names)){
      temp_name <- unique_names[i]
      temp_sub <- temp_traj[, grepl(temp_name, names(temp_traj))]
      temp_mean <- apply(temp_sub, 1, function(x) mean(x))
      mean_list[[i]] <- temp_mean
    }
    
    # get mean data 
    mean_dat <- as.data.frame(do.call('cbind', mean_list))
    
    # rename to identify columsn
    names(mean_dat) <- c('age_18m_3y', 'height_18m_3y','weight_18m_3y', 'no2', 'fev_3m_3y')
    
    # remove age, height, weight,  no2, and fev columns in temp_num
    temp_num <- temp_num[, !grepl('g1|g2|g3', names(temp_num))]
    
    # join temp_num with mean_dat
    temp_all <- as.data.frame(cbind(temp_num, mean_dat))
    
    var_names <- names(temp_all)
    # remove outcom from var_names
    result_list <-list()
  
    for(v in 1:length(var_names)){
      temp_var <- temp_all[,v]
      var_name <- names(temp_all)[v]
      mod_dat <- as.data.frame(cbind(cluster_group = outcome, var_name = temp_var))
      results_dat <- tidy(kruskal.test(mod_dat$cluster_group, mod_dat$var_name))
      results_dat$var_name <- var_name
      result_list[[v]] <- results_dat
    }
    

  }
  if(method == 'chi_squared'){
    # get outcome variable
    outcome <- temp_dat[, ncol(temp_dat)]
    
    # get all numeric data 
    factors <- unlist(lapply(temp_dat, is.factor))  
    temp_fac <- temp_dat[, factors]
    
    # get variable names
    var_names <- names(temp_fac)
    # remove outcom from var_names
    var_names <- var_names[var_names != 'clustering_label']
    result_list <-list()

    for(v in 1:length(var_names)){
      temp_var <- temp_fac[,v]
      var_name <- names(temp_fac)[v]
      mod_dat <- as.data.frame(cbind(cluster_group = outcome, var_name = temp_var))
      results_dat <- tidy(chisq.test(mod_dat$cluster_group, mod_dat$var_name))
      results_dat$var_name <- var_name
      result_list[[v]] <- results_dat
    }
    
  }
  
  # combine results list
  final_results <- do.call('rbind', result_list)
  return(final_results)
}


