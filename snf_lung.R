###### this script will run snf over lung function data
library(SNFtool)
library(tidyverse)
library(readr)
library(cluster)

# read in data 
temp_dat <- read.csv('lung_data.csv', 
                     na.strings = c("777", NA), 
                     stringsAsFactors = FALSE)

# determine which variables we can throw out - keep only the Z-score FEVs
# get the variables that have a z score transformation 
cols_temp <- colnames(temp_dat)[grepl('LCI|FEV|FVC', colnames(temp_dat))]

# only keep zscores
cols_remove <- cols_temp[!grepl('zscore', cols_temp)]

# remove those variables from the data set
temp_dat <- temp_dat[, !colnames(temp_dat) %in% cols_remove]

# remove StudyCenter column
temp_dat$StudyCenter <- NULL

# make column names lowerase 
colnames(temp_dat) <- tolower(colnames(temp_dat))

# recode factor variables to deal with blank space and NA
temp_dat$lbf_3m_status <- ifelse(temp_dat$lbf_3m_status == '', NA,
                                 ifelse(temp_dat$lbf_3m_status == 'unknown', NA, temp_dat$lbf_3m_status))

temp_dat$lbf_6m_status <- ifelse(temp_dat$lbf_6m_status == '', NA,
                                 ifelse(temp_dat$lbf_6m_status == 'unknown', NA, temp_dat$lbf_6m_status))

# determine the amount of missingness for each column
temp_cols <-as.data.frame(unlist(apply(temp_dat, 2, function(x)
  {
  round(length(which(is.na(x)))/nrow(temp_dat), 2)
  })))

temp_cols$variable <- rownames(temp_cols)
colnames(temp_cols)[1] <- 'percent_missing'
rownames(temp_cols) <- NULL


 # save data and manually enter data groups
write_csv(temp_cols, 'variable_percent_missing.csv')

# assign data type (Wall) to each variable 
temp_dat_cols <- names(temp_dat)

#remove subject type and move subject number (ID) to rownames, then remove
rownames(temp_dat) <- temp_dat$subjectnumber
temp_dat$subjectnumber <- temp_dat$subjecttype <- NULL

# look for time series variables and impute only if an NA is sandwhiched between two time 
# series variables with same value
# possible variables: (adclinic1y, adclinicc3y, lbf_3m_status, lbf_6m_status, smk_inside_3m, smk_inside_6m,
# 'smk_inside_12m, smk_inside_18m, smk_inside_24m, smk_inside_30m, smk_inside_36m)
temp <- as.data.frame(cbind(n2 = temp_dat$n2 , sf6 = temp_dat$sf6, three = temp_dat$smk_inside_3m, 
                            six = temp_dat$smk_inside_6m, twelve=temp_dat$smk_inside_12m,
                            eighteen = temp_dat$smk_inside_18, twenty_four = temp_dat$smk_inside_24m, 
                            thirty = temp_dat$smk_inside_30m, thirty_six =temp_dat$smk_inside_36))

# HERE- remember if the factor is only one level (like gender) then maybe you can just code it as 0 and 1 numeric
# HERE- anytime there is an ER visit, does that emply hospitalization? 
# HERE - how many layers can SNF handle
# create group vector to merge to column names
group_vector <- c('group_1_fac', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 
                  'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1_fac', 'group_1_fac', 
                  'group_2', 'group_2', 'group_2', 'group_2', 'group_2', 'group_2_fac', 'group_2', 'group_2', 
                  'group_2_fac', 'group_2_fac', 
                  'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 
                  'group_3', 'group_3',  'group_3', 'group_3', 'group_3', 'group_3', 
                  'group_4_fac', 'group_4_fac', 'group_4_fac', 'group_4_fac', 'group_4_fac', 'group_4_fac', 
                  'group_4_fac', 'group_4_fac', 
                  'group_5', 'group_5', 'group_5', 'group_5', 'group_5', 'group_5', 'group_5', 
                  'group_6','group_6', 'group_6', 'group_6', 'group_6',  'group_6', 'group_6', 'group_6', 'group_6',
                  'group_7', 'group_7', 'group_7', 
                  'group_8', 'group_8', 
                  'group_9', 'group_9', 'group_9',  'group_9', 'group_9', 'group_9', 'group_9', 'group_9', 
                  'group_9',  'group_9', 'group_9', 'group_9')

# add groupvector to names
names(temp_dat) <- paste0(names(temp_dat), '_', group_vector)

# 


# plot the distribution of NA
# ggplot(temp_cols, 
#        aes(reorder(variable, 
#                    -percent_missing), 
#            percent_missing)) +
#   geom_bar(stat = 'identity', 
#            fill = 'red', 
#            alpha = 0.4) +
#   labs(x = 'Variables', 
#        y = 'Percent missing') +
#   theme_gray(base_size = 10) +
#   theme(axis.text.x = element_text(angle=90)) 
#   
# 
# hist(temp_cols$percent_missing, 
#      xlab = 'Percent NA', 
#      main = 'Distribution of NAs', 
#      col = 'lightblue',
#      breaks = 15)
# 
# # start by removing variables with over 50%  NA
# temp_complete <- temp_dat[, -which(colMeans(is.na(temp_dat)) > 0.5)]
# 
# temp_full <- temp_complete[complete.cases(temp_complete), ]

###### 
# implement SNF
######

# create a list of data frames based on the types of data 
# maybe observable clinical data and numerical data

data <- temp_full
numClus = 4
sampleRows = TRUE


SNFClustering <- function(data, 
                          numClus = 4, 
                          sampleRows = TRUE,
                          dist_type,
                          multiple_data_types) {
  # Transpose data for the dist function
  # The dist function takes distances between rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- lapply(data, t)
  }
  
  # Calculate the distance between samples
  if(dist_type == 'normal') {
    distances <- as.matrix(dist(data))
  } else {
    distances <- as.matrix(daisy(data))
  }
  # Convert the distances to affinities
  affinities <- affinityMatrix(distances)
  
  if(multiple_data_types){
    # Fuse the affinity matrices
    fusedMatrix <- SNF(affinities)
    
    # Cluster the fused matrix
    labels <- spectralClustering(fusedMatrix, numClus)
  } else {
    # Cluster the fused matrix
    labels <- spectralClustering(affinities, numClus)
  }
  
  # add a column to data for the labels
  data$clustering_label <- labels
  
  return(data)
}

temp <- SNFClustering(data = temp_full, 
                      numClus = 4, 
                      sampleRows = TRUE, 
                      dist_type = 'daisy', 
                      multiple_data_types = FALSE)

