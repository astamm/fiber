# library(fields)
source("../Create_features_dataset/utility_functions.R")
# NEW 22/11/18

##############################################################################################
########################### DISTANCE - BASED ####################################

# This function returns the indexes of streamline identified as outliers using a distance-based approach.
# In particular for every streamline, the barycenter is projected on the xy plane (and then on the yz). 
# On this plane(s) outliers are identified using the basic distance-based approach (Knorr and Ng 1997)
# Then the same procedure is repeated using the spatial median of the streamlines
get_outliers_distance <- function (tract_side, percentage=0.01, separate_sp_bar = F)
{
  
  
  ##################### SPATIAL MEDIAN (sp)
  epsilon_sp=5
  dataset_Median = get_spatial_median(tract_side$Streamlines)
  dataset_MedianXY = dataset_Median[,-3]
  dataset_MedianYZ = dataset_Median[,-1]
  
  ### XY
  # distance.matrix_sp_XY <- rdist(dataset_MedianXY)
  # card_neighbourhood_sp_xy = apply(distance.matrix_sp_XY, 1, check_neighbourhood, epsilon_sp)
  
  distance.matrix_sp_XY <- dist(dataset_MedianXY, method = "euclidean")
  card_neighbourhood_sp_xy = check_neighbourhood(distance.matrix_sp_XY, epsilon_sp)
  
  n=dim(dataset_MedianXY)[1]
  
  outliers.sp_xy = which(card_neighbourhood_sp_xy/(n-1) < percentage)
  
  # outliers.sp_xy <- NULL
  # points.sp_xy <- NULL
  # for (i in 1:n) {
  #   count=-1
  #   for (j in 1:n) {
  #     if (distance.matrix_XY[i,j]<epsilon) 
  #       count=count+1
  #   }
  #   
  #   if (count/(n-1) <percentage){
  #     outliers.sp_xy <- c(outliers.sp_xy, i) #Punti che hanno meno di percentage punti a distanza minore di epsilon
  #     points.sp_xy=rbind(points.sp_xy, dataset_MedianXY[i,])
  #   }
  # }
  
  ### YZ
  # distance.matrix_sp_YZ <- rdist(dataset_MedianYZ)
  # card_neighbourhood_sp_yz = apply(distance.matrix_sp_YZ, 1, check_neighbourhood, epsilon_sp)
  
  n=dim(dataset_MedianYZ)[1]
  
  distance.matrix_sp_YZ <- dist(dataset_MedianYZ, method = "euclidean")
  card_neighbourhood_sp_yz = check_neighbourhood(distance.matrix_sp_YZ, epsilon_sp)
  
  
  # outliers.sp_yz <- NULL
  # points.sp_yz <- NULL
  # for (i in 1:n) {
  #   count=-1
  #   for (j in 1:n) {
  #     if (distance.matrixYZ[i,j]<epsilon) 
  #       count=count+1
  #   }
  #   
  #   if (count/(n-1) <percentage){
  #     outliers.sp_yz <- c(outliers.sp_yz, i) #Punti che hanno meno di percentage punti a distanza minore di epsilon
  #     points.sp_yz=rbind(points.sp_yz, dataset_MedianYZ[i,])
  #   }
  # }
  
  outliers.sp_yz = which(card_neighbourhood_sp_yz/(n-1) < percentage)
  
  #####################  BARYCENTER (bar)
  epsilon_bar=3
  
  dataset_Barycenter = get_barycenter(tract_side$Streamlines)
  dataset_BarycenterXY = dataset_Barycenter[,-3]
  dataset_BarycenterYZ = dataset_Barycenter[,-1]
  
  ### XY
  # distance.matrix_XY <- as.matrix(dist(dataset_BarycenterXY, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
  # distance.matrix_bar_XY <- rdist(dataset_BarycenterXY)
  # card_neighbourhood_bar_xy = apply(distance.matrix_bar_XY, 1, check_neighbourhood, epsilon_bar)
  
  n=dim(dataset_BarycenterXY)[1]
  
  distance.matrix_bar_XY <- dist(dataset_BarycenterXY, method = "euclidean")
  card_neighbourhood_bar_xy = check_neighbourhood(distance.matrix_bar_XY, epsilon_bar)
  
  outliers.bar_xy = which(card_neighbourhood_bar_xy/(n-1) < percentage)
  
  
  ### YZ
  
  # distance.matrix_bar_YZ <- rdist(dataset_BarycenterYZ)
  # card_neighbourhood_bar_yz = apply(distance.matrix_bar_YZ, 1, check_neighbourhood,epsilon_bar)
  n=dim(dataset_BarycenterYZ)[1]
  
  distance.matrix_bar_YZ <- dist(dataset_BarycenterYZ, method = "euclidean")
  card_neighbourhood_bar_yz = check_neighbourhood(distance.matrix_bar_YZ, epsilon_bar)
  
  outliers.bar_yz = which(card_neighbourhood_bar_yz/(n-1) < percentage)
  
  if (separate_sp_bar) {
    return(list(outliers.sp_xy=outliers.sp_xy, outliers.sp_yz=outliers.sp_yz, 
                outliers.bar_yz=outliers.bar_yz, outliers.bar_xy=outliers.bar_xy))
  }
  outliers = unique(c(outliers.sp_xy, outliers.sp_yz, outliers.bar_yz, outliers.bar_xy))
  return (sort(outliers))
}


check_neighbourhood = function (dist_object, epsilon) {
  n <- attr(dist_object, "Size")
  purrr::map_int(1:n, check_neighborhood_of_i, dist_object=dist_object, epsilon=epsilon, n=n)
}

check_neighborhood_of_i <- function(i, dist_object, epsilon, n) {
  # n= attr(dist_object, "Size")
  j_indices_before= NULL
  j_indices_after= NULL
  if (i!=1)
    j_indices_before = 1:(i-1)
  if (i!=n)
    j_indices_after = (i+1):n
  distance_indices_before = n * (j_indices_before - 1) - j_indices_before * (j_indices_before - 1) / 2 + i - j_indices_before
  distance_indices_after = n * (i - 1) - i * (i - 1) / 2 + j_indices_after - i
  dvalues_before = dist_object[distance_indices_before]
  dvalues_after = dist_object[distance_indices_after]
  return(sum(dvalues_before < epsilon)+sum(dvalues_after < epsilon))
}


# check_neighbourhood = function (row_dist_matrix, epsilon) {
#   return(sum(row_dist_matrix<epsilon)-1)
# }
# 
# check_neighbourhood_dist_cycle = function (dist_object, epsilon) {
#   n= attr(dist_object, "Size")
#   card=NULL
#   for (i in 1:n) {
#     count=0
#     for (j in 1:n) {
#       if (j>i)
#         dij =  dist_object[n*(i-1) - i*(i-1)/2 + j-i]   # Distanza tra riga i e j (con j maggiore di i e minore o uguale a n)
#       else if (j<i)
#         dij =  dist_object[n*(j-1) - j*(j-1)/2 + i-j]
#       if (i!=j && dij<epsilon) 
#          count=count+1;
#     }
#     card=c(card, count)
#   }
#   return (card)
# }

##############################################################################################
########################### DEPTH - BASED ####################################
library(depth)

# Returns the index of the streamline with the barycenter closest to the point outlier
index_min_distance <- function (outlier, dataset) {
  n=dim(dataset)[1]
  dist_vec=rdist.vec(dataset,cbind(rep(outlier[1], n), rep(outlier[2], n)))
  return(which.min(dist_vec))
}

# Given a set of points (spatial median or barycenters) identified as outliers, 
# returns the indexes of the corresponding streamlines
idx_outliers <- function( dataset, outliers) 
{
  idx= apply(as.matrix(outliers), 1, index_min_distance, dataset=dataset)
  return(idx)
  
}

# This function returns the indexes of streamline identified as outliers using a depth-based approach.
# In particular for every streamline, the barycenter is projected on the xy plane (and then on the yz). 
# On this plane(s) outliers are identified using the basic distance-based approach (Knorr and Ng 1997)
# Then the same procedure is repeated using the spatial median of the streamlines
get_outliers_depth <- function (tract_side, separate_sp_bar = F){
  
  #### Spatial median
  
  dataset_Median = get_spatial_median(tract_side$Streamlines)
  dataset_MedianXY = dataset_Median[,-3]
  dataset_MedianYZ = dataset_Median[,-1]
  
  # XY
  #isodepth(dataset_MedianXY, dpth =c(1,2,3), output = FALSE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = c('red','green', 'blue', 'orange', 'purple'))
  outliers_XY = isodepth(dataset_MedianXY, dpth =c(1,2,3), output = TRUE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = NULL)
  
  outliers1_XY = as.data.frame(outliers_XY$Contour1)
  outliers2_XY = as.data.frame(outliers_XY$Contour2)
  outliers3_XY = as.data.frame(outliers_XY$Contour3)
  
  idx1 = idx_outliers(dataset_MedianXY, outliers1_XY)
  idx2 = idx_outliers(dataset_MedianXY, outliers2_XY)
  idx3 = idx_outliers(dataset_MedianXY, outliers3_XY)
  
  # YZ
  #isodepth(dataset_MedianYZ, dpth =c(1,2,3, 4), output = FALSE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = c('red','green', 'blue', 'orange', 'purple'))
  outliers_YZ = isodepth(dataset_MedianYZ, dpth =c(1,2,3, 4), output = TRUE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = NULL)
  
  outliers1_YZ = as.data.frame(outliers_YZ$Contour1)
  outliers2_YZ = as.data.frame(outliers_YZ$Contour2)
  outliers3_YZ = as.data.frame(outliers_YZ$Contour3)
  outliers4_YZ = as.data.frame(outliers_YZ$Contour4)
  
  idx4 = idx_outliers(dataset_MedianYZ, outliers1_YZ)
  idx5 = idx_outliers(dataset_MedianYZ, outliers2_YZ)
  idx6 = idx_outliers(dataset_MedianYZ, outliers3_YZ)
  idx7 = idx_outliers(dataset_MedianYZ, outliers4_YZ)
  
  outliers_depth_Median = sort(unique(c(idx1, idx2, idx4, idx5, idx6, idx7)))
  
  #### Barycenter
  
  dataset_Barycenter = get_barycenter(tract_side$Streamlines)
  dataset_BarycenterXY = dataset_Barycenter[,-3]
  dataset_BarycenterYZ = dataset_Barycenter[,-1]
  
  # XY
  #isodepth(dataset_BarycenterXY, dpth =c(1,2,3), output = FALSE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = c('red','green', 'blue', 'orange', 'purple'))
  outliers_XY = isodepth(dataset_BarycenterXY, dpth =c(1,2,3), output = TRUE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = NULL)
  
  outliers1_XY = as.data.frame(outliers_XY$Contour1)
  outliers2_XY = as.data.frame(outliers_XY$Contour2)
  outliers3_XY = as.data.frame(outliers_XY$Contour3)
  
  idx1 = idx_outliers(dataset_BarycenterXY, outliers1_XY)
  idx2 = idx_outliers(dataset_BarycenterXY, outliers2_XY)
  idx3 = idx_outliers(dataset_BarycenterXY, outliers3_XY)
  
  # YZ
  #isodepth(dataset_BarycenterYZ, dpth =c(1,2,3, 4), output = FALSE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = c('red','green', 'blue', 'orange', 'purple'))
  outliers_YZ = isodepth(dataset_BarycenterYZ, dpth =c(1,2,3, 4), output = TRUE, twodim = TRUE,mustdith = FALSE, maxdith = 10, dithfactor = 10,trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X",ylab = "Y", zlab = "Tukey's depth", colcontours = NULL)
  
  outliers1_YZ = as.data.frame(outliers_YZ$Contour1)
  outliers2_YZ = as.data.frame(outliers_YZ$Contour2)
  outliers3_YZ = as.data.frame(outliers_YZ$Contour3)
  outliers4_YZ = as.data.frame(outliers_YZ$Contour4)
  
  idx4 = idx_outliers(dataset_BarycenterYZ, outliers1_YZ)
  idx5 = idx_outliers(dataset_BarycenterYZ, outliers2_YZ)
  idx6 = idx_outliers(dataset_BarycenterYZ, outliers3_YZ)
  idx7 = idx_outliers(dataset_BarycenterYZ, outliers4_YZ)
  
  outliers_depth_Barycenter = sort(unique(c(idx1, idx2, idx4, idx5, idx6, idx7)))
  
  if(separate_sp_bar) {
    return (list(outliers_depth_Barycenter = outliers_depth_Barycenter, 
                 outliers_depth_Median = outliers_depth_Median))
  }
  
  return (sort(unique(c(outliers_depth_Barycenter, outliers_depth_Median))))
  
}

