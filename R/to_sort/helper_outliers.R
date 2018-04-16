library(fields)

v <- matrix(runif(6000), nrow=2000, ncol=3)
microbenchmark::microbenchmark(
  a <- dist(v),
  b <- fields::rdist(v, compact = TRUE),
  c <- fields::rdist(v)
)

##############################################################################################
########################### DISTANCE - BASED ####################################

# This function returns the indexes of streamline identified as outliers using a distance-based approach.
# In particular for every streamline, the barycenter is projected on the xy plane (and then on the yz).
# On this plane(s) outliers are identified using the basic distance-based approach (Knorr and Ng 1997)
# Then the same procedure is repeated using the spatial median of the streamlines
get_outliers_distance <- function (tract)
{
  percentage=0.01

  ##################### SPATIAL MEDIAN (sp)
  epsilon_sp=5
  dataset_Median = get_spatial_median(tract$data)
  dataset_MedianXY = dataset_Median[,-3]
  dataset_MedianYZ = dataset_Median[,-1]

  ### XY
  # There could be a problem with memory while creating this matrix
  # distance.matrix_XY <- as.matrix(dist(dataset_MedianXY, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
  distance.matrix_sp_XY <- dist(dataset_MedianXY)
  n=dim(dataset_MedianXY)[1]

  card_neighbourhood_sp_xy = apply(as.matrix(distance.matrix_sp_XY), 1, check_neighbourhood, epsilon_sp)
  outliers.sp_xy = which(card_neighbourhood_sp_xy/(n-1) < percentage)

  ### YZ
  distance.matrix_sp_YZ <- dist(dataset_MedianYZ)
  n=dim(dataset_MedianYZ)[1]

  card_neighbourhood_sp_yz = apply(distance.matrix_sp_YZ, 1, check_neighbourhood, epsilon_sp)
  outliers.sp_yz = which(card_neighbourhood_sp_yz/(n-1) < percentage)

  #####################  BARYCENTER (bar)
  epsilon_bar=3

  dataset_Barycenter = get_barycenter(tract$data)
  dataset_BarycenterXY = dataset_Barycenter[,-3]
  dataset_BarycenterYZ = dataset_Barycenter[,-1]

  ### XY
  # distance.matrix_XY <- as.matrix(dist(dataset_BarycenterXY, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
  distance.matrix_bar_XY <- dist(dataset_BarycenterXY)
  n=dim(dataset_BarycenterXY)[1]


  card_neighbourhood_bar_xy = apply(distance.matrix_bar_XY, 1, check_neighbourhood, epsilon_bar)
  outliers.bar_xy = which(card_neighbourhood_bar_xy/(n-1) < percentage)


  ### YZ

  distance.matrix_bar_YZ <- dist(dataset_BarycenterYZ)
  n=dim(dataset_BarycenterYZ)[1]

  card_neighbourhood_bar_yz = apply(distance.matrix_bar_YZ, 1, check_neighbourhood,epsilon_bar)
  outliers.bar_yz = which(card_neighbourhood_bar_yz/(n-1) < percentage)

  outliers = unique(c(outliers.sp_xy, outliers.sp_yz, outliers.bar_yz, outliers.bar_xy))

  return (sort(outliers))
}


check_neighbourhood = function (row_dist_matrix, epsilon) {
  return(sum(row_dist_matrix<epsilon)-1)
}



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
get_outliers_depth <- function (tract){

  #### Spatial median

  dataset_Median = get_spatial_median(tract$data)
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

  dataset_Barycenter = get_barycenter(tract$data)
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

  return (sort(unique(c(outliers_depth_Barycenter, outliers_depth_Median))))

}
