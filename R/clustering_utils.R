
dataset_first_points = function(cst, ind_x, ind_y, ind_z){   # cst[[1]]= cst_left   cst[[2]]=cst_right
  dataset = NULL
  nl = length(cst$lhs$Streamlines)
  nr = length(cst$rhs$Streamlines)

  first_left = map_df(cst$lhs$Streamlines, slice,1) [c(ind_x, ind_y, ind_z)]   # Prendo le colonne x,y,z della prima riga
  first_right = map_df(cst$rhs$Streamlines, slice,1) [c(ind_x, ind_y, ind_z)]

  first_left[,1] = -first_left[,1]  # Proietto il tratto sinistro nel piano delle x positive
  tmp = rbind(first_left, first_right)
  
  side = rep(c("left","right"),c(nl,nr))
  patient = rep (as.numeric(cst$lhs$PatientId), nl+nr)
  dataset = cbind(tmp,side,patient)
  
  return(dataset)
}

add_cluster_column = function (data_first_pat_j, clusters, num_streamline_patients) {
  j = data_first_pat_j$patient[1]  # Indice del paziente: 0 per il malato e gli altri a seguire
  ################################### MANCA IL MALATO ###################################
  # if (j == 0 ) {
  #   num_previous_streamline = 0
  # }
  # else {
  #   num_previous_streamline = sum(num_streamline_patients[1:j])
  # }
  # clust = clusters [(num_previous_streamline+1) : (num_previous_streamline + num_streamline_patients[(j+1)])]
  # return (cbind(data_first_pat_j, clust))
  if (j == 1 ) {
    num_previous_streamline = 0
  }
  else {
    num_previous_streamline = sum(num_streamline_patients[1:(j-1)])
  }
  clust = clusters [(num_previous_streamline+1) : (num_previous_streamline + num_streamline_patients[(j)])]
  return (cbind(data_first_pat_j, clust))
  #######################################################################################
}

get_cluster_clara = function(X,n, PLOT = F){
  
  M <- colMeans(X[,1:3])
  S = cov(X[,1:3])
  
  PC <- princomp(X[,1:3])
  set.seed(1994)
  clara.obj = clara(PC$scores[,1],n)
  
  if(PLOT == T){
    open3d()
    plot3d(X[which(clara$clustering==1),], size=3, col=rainbow(n)[1], aspect = F)
    for (i in 2:n)
    {
      points3d(X[which(clara$clustering==i),], size = 3, col = rainbow(n)[i])
    }
    
  }
  return(clara.obj)
}

get_k_opt = function(data_big, n=22, treshold = 20, num_healthy_patients){
  for(i in 1:n){
    clara = get_cluster_clara(data_big,i)
    cluster = clara$clustering
    new_data = cbind(data_big, cluster)
    # We stop when we find a healthy patient for whom a cluster has a number of elements smaller than the threshold
    for(j in 1:num_healthy_patients){
      check_sx = filter(new_data, side == "left", patient == j)$cluster
      # check_sx = new_data[which(side == "left"  &&  patient == j) ,6]
      num_sx = table(check_sx)
      check_sx = (unique(check_sx))
      # check_dx = new_data[which(side == "right" &&  patient == j) ,6]
      check_dx = filter(new_data, side == "left", patient == j)$cluster
      num_dx = table(check_dx)
      check_dx = unique(check_dx)
      if( sum(is.element(1:i,check_sx))<i ||    # Check if there is an empty cluster
          sum(is.element(1:i,check_dx))<i || 
          min(num_sx)< treshold || 
          min(num_dx)< treshold ) 
        {
          return(i-1)
        }
    }
  }
  return (n)
}

choose_k_sick = function (data, k_opt){
  b <- w <- NULL
  for(k in (k_opt-3):(k_opt+3)){
    result.k <- kmeans(data, k)
    w <- c(w, sum(result.k$wit)) # Equivalent to c(w,result.k$tot.wit)
    b <- c(b, result.k$bet)
  }
  quartz()
  # x11()
  matplot((k_opt-3):(k_opt+3), w/(w+b), pch='', xlab='clusters', ylab='within/tot', main='Choice of k', ylim=c(0,1))
  lines((k_opt-3):(k_opt+3), w/(w+b), type='b', lwd=2)
  
  quartz()
  # x11()
  matplot((k_opt-3):(k_opt+3), b/(w+b), pch='', xlab='clusters', ylab='between/tot', main='Choice of k', ylim=c(0,1))
  lines((k_opt-3):(k_opt+3), b/(w+b), type='b', lwd=2)
}

reproject_x = function (data_patient_j) {
  names_col = names (data_patient_j)
  ind_x = which(names_col=="x")
  ind_side = which(names_col=="side")
  data_patient_j = apply (data_patient_j, 1, change_sign_x_if_left, ind_x=ind_x, ind_side=ind_side)
  return (as.data.frame(t(data_patient_j)))
}

change_sign_x_if_left = function (row_data_patient_j, ind_x, ind_side) {
  if (row_data_patient_j[ind_side] == "left") 
    row_data_patient_j[ind_x] = - as.numeric(row_data_patient_j[ind_x])
  return(row_data_patient_j)
}

merge_left_right_features = function (features_patient_j) {
  return (rbind(features_patient_j[[1]], features_patient_j[[2]]))
}

get_reduced_tot = function(features, mean_left, sd_left, mean_right, sd_right){
  # features needs to have patient and cluster columns
  
  # left
  features_left = features[features$side == "left",]
  fac_left = factor(features_left$clust)
  n_left = length(levels(fac_left))
  tmp_sx = apply(features_left[,1:33], 2, tapply, fac_left, mean)
  
  tmp_var_sx = split(as.data.frame(scale(features_left[,1:33])), features_left$clust)
  var_sx = map(tmp_var_sx, cov)
  
  side_sx = rep("left", n_left)
  clust = levels(fac_left)
  patient =  rep(features_left$patient[1], n_left)
  centroids_left = data.frame(tmp_sx, side = side_sx, patient = patient,  clust = clust)
  
  # right
  features_right = features[features$side == "right",]
  fac_right = factor(features_right$clust)
  n_right = length(levels(fac_right))
  tmp_dx = apply(features_right[,1:33], 2, tapply, fac_right, mean)
  
  tmp_var_dx = split(as.data.frame(scale(features_right[,1:33])), features_right$clust)
  var_dx = map(tmp_var_dx, cov)
  
  side_dx = rep("right", n_right)
  clust = levels(fac_right)
  patient =  rep(features_right$patient[1], n_right)
  centroids_right = data.frame(tmp_dx, side = side_dx, patient = patient,  clust = clust)
  
  centroids = rbind(centroids_left, centroids_right)
  rownames(centroids) = 1:dim(centroids)[1]
  return(list(centroids = as.data.frame(centroids), var_sx = var_sx, var_dx = var_dx))
}

distance_from_centroid = function(features_row, centroid) {
  return( sqrt (sum ( ( as.numeric(features_row[1:33])- as.numeric(centroid[1:33]) )^2 ) ))
}




num_of_streamline_patient = function (cst) {
  return (num_of_streamline_tract(cst$lhs) + num_of_streamline_tract(cst$rhs))
}

num_of_streamline_tract = function (tract) {
  return (length(tract$Streamlines))
}

num_of_streamline_tract_sx_dx = function (tract_sx_dx) {  
  # tract_sx_dx[[1]]= tract_left  tract_sx_dx[[2]]= tract_right
  return (c(length(tract_sx_dx[[1]]$data), length(tract_sx_dx[[2]]$data)))
}

num_of_streamline_from_tensor_list_cst = function (tensor_list_cst) { 
  return (c(length(tensor_list_cst$lhs), length(tensor_list_cst$rhs)))
}
