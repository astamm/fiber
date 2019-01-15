# ________________________________________________________________________________________________
# Position information


# Barycenter:

# NEW 22/11/18
get_barycenter_streamline <- function(streamline, validate = TRUE) {
  # Given a streamline, it returns the 3D coordinates of the barycenter of the curve
  if (validate) {
    if (!is_streamline(streamline)) 
      stop("The input dataset is not of class streamline.")
  }
  streamline_coords = select(streamline, x, y, z)
  bar=colMeans(streamline_coords)
  bar=data_frame(x_Bar=bar[1], y_Bar=bar[2], z_Bar=bar[3])
  return (bar)
}

# NEW 22/11/18
get_barycenter <- function (streamlist) {
  # Creates the feature Barycenter of a give tract, whose streamlines are collected in data
  
  # bar=NULL
  # for (i in 1:length(streamlist)) {
  #   bar=rbind(bar, get_barycenter_streamline(streamlist[[i]]))
  # }
  bar = map_df(streamlist, get_barycenter_streamline)
  return (bar)
}



# Spatial Median:

# NEW 22/11/18
get_sp_median_streamline = function(streamline, validate = TRUE) {
  # Given a streamline, it returns the 3D coordinates of the spatial median of the curve
  if (validate) {
    if (!is_streamline(streamline)) 
      stop("The input dataset is not of class streamline.")
  }
  streamline_coords = select(streamline, x, y, z)
  spatial_med = Gmedian(streamline_coords)
  spatial_med = data_frame(x_SpMed=spatial_med[1], y_SpMed=spatial_med[2], z_SpMed=spatial_med[3])
  return(spatial_med)
}

# NEW 22/11/18
get_spatial_median = function (streamlist){
  # Creates the feature Spatial Median of a give tract, whose streamlines are collected in data
  # sp_median = NULL
  # for (i in 1:length(streamlist)){
  #   sp_median = rbind(sp_median,get_sp_median_streamline(streamlist[[i]][,2:4]))
  # }
  # sp_median = data.frame(x_SpMed = sp_median[,1],y_SpMed = sp_median[,2],z_SpMed = sp_median[,3] )
  sp_median = map_df(streamlist, get_sp_median_streamline)
  return (sp_median)
}


# Distance spatial_median - barycenter

# distance_centers <- function(data, barycenter, spatial_median){
#   # Returns the distance between the barycenters and the spatial medians of the streamlines
#   # of a given tract
#   dist = NULL
#   for (i in 1:length(data))
#   {
#     deltaX = (barycenter[i,1]-spatial_median$x_SpMed[i])^2
#     deltaY = (barycenter[i,2]-spatial_median$y_SpMed[i])^2
#     deltaZ = (barycenter[i,3]-spatial_median$z_SpMed[i])^2
#     dist = c(dist, sqrt(deltaX+deltaY+deltaZ))
#   }
#   return (dist)
# }






# MD_info = list (n_sector, c(interior break points)) 
# intervals: [break_points[i], break_points[i+1])
create_dataset_new = function(tract, side_label, MD_info, RD_info, AD_info, FA_info, standardized = T)
{
  features = map_df(tract$Streamlines, compute_streamline_features, MD_info, RD_info, AD_info, FA_info)
  
  n = length(tract$Streamlines)
  patient = rep (as.numeric(tract$PatientId), n)
  side = rep(side_label,n)
  
  if (standardized == T){
    features = scale(features)
    features = as.data.frame(features)
  }
  
  features = mutate(features, side, patient)
  
  return(features)
}


compute_streamline_features = function (streamline, MD_info, RD_info, AD_info, FA_info) {
  # curv_torsNOSTRA = get_curvature_torsion_fda3D_lambda_fixed(streamline, lambda=lambda_opt)
  # Curvature
  curvature = get_curvature(streamline)$curvature
  curv_max = max(curvature)      
  curv_mean = mean(curvature)
  curv_sd = sd(curvature)
  # curv_maxNOSTRA = curv_torsNOSTRA[1]     
  # curv_meanNOSTRA = curv_torsNOSTRA[2]  
  # curv_sdNOSTRA = curv_torsNOSTRA[3]  
  # curv_maxDIFF= curv_max- curv_torsNOSTRA[1]
  # curv_meanDIFF = curv_mean - curv_torsNOSTRA[2]
  # curv_sdDIFF = curv_sd - curv_torsNOSTRA[3]
  
  # Torsion
  torsion = get_torsion(streamline)$torsion
  tors_max = max(torsion)      
  tors_mean = mean(torsion)
  tors_sd = sd(torsion)
  # tors_maxNOSTRA = curv_torsNOSTRA[4]     
  # tors_meanNOSTRA = curv_torsNOSTRA[5]  
  # tors_sdNOSTRA = curv_torsNOSTRA[6] 
  # tors_maxDIFF = tors_max - curv_torsNOSTRA[4]
  # tors_meanDIFF = tors_mean - curv_torsNOSTRA[5]
  # tors_sdDIFF = tors_sd - curv_torsNOSTRA[6]
  
  clength = get_curvilinear_length(streamline, validate = F)
  elength = get_euclidean_length(streamline, validate = F)
  sinuosity = clength/elength
  
  barycenter = get_barycenter_streamline(streamline)
  x_barycenter = barycenter[1]
  y_barycenter = barycenter[2]
  z_barycenter = barycenter[3]
  spatial_median = get_sp_median_streamline(streamline)
  x_SpMed = spatial_median[1]
  y_SpMed = spatial_median[2]
  z_SpMed = spatial_median[3]
  dist_cent = norm_vec2(t(barycenter-spatial_median))
  
  MD_sectors = diffusion_index_sectors(streamline, "md", MD_info)
  RD_sectors = diffusion_index_sectors(streamline, "rd", RD_info)
  AD_sectors = diffusion_index_sectors(streamline, "ad", AD_info)
  FA_sectors = diffusion_index_sectors(streamline, "fa", FA_info)
  
  
  return(data.frame(curv_max, curv_mean, curv_sd, 
                    # curv_maxNOSTRA, curv_meanNOSTRA, curv_sdNOSTRA,
                    tors_max, tors_mean, tors_sd, 
                    # tors_maxNOSTRA, tors_meanNOSTRA, tors_sdNOSTRA,
                    # curv_maxDIFF, curv_meanDIFF, curv_sdDIFF,
                    # tors_maxDIFF, tors_meanDIFF, tors_sdDIFF,
                    clength, elength, sinuosity,
                    x_barycenter, y_barycenter, z_barycenter,
                    x_SpMed,  y_SpMed,  z_SpMed, 
                    dist_cent,
                    FA_sectors,
                    MD_sectors, 
                    RD_sectors,
                    AD_sectors
  )
  )
}


diffusion_index_sectors = function(streamline, name_variable, info_diffusion_index){
  diffusion_vector = as_vector(streamline[name_variable])
  n_sectors = info_diffusion_index[[1]]
  if(n_sectors == 1) return (mean(diffusion_vector))
  break_points = info_diffusion_index[[2]]
  start = 1
  end = break_points[1]-1
  sectors = c()
  sectors = cbind(sectors,mean(diffusion_vector[start:end]))
  colnames(sectors)[1] = paste(name_variable,"sector", 1, sep = "_")
  
  i=1;
  while(i < (n_sectors-1)){
    start = break_points[i] 
    end = break_points[i+1]-1
    sectors = cbind(sectors, mean(diffusion_vector[start:end]))
    colnames(sectors)[i+1] = paste(name_variable,"sector",i+1, sep = "_")
    i=i+1
  }
  sectors = cbind(sectors, mean(diffusion_vector[break_points[n_sectors-1]:length(diffusion_vector)]))
  colnames(sectors)[n_sectors] = paste(name_variable,"sector",n_sectors, sep = "_")
  
  return (sectors)
}


norm_vec2 <- function(x){sqrt(crossprod(x))}
