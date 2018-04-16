create_dataset = function(tract, side_label, standardized = T)
{
  n = length(tract$data)
  side = rep(side_label,n)
  barycenter = get_barycenter(tract$data)
  spatial_median = get_spatial_median(tract$data)
  dist = distance_centers(tract$data, barycenter, spatial_median)
  
  #cat("Every 'o' is a streamline \n")
  temp = map(tract$data, get_curvature_torsion_fda3D)       
  temp_matrix = NULL
  for (i in 1: length(temp)) {
    temp_matrix=rbind(temp_matrix, temp[[i]])
  }
  
  features <- tract %>%
    as_tibble() %>%
    transmute(
      
      curv_max = temp_matrix[,1],        
      curv_mean = temp_matrix[,2],
      curv_sd = temp_matrix[,3],
      tors_max = temp_matrix[,4],
      tors_mean = temp_matrix[,5],
      tors_sd = temp_matrix[,6],
      
      clength = map_dbl(data, get_curvilinear_length, validate = F),
      elength = map_dbl(data, get_euclidean_length, validate = F),
      sinuosity = map_dbl(data, get_sinuosity, validate = F),
      x_barycenter = barycenter[,1],
      y_barycenter = barycenter[,2],
      z_barycenter = barycenter[,3],
      x_SpMed = spatial_median$x_SpMed,
      y_SpMed = spatial_median$y_SpMed,
      z_SpMed = spatial_median$z_SpMed,
      dist_cent = dist
    )
  
  MD_sectors = get_MD(tract)
  RD_sectors = get_RD(tract)
  AD_sectors = get_AD(tract)
  FA_sectors = get_FA(tract)
  
  features = mutate(features,
                    FA_sector1 = FA_sectors[,1],
                    FA_sector2 = FA_sectors[,2],
                    FA_sector3 = FA_sectors[,3],
                    FA_sector4 = FA_sectors[,4],
                    FA_sector5 = FA_sectors[,5],
                    MD_sector1 = MD_sectors[,1],
                    MD_sector2 = MD_sectors[,2],
                    MD_sector3 = MD_sectors[,3],
                    MD_sector4 = MD_sectors[,4],
                    RD_sector1 = RD_sectors[,1],
                    RD_sector2 = RD_sectors[,2],
                    RD_sector3 = RD_sectors[,3],
                    RD_sector4 = RD_sectors[,4],
                    AD_sector1 = AD_sectors[,1],
                    AD_sector2 = AD_sectors[,2],
                    AD_sector3 = AD_sectors[,3],
                    AD_sector4 = AD_sectors[,4]
  )
  
  if (standardized == T){
    features = scale(features)
    features = as.data.frame(features)
  }
  
  features = mutate(features, side)
  
  return(features)
}




MD_sector = function(streamline, sector){
  if (sector == 1)
  {
    return (mean(streamline[1:25,7]))
  }
  if (sector == 2)
  {
    return (mean(streamline[26:35,7]))
  }
  if (sector == 3)
  {
    return (mean(streamline[36:39,7]))
  }
  if (sector == 4)
  {
    return (mean(streamline[40:50,7]))
  }
  
}
get_MD <- function (tract)
{
  # Primo settore (0-25)
  sector1 = map_dbl(tract$data, MD_sector, sector = 1)
  
  # Secondo settore (26-35)
  sector2 = map_dbl(tract$data, MD_sector, sector = 2)
  
  # Terzo settore (35-39)
  sector3 = map_dbl(tract$data, MD_sector, sector = 3)
  
  # Quarto settore (40-50)
  sector4 = map_dbl(tract$data, MD_sector, sector = 4)
  
  return (cbind(sector1, sector2, sector3, sector4))
}

RD_sector = function(streamline, sector){
  if (sector == 1)
  {
    return (mean(streamline[1:16,6]))
  }
  if (sector == 2)
  {
    return (mean(streamline[17:32,6]))
  }
  if (sector == 3)
  {
    return (mean(streamline[33:36,6]))
  }
  if (sector == 4)
  {
    return (mean(streamline[37:50,6]))
  }
  
}
get_RD <- function (tract)
{
  # Primo settore (1-16)
  sector1 = map_dbl(tract$data, RD_sector, sector = 1)
  
  # Secondo settore (17-32)
  sector2 = map_dbl(tract$data, RD_sector, sector = 2)
  
  # Terzo settore (33-36)
  sector3 = map_dbl(tract$data, RD_sector, sector = 3)
  
  # Quarto settore (37-50)
  sector4 = map_dbl(tract$data, RD_sector, sector = 4)
  
  return (cbind(sector1, sector2, sector3, sector4))
}

AD_sector = function(streamline, sector){
  if (sector == 1)
  {
    return (mean(streamline[1:17,5]))
  }
  if (sector == 2)
  {
    return (mean(streamline[18:28,5]))
  }
  if (sector == 3)
  {
    return (mean(streamline[29:34,5]))
  }
  if (sector == 4)
  {
    return (mean(streamline[35:50,5]))
  }
  
}
get_AD <- function (tract)
{  # Primo settore (1-17)
  sector1 = map_dbl(tract$data, AD_sector, sector = 1)
  
  # Secondo settore (18-28)
  sector2 = map_dbl(tract$data, AD_sector, sector = 2)
  
  # Terzo settore (29-34)
  sector3 = map_dbl(tract$data, AD_sector, sector = 3)
  
  # Quarto settore (35-50)
  sector4 = map_dbl(tract$data, AD_sector, sector = 4)
  
  return (cbind(sector1, sector2, sector3, sector4))
}


FA_sector = function(streamline, sector){
  if (sector == 1)
  {
    return (mean(streamline[1:10,8]))
  }
  if (sector == 2)
  {
    return (mean(streamline[11:18,8]))
  }
  if (sector == 3)
  {
    return (mean(streamline[19:28,8]))
  }
  if (sector == 4)
  {
    return (mean(streamline[29:34,8]))
  }
  if (sector == 5)
  {
    return (mean(streamline[35:50,8]))
  }
  
}
get_FA <- function (tract)
{
  # Primo settore (1-10)
  sector1 = map_dbl(tract$data, FA_sector, sector = 1)
  
  # Secondo settore (11-18)
  sector2 = map_dbl(tract$data, FA_sector, sector = 2)
  
  # Terzo settore (19-28)
  sector3 = map_dbl(tract$data, FA_sector, sector = 3)
  
  # Quarto settore (29-34)
  sector4 = map_dbl(tract$data, FA_sector, sector = 4)
  
  # Quinto settore (35-50)
  sector5 = map_dbl(tract$data, FA_sector, sector = 5)
  
  return (cbind(sector1, sector2, sector3, sector4, sector5))
}

# I settori sono stati scelti in base ai risultati forniti da HICA

label_left_right <- function(data)
{
  label = NULL
  for (i in 1:length(data))
  {
    if ((data[[i]]$x[1]) < 0)
    {
      label=c(label, "left")
    }
    else 
    {
      label=c(label, "right")
    }
  }
  return (label)
}


# save(create_dataset, create_reduced_dataset, get_reduced_dataset, file="functions_create_dataset.RData")
