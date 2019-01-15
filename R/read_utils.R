library(fiber)
library(purrr)

# source("utility_functions.R")
source("helper_outliers.R")

read_csv = function(folder, case, scan){
  setwd(paste(getwd(),folder, sep=""))
  setwd(paste(getwd(),"/DIFF_30DIR_CST", sep=""))
  
  cst <- fiber::read_tract("DIFF_30DIR_CST_0.csv", name = "CST", case = case, scan = scan, biomarkers = 'SingleTensor')
  return(cst)
}

remove_point_outliers = function(stream, l_rd = 0.00001, u_rd = 0.001, l_ad = 0.001, u_ad = 0.0024){
  microstruct = map_dfr(stream$t, get_microstructure)
  idx = which((microstruct$AD<l_ad) | (microstruct$AD>u_ad) | (microstruct$RD< l_rd) | (microstruct$RD>u_rd)) 
  if(!is_empty(idx)){
    stream = stream[-idx,]
  }
  return(streamline_cast(stream))
}

# Non posso usare il constructor di tract perch? non ? aggiornato
divide_cst = function(cst){
  # cst_left$data <- cst$data[map_lgl(cst$data, ~ (.$x[1] < 0))]
  # cst_right$data <- cst$data[map_lgl(cst$data, ~ (.$x[1] > 0))]
  
  cst_left <- cst_right <- cst
  cst_left$Streamlines <- cst$Streamlines[map_lgl(cst$Streamlines, ~ (.$x[1] < 0))]
  cst_right$Streamlines <- cst$Streamlines[map_lgl(cst$Streamlines, ~ (.$x[1] > 0))]
  # cst_left = tract(name = cst$TractName[1], case = cst$PatientId[1], scan = cst$ScanId[1], side = "left", data = data_left) 
  # cst_right = tract(name = cst$name[1], case = cst$case[1], scan = cst$scan[1], side = "right", data = data_right) 
  cst_left$Streamlines <- map(cst_left$Streamlines, streamline_cast)
  cst_right$Streamlines <- map(cst_right$Streamlines, streamline_cast)
  
  return (list(lhs = cst_left, rhs = cst_right))
}

streamline_cast = function(stream) {
  return(streamline(s=stream$s, x=stream$x, y=stream$y, z=stream$z, t=stream$t))
}

library(dplyr)

reparametrize_streamline = function(streamline, n_points) {
  s = modelr::seq_range(streamline$s, n = n_points)
  x = approx(streamline$s, streamline$x, xout = s)$y
  y = approx(streamline$s, streamline$y, xout = s)$y
  z = approx(streamline$s, streamline$z, xout = s)$y
  t = approx_tensors(streamline$s, streamline$t, xout = s)$y
  
  biomarkers_df = map_df(t, get_microstructure)
  return(streamline(s,x,y,z,t, ad=biomarkers_df$AD, rd=biomarkers_df$RD, md=biomarkers_df$MD, fa=biomarkers_df$FA))
}

reparametrize_with_interpolation = function(tract, n=50){
  if(!is_tract(tract)){
    stop("A tract is needed")
  }
  n_points= n
  tract$Streamlines = map(tract$Streamlines, reparametrize_streamline, n_points=n_points)
  
  return(tract)
}