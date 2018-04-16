library(tidyverse)
library(fdatractography)
library(Gmedian)
library(fda)

# ________________________________________________________________________________________________
# Getter diffusivity information

get_axial_diffusivity <- function(tensor)
{
  eigenvalues = eigen(tensor, symmetric = TRUE)$values
  lambda_max = eigenvalues[1]
  return (lambda_max)
}

get_axial_diffusivity_streamline <- function(streamline, validate = TRUE)
{
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }
  return (map_dbl(streamline$diffusion, get_axial_diffusivity))
}


get_radial_diffusivity <- function(tensor)
{
  eigenvalues = eigen(tensor, symmetric = TRUE)$values
  return ((eigenvalues[2]+eigenvalues[3])/2)
}

get_radial_diffusivity_streamline <- function(streamline, validate = TRUE)
{
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }

  return (map_dbl(streamline$diffusion, get_radial_diffusivity))
}

get_mean_diffusivity_streamline <- function(streamline, validate = TRUE)
{
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }

  return (map_dbl(streamline$diffusion, get_mean_diffusivity))
}

get_fractional_anisotropy_streamline <- function(streamline, validate = TRUE)
{
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }

  return (map_dbl(streamline$diffusion, get_fractional_anisotropy))
}





# ________________________________________________________________________________________________
# Position information


# Barycenter:

get_barycenter_streamline <- function(streamline) {
  # Given a streamline, it returns the 3D coordinates of the barycenter of the curve
  bar=c(mean(streamline$x), mean(streamline$y), mean(streamline$z))
  return (bar)
}

get_barycenter <- function (data) {
  # Creates the feature Barycenter of a give tract, whose streamlines are collected in data

  bar=NULL
  for (i in 1:length(data)) {
    bar=rbind(bar, get_barycenter_streamline(data[[i]]))
  }
  return (bar)
}



# Spatial Median:

get_sp_median_streamline = function(streamline) {
  # Given a streamline, it returns the 3D coordinates of the spatial median of the curve
  spatial_med = Gmedian(streamline)
  return(spatial_med)
}

get_spatial_median = function (data){
  # Creates the feature Spatial Median of a give tract, whose streamlines are collected in data
  sp_median = NULL
  for (i in 1:length(data)){
    sp_median = rbind(sp_median,get_sp_median_streamline(data[[i]][,2:4]))
  }
  sp_median = data.frame(x_SpMed = sp_median[,1],y_SpMed = sp_median[,2],z_SpMed = sp_median[,3] )
  return (sp_median)
}


# Distance spatial_median - barycenter

distance_centers <- function(data, barycenter, spatial_median){
  # Returns the distance between the barycenters and the spatial medians of the streamlines
  # of a given tract
  dist = NULL
  for (i in 1:length(data))
  {
    deltaX = (barycenter[i,1]-spatial_median$x_SpMed[i])^2
    deltaY = (barycenter[i,2]-spatial_median$y_SpMed[i])^2
    deltaZ = (barycenter[i,3]-spatial_median$z_SpMed[i])^2
    dist = c(dist, sqrt(deltaX+deltaY+deltaZ))
  }
  return (dist)
}


#### Curvature and torsion
# FIXED LAMBDA
# Evaluate curvature and torsion using splines of order m
get_curvature_torsion_fda3D_lambda_fixed = function (streamline, m=8, Lfdobj=4, lambda) {
  cat("o")
  s <- streamline$s

  res <- get_functional_representation_lambda_fixed(streamline, norder = m, Lfdobj = Lfdobj, lambda)
  Obs1_left = eval.fd(s, as.fd(res$fd), Lfd=1)
  Obs2_left = eval.fd(s, as.fd(res$fd), Lfd=2)
  Obs3_left = eval.fd(s, as.fd(res$fd), Lfd=3)

  # Calcolo vettore curvature
  curvatura = (sqrt(  (Obs2_left[,3]*Obs1_left[,2] - Obs2_left[,2]*Obs1_left[,3])^2 +
                        (Obs2_left[,1]*Obs1_left[,3] - Obs2_left[,3]*Obs1_left[,1])^2 +
                        (Obs2_left[,2]*Obs1_left[,1] - Obs2_left[,1]*Obs1_left[,2])^2) )   /
    (Obs1_left[,1]^2 + Obs1_left[,2]^2 + Obs1_left[,3]^2)^(3/2)
  torsione = (  Obs3_left[,1]*(Obs1_left[,2]*Obs2_left[,3] - Obs2_left[,2]*Obs1_left[,3]) +
                  Obs3_left[,2]*(Obs2_left[,1]*Obs1_left[,3] - Obs1_left[,1]*Obs2_left[,3]) +
                  Obs3_left[,3]*(Obs1_left[,1]*Obs2_left[,2] - Obs2_left[,1]*Obs1_left[,2]))   /
    (  (Obs1_left[,2]*Obs2_left[,3] - Obs2_left[,2]*Obs1_left[,3])^2 +
         (Obs2_left[,1]*Obs1_left[,3] - Obs1_left[,1]*Obs2_left[,3])^2 +
         (Obs1_left[,1]*Obs2_left[,2] - Obs2_left[,1]*Obs1_left[,2])^2 )
  return (c(max(curvatura), mean(curvatura), sd(curvatura), max(torsione), mean(torsione), sd(torsione), res$lambda))
}

get_functional_representation_lambda_fixed <- function(streamline,norder=8, Lfdobj=4, lambda) {
  s <- streamline$s   # Uso come breaks i punti stessi
  basis <- create.bspline.basis(s, norder = norder)
  Y <- array(dim = c(length(s), 1, 3))
  Y[, , 1] <- streamline$x
  Y[, , 2] <- streamline$y
  Y[, , 3] <- streamline$z

  optFDPar <- fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = lambda)
  res <- smooth.basis(s, Y, optFDPar)

  list(fd = res, lambda = lambda)
}


# VARIABLE LAMBDA
# Choose the optimum lambda according to generalized cross validation
fda3D_evaluate_lambda = function (streamline, m=8, Lfdobj=4, lambda0=10) {
  cat("o")
  s <- streamline$s

  res <- get_functional_representation_evaluate_lambda(streamline, norder = m, Lfdobj = Lfdobj, lambda0 = lambda0)
  Obs1_left = eval.fd(s, as.fd(res$fd), Lfd=1)
  Obs2_left = eval.fd(s, as.fd(res$fd), Lfd=2)
  Obs3_left = eval.fd(s, as.fd(res$fd), Lfd=3)

  return (res$lambda)
}

get_functional_representation_evaluate_lambda <- function(streamline,norder=8, Lfdobj=4, lambda0=10) {
  s <- streamline$s   # Uso come breaks i punti stessi
  basis <- create.bspline.basis(s, norder = norder)
  Y <- array(dim = c(length(s), 1, 3))
  Y[, , 1] <- streamline$x
  Y[, , 2] <- streamline$y
  Y[, , 3] <- streamline$z

  # lambda = optim(lambda0, gcv_lambda, basis=basis, s=s, Y=Y, Lfdobj=Lfdobj, method="Brent",lower = 10^(-6), upper = 50)
  lambda = optimize (gcv_lambda, c(10^(-6),50), basis=basis, s=s, Y=Y, Lfdobj=Lfdobj, maximum = FALSE)
  # lambda = optim(lambda0, gcv_lambda, basis=basis, Lfdobj=Lfdobj)

  optFDPar <- fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = lambda$minimum)
  res <- smooth.basis(s, Y, optFDPar)

  list(fd = res, lambda = lambda$minimum)
}

gcv_lambda = function (lambda, basis, Lfdobj, s, Y) {
  functionalPar <- fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = lambda)
  return (mean(smooth.basis(s, Y, functionalPar)$gcv))
}
