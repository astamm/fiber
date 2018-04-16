library(tidyverse)
library(fdatractography)
library(Gmedian)
library(fda)

# ________________________________________________________________________________________________
# Getter diffusivity information

get_axial_diffusivity_streamline <- function(streamline, validate = TRUE)
{
  if (validate) {
    if (!is_streamline(streamline))
      stop("The input dataset is not of class streamline.")
  }
  return (map_dbl(streamline$diffusion, get_axial_diffusivity))
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
# Evaluate curvature and torsion using splines of order m, choosing the optimum lambda according
# to generalized cross validation
get_curvature_torsion_fda3D_vec = function (streamline, m = 6L, Lfdobj = 4L, ubound = 100) {
  s <- streamline$s
  res <- get_functional_representation(streamline, norder = m, Lfdobj = Lfdobj, ubound = ubound)
  deriv1 <- fda::eval.fd(s, fda::as.fd(res$fd), Lfd = 1L)
  deriv2 <- fda::eval.fd(s, fda::as.fd(res$fd), Lfd = 2L)
  deriv3 <- fda::eval.fd(s, fda::as.fd(res$fd), Lfd = 3L)

  # Calcolo vettore curvature
  curv <- (sqrt(  (deriv2[,3]*deriv1[,2] - deriv2[,2]*deriv1[,3])^2 +
                        (deriv2[,1]*deriv1[,3] - deriv2[,3]*deriv1[,1])^2 +
                        (deriv2[,2]*deriv1[,1] - deriv2[,1]*deriv1[,2])^2) )   /
    (deriv1[,1]^2 + deriv1[,2]^2 + deriv1[,3]^2)^(3/2)
  tors <- (  deriv3[,1]*(deriv1[,2]*deriv2[,3] - deriv2[,2]*deriv1[,3]) +
                  deriv3[,2]*(deriv2[,1]*deriv1[,3] - deriv1[,1]*deriv2[,3]) +
                  deriv3[,3]*(deriv1[,1]*deriv2[,2] - deriv2[,1]*deriv1[,2]))   /
    (  (deriv1[,2]*deriv2[,3] - deriv2[,2]*deriv1[,3])^2 +
         (deriv2[,1]*deriv1[,3] - deriv1[,1]*deriv2[,3])^2 +
         (deriv1[,1]*deriv2[,2] - deriv2[,1]*deriv1[,2])^2)

  list(curvature = curv, torsion = tors, lambda = res$lambda)
}

get_curvature_torsion_fda3D = function (streamline, m = 6L, Lfdobj = 4L, ubound = 100) {
  l <- get_curvature_torsion_fda3D_vec(streamline, m, Lfdobj, ubound)
  c(
    max(l$curvature), mean(l$curvature), sd(l$curvature),
    max(l$torsion), mean(l$torsion), sd(l$torsion)
  )
}

# Helper functions
get_functional_representation <- function(streamline, norder = 6L, Lfdobj = 4L, ubound = 100) {
  s <- streamline$s
  basis <- fda::create.bspline.basis(s, norder = norder)
  Y <- array(dim = c(length(s), 1L, 3L))
  Y[, , 1] <- streamline$x
  Y[, , 2] <- streamline$y
  Y[, , 3] <- streamline$z
  opt <- optimise(gcv_lambda, c(1e-5, ubound), basis = basis, Lfdobj = Lfdobj, Y = Y, s = s)
  optLambda <- opt$minimum
  optFDPar <- fda::fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = optLambda)
  res <- fda::smooth.basis(s, Y, optFDPar)
  list(fd = res, lambda = optLambda)
}

gcv_lambda <- function(lambda, basis, Lfdobj, Y, s) {
  functionalPar <- fda::fdPar(fdobj = basis, Lfdobj = Lfdobj, lambda = lambda)
  mean(fda::smooth.basis(s, Y, functionalPar)$gcv)
}

get_curvature_torsion_fda3D_vec(cst10_left$data[[1]])

system.time(
  l <- map(cst10_left$data[1:1000], get_curvature_torsion_fda3D_vec, ubound = 25)
)
median(unlist(transpose(l)$lambda))
mean(unlist(transpose(l)$lambda))
cc <- get_curvature(cst10_left$data[[1]] %>% as_streamline())
plot(l[[1]]$curvature, type = "l")
plot(cc$s, cc$k, type = "l")
