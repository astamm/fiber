# It receives as parameter either the left or the right hand-side of the tract and the desired number
# of points. It returns the reparametrized tract interpolating the tensors (with linear invariant interpolation)
# and computing, after the tensor interpolation, the four diffusivity features (AD, RD, MD, FA)

source("read_tract_modificata.R")
reparametrize_with_interpolation = function(tract, n=50){
  if(!is_tract(tract)){
    stop("A tract is needed")
  }
  
  tract = get_invariants(tract)
  n_points= n
  
  for (i in 1:length(tract$data)) {
    s = modelr::seq_range(tract$data[[i]]$s, n = n_points)
    x = approx(tract$data[[i]]$s, tract$data[[i]]$x, xout = s)$y
    y = approx(tract$data[[i]]$s, tract$data[[i]]$y, xout = s)$y
    z = approx(tract$data[[i]]$s, tract$data[[i]]$z, xout = s)$y
    K1= approx(tract$data[[i]]$s, tract$data[[i]]$K1, xout = s)$y
    R2= approx(tract$data[[i]]$s, tract$data[[i]]$R2, xout = s)$y
    R3= approx(tract$data[[i]]$s, tract$data[[i]]$R3, xout = s)$y
    
    tract$data[[i]] = data.frame(s,x,y,z,K1,R2,R3)
  }

  tract = get_biomarkers(tract, n_points)
  return(tract)
  
}

## HELPER FUNCTIONS:

# Given the tract ite returns the 3 invariants
get_invariants <- function (tract) {
  for (i in 1: length(tract$data)) {
    K1 = map_dbl(tract$data[[i]]$diffusion, evaluate_K1)
    R2 = map_dbl(tract$data[[i]]$diffusion, evaluate_R2)
    R3 = map_dbl(tract$data[[i]]$diffusion, evaluate_R3)
    
    tract$data[[i]] = mutate(tract$data[[i]], K1, R2, R3)
  }
  return (tract)
}

evaluate_K1 <- function (diff_tensor) {
  K1 = psych::tr(diff_tensor)
  return (K1)
}

evaluate_R2 <- function (diff_tensor) {
  Dtilde = diff_tensor - psych::tr(diff_tensor)*diag(3)/3
  
  R2 = sqrt((3/2))*norm(Dtilde, "F")/norm(diff_tensor, "F")

  
  return (R2)
}

evaluate_R3 <- function (diff_tensor) {
  Dtilde = diff_tensor - psych::tr(diff_tensor)*diag(3)/3
  R3 = 3*sqrt(6)*det(Dtilde/norm(Dtilde, "F"))
  return (R3)
}

# It computes the four diffusivity features
get_biomarkers <- function (tract, n_points) {
  for (i in 1: length(tract$data)) {
    eigenvalues = NULL
    for (j in 1:n_points)
      eigenvalues = rbind(eigenvalues, evaluate_eigenvalues(tract$data[[i]]$K1[j],tract$data[[i]]$R2[j],tract$data[[i]]$R3[j]))
    
    tract$data[[i]] = tract$data[[i]]%>%
                        dplyr::mutate( 
                             lambda1 = eigenvalues[,1], 
                             lambda2 = eigenvalues[,2], 
                             lambda3 = eigenvalues[,3]) 
    
    tract$data[[i]] = tract$data[[i]][,-c(5,6,7)]
    tract$data[[i]] = add_biomarkers(tract$data[[i]])

  }
  
  return (tract)
}

evaluate_eigenvalues <- function (K1, R2, R3) {
  Coeff1 = 1/3*K1
  Coeff2 = 2*K1*R2/(3*sqrt(3-2*(R2)^2))
  lambda1 = Coeff1 + Coeff2*cos(acos(R3)/3)
  lambda2 = Coeff1 + Coeff2*cos((acos(R3)-2*pi)/3)
  lambda3 = Coeff1 + Coeff2*cos((acos(R3)+2*pi)/3)
  return (c(lambda1, lambda2, lambda3))
}

add_biomarkers <- function(data) {
    m = (data$lambda1+ data$lambda2 + data$lambda3)/3
    data = data %>%
      dplyr::mutate(
        AD = lambda1,
        RD = (lambda2 + lambda3)/3,
        MD = m,
        FA = sqrt(1.5 * ((lambda1-m)^2+(lambda2-m)^2+(lambda3-m)^2) /((lambda1)^2+(lambda2)^2+(lambda3)^2))
      )
    data = data[,-c(5,6,7)]
    return(data)
}




## ESEMPIO DI UTILIZZO
#  setwd("C:/Users/User/Politecnico di Milano/Aymeric Stamm - fdatractography")
# # setwd("C:/Users/Vale/Politecnico di Milano/Aymeric Stamm - fdatractography")
# # setwd("/Users/ILARIASARTORI/Politecnico di Milano/Aymeric Stamm - fdatractography")
# cst05 = read_csv("/06005", name = "patient05", case = "05", scan = "01") 
# cst_left <- cst_right <- cst05
# cst_left$data <- cst05$data[map_lgl(cst05$data, ~ (.$x[1] < 0))]
# cst_right$data <- cst05$data[map_lgl(cst05$data, ~ (.$x[1] > 0))]
# # 
# cst_left = reparametrize_with_interpolation(cst_left)
# cst_right = reparametrize_with_interpolation(cst_right)




