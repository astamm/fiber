setwd("/Users/ILARIASARTORI/Google drive/Progetto_StatApp/File per paper")
setwd("C:/Users/User/Google Drive/Progetto_StatApp/File per paper")
source ("utility_functions.R")
source("helper_outliers.R")

##################################################################################
#######                  GRAFICI                      #########################
##################################################################################
setwd("/Users/ILARIASARTORI/Google drive/Progetto_StatApp/File per paper/RData_tratti_riparametrizzati")


load("cst01_reparametrized.RData")
tract = cst01_left

####################################################################
## far girare il corpo di GET_OUTLIERS_DISTANCE

library(rgl)
open3d()
primo=1
for (i in 1:length(tract$data)) {
  x= tract$data[[i]]$x
  y= tract$data[[i]]$y
  z= tract$data[[i]]$z
  n=nrow(tract$data[[i]])

  X=cbind(x,y,z)

  color = NULL
  for (j in 1:n) {
    if (is.element(i, unique(c(outliers.bar_yz, outliers.bar_xy))))
      color[j]= 'green'
    else if (is.element(i, unique(c(outliers.sp_xy, outliers.sp_yz))))
      color[j]= 'blue'
    else
      color[j]= 'gray'
  }

  lines3d(X[which(color=='gray'),], asp=1, size=0.5, col="gray")
  lines3d(X[which(color=='blue'),], asp=1, size=0.5, col="blue")
  lines3d(X[which(color=='green'),], asp=1, size=0.5, col="green")

  if (primo) {
    axes3d()
    title3d(xlab='x', ylab='y', zlab='z')
    primo=0
  }
}

# Per vedere bene la differenza, fare anche un altro grafico in cui si evidenziano le
# streamline solo classificate outliers per il baricentro (basta invertire la classificazione
# per colore nell'if)
open3d()
primo=1
for (i in 1:length(cst01_left$data)) {
  x= cst01_left$data[[i]]$x
  y= cst01_left$data[[i]]$y
  z= cst01_left$data[[i]]$z
  n=nrow(cst01_left$data[[i]])

  X=cbind(x,y,z)

  color = NULL
  for (j in 1:n) {
    if (is.element(i, unique(c(outliers.sp_xy, outliers.sp_yz))))
      color[j]= 'blue'
    else if (is.element(i, unique(c(outliers.bar_yz, outliers.bar_xy))))
      color[j]= 'green'
    else
      color[j]= 'gray'
  }

  lines3d(X[which(color=='gray'),], asp=1, size=0.5, col="gray")
  lines3d(X[which(color=='blue'),], asp=1, size=0.5, col="blue")
  lines3d(X[which(color=='green'),], asp=1, size=0.5, col="green")

  if (primo) {
    axes3d()
    title3d(xlab='x', ylab='y', zlab='z')
    primo=0
  }
}

# Mi sembrerebbe che le streamline outliers dal punto di vista della spatial median siano
# quelle che hanno la parte centrale abbastanza allineata, ma piu' ci avviciniamo alla corteccia
# piu sono strane: il "ramo" e' troppo basso o troppo alto rispetto alle altre streamline,
# o magari e' molto piu' lungo
# Le streamline outliers solo per il baricentro invece forse sono quelle che hanno anche il tronco un po'
# spostato

####################################################################
## far girare il corpo di GET_OUTLIERS_DEPTH
library(rgl)
open3d()
primo=1
for (i in 1:length(tract$data)) {
  x= tract$data[[i]]$x
  y= tract$data[[i]]$y
  z= tract$data[[i]]$z
  n=nrow(tract$data[[i]])

  X=cbind(x,y,z)

  color = NULL
  for (j in 1:n) {
    if (is.element(i, outliers_depth_Barycenter))
      color[j]= 'green'
    else if (is.element(i, outliers_depth_Median))
      color[j]= 'blue'
    else
      color[j]= 'gray'
  }

  lines3d(X[which(color=='gray'),], asp=1, size=0.5, col="gray")
  lines3d(X[which(color=='blue'),], asp=1, size=0.5, col="blue")
  lines3d(X[which(color=='green'),], asp=1, size=0.5, col="green")

  if (primo) {
    axes3d()
    title3d(xlab='x', ylab='y', zlab='z')
    primo=0
  }
}


##################################################################################
#######                  SALVATAGGIO                      #########################
##################################################################################
# setwd("C:/Users/Vale/Google Drive/Progetto_StatApp/File per paper/Rdata_tratti_riparametrizzati")
setwd("/Users/ILARIASARTORI/Google drive/Progetto_StatApp/File per paper/RData_tratti_riparametrizzati")
load("cst04_reparametrized.RData")

######## DISTANCE - BASED
# Left
system.time(
  outliers_sx <- get_outliers_distance(cst10_left)
)
cst04_left$data = cst04_left$data[-outliers_sx]
# length(cst04_left$data)
# Right
outliers_dx = get_outliers_distance(cst04_right)
cst04_right$data = cst04_right$data[-outliers_dx]

######## DEPTH - BASED
# Left
outliers_sx = get_outliers_depth(cst04_left)
cst04_left$data = cst04_left$data[-outliers_sx]
# length(cst04_left$data)
# Right
outliers_dx = get_outliers_depth(cst04_right)
cst04_right$data = cst04_right$data[-outliers_dx]


# setwd("C:/Users/Vale/Google Drive/Progetto_StatApp/File per paper/Rdata_tratti_riparametrizzati_no_outliers")
setwd("/Users/ILARIASARTORI/Google drive/Progetto_StatApp/File per paper/RData_tratti_riparametrizzati_no_outliers")
save(cst04_left, cst04_right, file = "cst04_rep_no_outliers.RData")

## Verifica
# load("cst04_rep_no_outliers.RData")
# length(cst04_left$data)
#
# rm(cst01_left)
# rm(cst01_right)
# rm(cst02_left)
# rm(cst02_right)
# rm(cst03_left)
# rm(cst03_right)
# rm(cst04_left)
# rm(cst04_right)

