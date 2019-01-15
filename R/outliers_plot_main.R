source("dataset_utils.R")
source("outliers_utils.R")
source("read_utils.R")

###############################################################################
#######                  GRAFICI                      #########################
###############################################################################

### 1) Read tract and remove points whith RD and AD out of range

# Lettura tratto
setwd("/Users/ILARIASARTORI/Desktop/")
cst = read_csv("/06001", case = "001", scan = "001") 

# Pulizia rispetto a RD e AD
cst$Streamlines = map(cst$Streamlines, remove_point_outliers) 

# Divide left and right
cst = divide_cst(cst)

### 2) Plots
tract = cst$lhs

####################################################################################
##################################### DISTANCE #####################################
####################################################################################
outliers = get_outliers_distance(tract, separate_sp_bar=T)
outliers.sp_xy = outliers$outliers.sp_xy 
outliers.sp_yz = outliers$outliers.sp_yz 
outliers.bar_yz = outliers$outliers.bar_yz 
outliers.bar_xy = outliers$outliers.bar_xy

library(rgl)
primo=1
for (i in 1:length(tract$Streamlines)) {
  if (primo) {
    if (is.element(i, unique(c(outliers.bar_yz, outliers.bar_xy)))) {
      plot(tract$Streamlines[[i]], col='green')
    }
    else if (is.element(i, unique(c(outliers.sp_xy, outliers.sp_yz)))){
      plot(tract$Streamlines[[i]], col='blue')
    }
    else
      plot(tract$Streamlines[[i]], col='gray')
    axes3d()
    title3d(xlab='x', ylab='y', zlab='z')
    primo=0
  }
  if (is.element(i, unique(c(outliers.bar_yz, outliers.bar_xy)))) {
    plot(tract$Streamlines[[i]], col='green', new_window=FALSE)
  }
  else if (is.element(i, unique(c(outliers.sp_xy, outliers.sp_yz)))){
    plot(tract$Streamlines[[i]], col='blue', new_window=FALSE)
  }
  else
    plot(tract$Streamlines[[i]], col='gray', new_window=FALSE)
}


# Per vedere bene la differenza, fare anche un altro grafico in cui si evidenziano le 
# streamline solo classificate outliers per il baricentro (basta invertire la classificazione
# per colore nell'if)
primo=1
for (i in 1:length(tract$Streamlines)) {
  if (primo) {
    if (is.element(i, unique(c(outliers.sp_xy, outliers.sp_yz)))){
      plot(tract$Streamlines[[i]], col='blue')
    }
    else if (is.element(i, unique(c(outliers.bar_yz, outliers.bar_xy)))) {
      plot(tract$Streamlines[[i]], col='green')
    }
    else
      plot(tract$Streamlines[[i]], col='gray')
    axes3d()
    title3d(xlab='x', ylab='y', zlab='z')
    primo=0
  }
  if (is.element(i, unique(c(outliers.sp_xy, outliers.sp_yz)))){
    plot(tract$Streamlines[[i]], col='blue', new_window=FALSE)
  }
  else if (is.element(i, unique(c(outliers.bar_yz, outliers.bar_xy)))) {
    plot(tract$Streamlines[[i]], col='green', new_window=FALSE)
  }
  else
    plot(tract$Streamlines[[i]], col='gray', new_window=FALSE)
}

# Mi sembrerebbe che le streamline outliers dal punto di vista della spatial median siano
# quelle che hanno la parte centrale abbastanza allineata, ma piu' ci avviciniamo alla corteccia
# piu sono strane: il "ramo" e' troppo basso o troppo alto rispetto alle altre streamline, 
# o magari e' molto piu' lungo
# Le streamline outliers solo per il baricentro invece forse sono quelle che hanno anche il tronco un po'
# spostato

####################################################################################
##################################### DEPTH #####################################
####################################################################################
library(fields)
outliers = get_outliers_depth(tract, separate_sp_bar=T)
outliers_depth_Barycenter = outliers$outliers_depth_Barycenter
outliers_depth_Median = outliers$outliers_depth_Median

primo=1
for (i in 1:length(tract$Streamlines)) {
  if (primo) {
    if (is.element(i, outliers_depth_Barycenter)){
      plot(tract$Streamlines[[i]], col='green')
    }
    else if (is.element(i, outliers_depth_Median)) {
      plot(tract$Streamlines[[i]], col='blue')
    }
    else
      plot(tract$Streamlines[[i]], col='gray')
    axes3d()
    title3d(xlab='x', ylab='y', zlab='z')
    primo=0
  }
  if (is.element(i, outliers_depth_Barycenter)){
    plot(tract$Streamlines[[i]], col='green', new_window=FALSE)
  }
  else if (is.element(i, outliers_depth_Median)) {
    plot(tract$Streamlines[[i]], col='blue', new_window=FALSE)
  }
  else
    plot(tract$Streamlines[[i]], col='gray', new_window=FALSE)
}

# Per vedere bene la differenza, fare anche un altro grafico in cui si evidenziano le 
# streamline solo classificate outliers per il baricentro (basta invertire la classificazione
# per colore nell'if)
primo=1
for (i in 1:length(tract$Streamlines)) {
  if (primo) {
    if (is.element(i, outliers_depth_Median)) {
      plot(tract$Streamlines[[i]], col='blue')
    }
    else if (is.element(i, outliers_depth_Barycenter)){
      plot(tract$Streamlines[[i]], col='green')
    }
    else
      plot(tract$Streamlines[[i]], col='gray')
    axes3d()
    title3d(xlab='x', ylab='y', zlab='z')
    primo=0
  }
  if (is.element(i, outliers_depth_Median)) {
    plot(tract$Streamlines[[i]], col='blue', new_window=FALSE)
  }
  else if (is.element(i, outliers_depth_Barycenter)){
    plot(tract$Streamlines[[i]], col='green', new_window=FALSE)
  }
  else
    plot(tract$Streamlines[[i]], col='gray', new_window=FALSE)
}

