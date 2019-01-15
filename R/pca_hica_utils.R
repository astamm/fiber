# MD
get_patient_MD = function(tract) {
  patient_MD_left = sapply(tract$lhs$Streamlines, get_streamline_MD)  # Restituisce una matrice (cbind delle MD)
  patient_MD_right = sapply(tract$rhs$Streamlines, get_streamline_MD) 
  return(cbind(patient_MD_left, patient_MD_right))
}

get_streamline_MD = function(streamline) {
  return(streamline$md)
}

# AD
get_patient_AD = function(tract) {
  patient_AD_left = sapply(tract$lhs$Streamlines, get_streamline_AD)  # Restituisce una matrice (cbind delle MD)
  patient_AD_right = sapply(tract$rhs$Streamlines, get_streamline_AD) 
  return(cbind(patient_AD_left, patient_AD_right))
}

get_streamline_AD = function(streamline) {
  return(streamline$ad)
}

# RD
get_patient_RD = function(tract) {
  patient_RD_left = sapply(tract$lhs$Streamlines, get_streamline_RD)  # Restituisce una matrice (cbind delle MD)
  patient_RD_right = sapply(tract$rhs$Streamlines, get_streamline_RD) 
  return(cbind(patient_RD_left, patient_RD_right))
}

get_streamline_RD = function(streamline) {
  return(streamline$rd)
}

# FA
get_patient_FA = function(tract) {
  patient_FA_left = sapply(tract$lhs$Streamlines, get_streamline_FA)  # Restituisce una matrice (cbind delle MD)
  patient_FA_right = sapply(tract$rhs$Streamlines, get_streamline_FA) 
  return(cbind(patient_FA_left, patient_FA_right))
}

get_streamline_FA = function(streamline) {
  return(streamline$fa)
}

# PCA functions

plot_pca = function (features, pc.features) {
  x11()
  layout(matrix(c(2,3,1,3),2,byrow=T))
  # variance of PC
  barplot(pc.features$sdev^2, las=2, main='Principal Components',ylim=c(0,max(pc.features$sdev^2)), ylab='Variances')  #  ylim=c(0,0.3), ## MD
  # original variances 
  abline(h=1, col='blue')
  barplot(sapply(as.data.frame(features),sd)^2, las=2, main='Original Variables',ylim=c(0,max(pc.features$sdev^2)), ylab='Variances')  #  ylim=c(0,0.1) ## MD
  # proportion of explained variance
  plot(cumsum(pc.features$sdev^2)/sum(pc.features$sde^2), type='b', axes=F, xlab='number of components', 
       ylab='contribution to the total variability', ylim=c(0,1))
  abline(h=0.6, lty=2, col='blue')
  box()
  axis(2,at=0:10/10,labels=0:10/10)
  axis(1,at=1:ncol(features),labels=1:ncol(features),las=1)
  
  ## loadings
  load.features    <- pc.features$loadings
  
  # Bar plot dei loadings PC; 
  x11()
  par( mfrow = c(4,1), mar=c(0,3,0,3))
  for(i in 1:4) barplot(load.features[,i], ylim = c(-0.5, 0.5))
}



# HICA functions
get_HICA_loadings <- function (dataset, K, title, coda = FALSE)
{
  
  if (coda)
  {
    dataset = dataset[,26:50]
    lev = 24
  }
  
  else 
  {
    lev = 49
  }
  
  
  library(fastHICA)
  
  basis <- basis_hica(as.matrix(dataset))
  energy <- energy_hica(basis, maxcomp = K)
  hica <- extract_hica(energy, comp = K, level=lev)
  
  
  hica.load <- hica$C
  rownames(hica.load) = as.character(c(1:50))
  
  x11()
  par(mar = c(1,4,0,2), oma=c(0,0,3,0), mfrow = c(K,1))
  for(i in 1:K)
  {
    barplot(hica.load[,i], ylim = c(-1, 1), col=rainbow(lev+1))
    abline(h=0)
  }
  title(paste("HICA for ", title), outer = T)  
  
}



