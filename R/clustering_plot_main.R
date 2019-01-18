library(fields)
library(rgl)
library(fiber)

source("clustering_plot_utils.R")
source("clustering_utils.R")

setwd("~/Desktop/Politecnico di Milano/Luca Torriani - Project StatApp/RData")
setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")
setwd("C:/Users/User/OneDrive - Politecnico di Milano/Project StatApp/RData")

load("features_patients_9.RData")
load("features_patients_reduced_9.RData")
load("cst_list.RData")

################################################################################################
######################################### Plot variance ########################################
################################################################################################

# Create plots variance intracluster patient per patient

# # No standardization
# id_pat = 8
# i = 3
# side = "left"
# quartz()
# tmp= dplyr::filter(features_list_sxdx[[id_pat]], (features_list_sxdx[[id_pat]]$clust==i & features_list_sxdx[[id_pat]]$side==side))
# image.plot(cov(tmp[,1:33]),axes=F, main = paste("Patient",id_pat, "Cluster",i, "Side", side))

features_rescaled_unito = purrr::map(features_list_sxdx, rescale, mean_tot, sd_tot)
features_rescaled = purrr::map(features_rescaled_unito, split_left_right)

# Standardization:feature_pat_i_clu_j_scaled = (feature_pat_i_cluj- mean(feature_pat_i))/ sd(feature_pat_i) NON PIU VERO
var_sx_features_reduced = purrr::map(features_rescaled, compute_variance, "left")
var_dx_features_reduced = purrr::map(features_rescaled, compute_variance, "right")
for (side in c("left", "right")) {
  for (id_pat in 1:20) {
    z_max_left = max(map_dbl(var_sx_features_reduced[[id_pat]],max))
    z_min_left =  min(map_dbl(var_sx_features_reduced[[id_pat]],min))
    z_max_right = max(map_dbl(var_dx_features_reduced[[id_pat]],max))
    z_min_right = min(map_dbl(var_dx_features_reduced[[id_pat]],min))
    z_max = max(z_max_left, z_max_right)
    z_min = min(z_min_left, z_min_right)
    z_lim=c(z_min,z_max)

    pdf(paste0('Variance_intracluster-Pat_',id_pat,'-Side_',side,'.pdf'))
    par(mfrow=c(3,3), mai = c(0.1,0.1,0.3,0.1), oma=c(1,1,4,1))
    for (i in 1:9){
      if(side=="left") image.plot(var_sx_features_reduced[[id_pat]][[i]], zlim=z_lim, axes=F, main = paste("Cluster",i))
      else image.plot(var_dx_features_reduced[[id_pat]][[i]], zlim=z_lim, axes=F, main = paste("Cluster",i))
    }
    mtext(paste("Patient",id_pat, ", ", "Side ", side), side = 3, line = 1, outer = TRUE, font=2)
    dev.off()
  }
}


######################################################################################
################################### First plot #######################################
######################################################################################

indexes_plot1 = get_cluster_patient_true_mean_indexes_healthy (features_reduced, features_list_sxdx)
plot_from_indexes (cst_list, indexes_plot1)



######################################################################################
################################## Second plot #######################################
######################################################################################
indexes_plot2 = get_cluster_true_mean_indexes_healthy (features_reduced, features_list_sxdx, mean_left=mean_left, sd_left=sd_left, mean_right=mean_right, sd_right=sd_right)
plot_from_indexes (cst_list, indexes_plot2$final_indexes)


# Create plots variance intracluster

# lista di 9 elementi (ogni cluster) divisi a loroa volta in lhs e rhs contententi tutte le streamline del cluster j unendo tutti
# i pazienti
streamlines_per_cluster = get_streamlines_per_clusters_list(features_rescaled_unito)
streamlines_per_cluster_left = purrr::map(streamlines_per_cluster, extract_data_sx)
streamlines_per_cluster_right = purrr::map(streamlines_per_cluster, extract_data_dx)

var_sx = purrr::map(streamlines_per_cluster_left, cov_cluster,1:33)
var_dx = purrr::map(streamlines_per_cluster_right, cov_cluster,1:33)

z_max_left = max(map_dbl(var_sx,max))
z_min_left = min(map_dbl(var_sx,min))
z_max_right = max(map_dbl(var_dx,max))
z_min_right = min(map_dbl(var_dx,min))
z_max = max(z_max_left, z_max_right)
z_min = min(z_min_left, z_min_right)
z_lim=c(z_min,z_max)

for (side in c("left", "right")) {
  pdf(paste0('Variance_intracluster','-Side_',side,'.pdf'))
  par(mfrow=c(3,3), mai = c(0.1,0.1,0.3,0.1), oma=c(1,1,4,1))
  for (i in 1:9){
    if(side=="left") image.plot(var_sx[[i]], zlim=z_lim, axes=F, main = paste("Cluster",i))
    else image.plot(var_dx[[i]], zlim=z_lim, axes=F, main = paste("Cluster",i))
  }
  mtext(paste("Side ", side), side = 3, line = 1, outer = TRUE, font=2)
  dev.off()
}

# PROBLEMA : CURVATURA HA DEI PUNTI STRANI, VEDIAMOLI
library(dplyr)
library(plyr)

df_big <- ldply (features_list_sxdx, data.frame)

for(i in 19:19){
  pdf(file = paste0(colnames(df)[i],".pdf", sep=""))
  par(mfrow = c(2,1))

  plot(df[,i], main = colnames(df)[i] )
  boxplot(df[,i])
  dev.off()
}




######################################################################################
################################### Third plot #######################################
######################################################################################
indexes_plot3 = get_cluster_true_mean_indexes_healthy (features_reduced, features_list_sxdx, mean_left=mean_left, sd_left=sd_left, mean_right=mean_right, sd_right=sd_right)

indexes = indexes_plot3$final_indexes
labels_df = do.call(rbind, map(indexes,unlist))
df_left = labels_df[,c(1,2)]
df_right = labels_df[,c(3,4)]
stream_left = apply(df_left, 1, get_streamline, side = "left", cst_list)
stream_right = apply(df_right, 1, get_streamline, side = "right", cst_list)
streams = do.call(c,list(stream_left,stream_right))

colours =rainbow(9)
plot(streams[[1]], plot_microstructure = TRUE, col = colours[1], scale = 6)
j= 2
for (i in 2:length(streams)){
  if(j==(length(streams)/2+1)) j = 1
  plot(streams[[i]], plot_microstructure = TRUE, new_window = F, col=colours[j], scale = 6)
  j =j+1
}



