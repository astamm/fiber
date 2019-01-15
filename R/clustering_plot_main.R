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

var_sx = indexes_plot2$var_sx
var_dx = indexes_plot2$var_dx

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



