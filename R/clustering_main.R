library(tidyverse)
library(cluster)
library(rgl)
library (fiber)
library(dplyr)

######################################################################################
################################### Load tracts ######################################
######################################################################################

setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")
load("cst_list.RData")

source("~/Desktop/OneDrive - Politecnico di Milano/CST_atlas/Clustering/helper_cluster.R")
source("/Users/ILARIASARTORI/Desktop/Poli/IVanno/CST_atlas/Clustering/helper_cluster.R")

# For each patient, get the number of streamlines left+right sides
num_streamline_patients = map_dbl (cst_list, num_of_streamline_patient)


################################################################################################
######################### Create dataset points on the cortex punti ############################
################################################################################################
names_col = names(cst_list[[1]]$lhs$Streamlines[[1]])
ind_x = which(names_col=="x")
ind_y = which(names_col=="y")
ind_z = which(names_col=="z")
data_first_points = map(cst_list, dataset_first_points, ind_x, ind_y, ind_z)
# Mappa di 20 elementi. Ogni elemento è il dataset dei punti iniziali delle streamline di quel paziente

save(data_first_points, file = "data_first_points.RData")


################################################################################################
######################################## Clustering ############################################
################################################################################################

#### Load dataset points on cortex 
load("data_first_points.RData")
data_first_points_healthy_matrix = map_df (data_first_points, rbind)


#### Cluster sani
n=50
threshold = 20
num_healthy_patients = 20
k_opt = get_k_opt(data_first_points_healthy_matrix, n, threshold, num_healthy_patients) ## 9
clara = get_cluster_clara(data_first_points_healthy_matrix, k_opt)  # prima di clara c'è un set.seed(1994)
cluster_healthy = clara$clustering

# Choiche of the threshold
k_opt_vec=NULL
threshold_vec = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70)
for (threshold in threshold_vec) {
  k_opt_vec = c(k_opt_vec, get_k_opt(data_first_points_healthy_matrix, n, threshold, num_healthy_patients))
}

x11()
plot(threshold_vec, k_opt_vec, pch = 19, xlab = "Threshold", ylab = "k_opt", main = "k_opt - threshold")
lines(threshold_vec, k_opt_vec)

####  In data_first_points aggiungo la colonna con i cluster di appartenenza
data_first_points = map(data_first_points, add_cluster_column, clusters = cluster_healthy, num_streamline_patients=num_streamline_patients)

#### Riproietto il tratto sinistro nel piano delle x negative
data_first_points = map(data_first_points, reproject_x)

################################################################################################
############################### Creo features_patients_9.RData ###############################
################################################################################################
### Loading features
load("features_list.RData")

# # Lista di o elementi.
# # Ogni elemento rappresenta un paziente
# # Ogni elemento (i.e paziente) è rappresentato da una lista di 2 elementi: le features sinistre e le destre
#######################################################################################

# Lista di 20 elementi. Ogni elemento è un dataframe con le features delle varie streamline di quel paziente,
# prima quelle sinistre e poi le destre
features_list_sxdx = map (features_list, merge_left_right_features)

# Aggiungo gli indici del cluster
features_list_sxdx = map (features_list_sxdx, add_cluster_column, clusters = cluster_healthy, num_streamline_patients=num_streamline_patients)  # features_final

save (features_list_sxdx, mean_left, sd_left, mean_right, sd_right, file="features_patients_9.RData")
#######################################################################################


################################################################################################
########################## Creo features_patients_reduced_9.RData ############################
################################################################################################
setwd("~/Desktop/Politecnico di Milano/Luca Torriani - Project StatApp/RData")
load("features_patients_9.RData")

# Calcolo i centroidi
features_reduced_tmp = purrr::map (features_list_sxdx, get_reduced_tot, mean_left=mean_left, sd_left=sd_left, mean_right=mean_right, sd_right=sd_right)

features_reduced = map(features_reduced_tmp, extract_centroid) # Centroids
var_sx_features_reduced = map(features_reduced_tmp, extract_var_sx)
var_dx_features_reduced = map(features_reduced_tmp, extract_var_dx)

save (features_reduced, var_sx_features_reduced, var_dx_features_reduced, file="features_patients_reduced_9.RData")



