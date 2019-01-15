library(tidyverse)
library(devtools)
library(readr)
library(boot)
library(gtools)
library(Gmedian)
library(sfsmisc)
library(cluster)
library(fiber)
library(fda)
library(purrr)


source("dataset_utils.R")

setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")
setwd("C:/Users/User/OneDrive - Politecnico di Milano/Project StatApp/RData")

load("cst_list.RData")
# load("cst00_rep_no_outliers.RData")


#########################################################################################
################################ CURVATURA TORSIONE SUE #################################
#########################################################################################
cst01 = cst_list[[1]]
# lambda_opt = 3.04284112980894
MD_info = list(4, c(27,36,40))
AD_info = list(4, c(14,27,34))
RD_info = list(4, c(21,31,37))
FA_info = list(5, c(9,21,29,35))

cst01_left_features = create_dataset_new (cst01$lhs, "left", MD_info = MD_info, RD_info = RD_info, AD_info=AD_info,
                                           FA_info = FA_info, standardized = F) 



cst01_right_features = create_dataset_new (cst01$rhs, "right", MD_info = MD_info, RD_info = RD_info, AD_info=AD_info,
                                            FA_info = FA_info, standardized = F)

setwd("C:/Users/User/OneDrive - Politecnico di Milano/Project StatApp/RData/features_rep_no_std")
save(cst01_left_features, cst01_right_features, file = "cst_01_features.RData")

# Creation of a big list
setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")

load("cst01_features.RData")
load("cst02_features.RData")
load("cst03_features.RData")
load("cst04_features.RData")
load("cst05_features.RData")
load("cst06_features.RData")
load("cst07_features.RData")
load("cst08_features.RData")
load("cst09_features.RData")
load("cst10_features.RData")
load("cst11_features.RData")
load("cst12_features.RData")
load("cst13_features.RData")
load("cst14_features.RData")
load("cst15_features.RData")
load("cst16_features.RData")
load("cst17_features.RData")
load("cst18_features.RData")
load("cst19_features.RData")
load("cst20_features.RData")


features_list = list( cst01_features, 
                      cst02_features, 
                      cst03_features, 
                      cst04_features, 
                      cst05_features, 
                      cst06_features,
                      cst07_features, 
                      cst08_features, 
                      cst09_features,
                      cst10_features,
                      cst11_features,
                      cst12_features,
                      cst13_features,
                      cst14_features,
                      cst15_features,
                      cst16_features,
                      cst17_features,
                      cst18_features,
                      cst19_features,
                      cst20_features)
select_left = function(features_pat) {
  return(features_pat$lhs)
}
select_right = function(features_pat) {
  return(features_pat$rhs)
}

features_left = map_dfr(features_list, select_left)
mean_left = colMeans(features_left[,1:33])
sd_left = apply(features_left[,1:33], 2, sd)
features_right = map_df(features_list, select_right)
mean_right = colMeans(features_right[,1:33])
sd_right = apply(features_right[,1:33], 2, sd)
setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")

save(features_list, mean_left, sd_left, mean_right, sd_right, file = "features_list.RData")
save(features_list, mean_left, sd_left, mean_right, sd_right, file = "features_list.RData")


