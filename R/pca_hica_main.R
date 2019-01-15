library(tidyverse)
library(devtools)
library(readr)
library(fdatractography)
library(boot)
library(gtools)
library(Gmedian)
library(sfsmisc)
library(cluster)

source("pca_hica_utils.R")
setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")
load("cst_list.RData")

#### MD
tract_tot_MD_list = map(cst_list, get_patient_MD)
tract_tot_MD = do.call(cbind.data.frame, tract_tot_MD_list)

tract_tot_MD=t(tract_tot_MD)
save(tract_tot_MD, file="tract_tot_MD.RData")

#### AD
tract_tot_AD_list = map(cst_list, get_patient_AD)
tract_tot_AD = do.call(cbind.data.frame, tract_tot_AD_list)

tract_tot_AD=t(tract_tot_AD)
save(tract_tot_AD, file="tract_tot_AD.RData")

#### RD
tract_tot_RD_list = map(cst_list, get_patient_RD)
tract_tot_RD = do.call(cbind.data.frame, tract_tot_RD_list)

tract_tot_RD=t(tract_tot_RD)
save(tract_tot_RD, file="tract_tot_RD.RData")


setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")
load("tract_tot_FA.RData")
load("tract_tot_MD.RData")
load("tract_tot_RD.RData")
load("tract_tot_AD.RData")


# PCA analysis

### FA
pc.featuresFA <- princomp(tract_tot_FA, scores=T)
summary(pc.featuresFA)
plot_pca(tract_tot_FA, pc.featuresFA)

### MD
pc.featuresMD <- princomp(tract_tot_MD, scores=T)
summary(pc.featuresMD)
plot_pca(tract_tot_MD, pc.featuresMD)

### RD
pc.featuresRD <- princomp(tract_tot_RD, scores=T)
summary(pc.featuresRD)
plot_pca(tract_tot_RD, pc.featuresRD)

### AD
pc.featuresAD <- princomp(tract_tot_AD, scores=T)
summary(pc.featuresAD)
plot_pca(tract_tot_AD, pc.featuresAD)


# HICA anlysis
##### MD
# Plot loadings
get_HICA_loadings(tract_tot_MD, 4, "Mean Diffusivity tracts")

##### AD
# Plot loadings
get_HICA_loadings(tract_tot_AD, 4, "Axial Diffusivity tracts")


##### RD
# Plot loadings
get_HICA_loadings(tract_tot_RD, 4, "Radial Diffusivity tracts")


##### FA
# Plot loadings
get_HICA_loadings(tract_tot_FA, 4, "Fractional Anisotropy tracts")

