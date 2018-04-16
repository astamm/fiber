library(tidyverse)
library(devtools)
library(readr)
library(fdatractography)
library(boot)
library(gtools)
library(Gmedian)
library(sfsmisc)
library(cluster)
library(fda)

setwd("/Users/ILARIASARTORI/Google drive/Progetto_StatApp/File per paper")
# setwd("C:/Users/Vale/Google Drive/Progetto_StatApp/File per paper")
setwd("C:/Users/User/Google Drive/Progetto_StatApp/File per paper")

source("utility_functions.R")
source("create_dataset.R")
source("tensor_variables.R")
# setwd("C:/Users/Vale/Google Drive/Progetto_StatApp/File per paper/Rdata_tratti_riparametrizzati_no_outliers")
setwd("/Users/ILARIASARTORI/Google drive/Progetto_StatApp/File per paper/RData_tratti_riparametrizzati_no_outliers")
setwd("C:/Users/User/Google Drive/Progetto_StatApp/File per paper/RData_tratti_riparametrizzati_no_outliers")


load("cst00_rep_no_outliers.RData")
load("cst01_rep_no_outliers.RData")
load("cst02_rep_no_outliers.RData")
load("cst03_rep_no_outliers.RData")
load("cst04_rep_no_outliers.RData")
load("cst05_rep_no_outliers.RData")
load("cst06_rep_no_outliers.RData")
load("cst07_rep_no_outliers.RData")
load("cst08_rep_no_outliers.RData")
load("cst09_rep_no_outliers.RData")
load("cst10_rep_no_outliers.RData")
load("cst11_rep_no_outliers.RData")
load("cst12_rep_no_outliers.RData")
load("cst13_rep_no_outliers.RData")
load("cst14_rep_no_outliers.RData")
load("cst15_rep_no_outliers.RData")
load("cst16_rep_no_outliers.RData")
load("cst17_rep_no_outliers.RData")
load("cst18_rep_no_outliers.RData")
load("cst19_rep_no_outliers.RData")
load("cst20_rep_no_outliers.RData")

cst00__left_features = create_dataset (cst00_left) 
cst00__right_features = create_dataset (cst00_right)

# Prova quello left! (se hai voglia prova a dargli un'occhiata)
cst01__left_features = create_dataset (cst01_left) 
cst01__right_features = create_dataset (cst01_right)

cst02__left_features = create_dataset (cst02_left) 
cst02__right_features = create_dataset (cst02_right)

cst03__left_features = create_dataset (cst03_left) 
cst03__right_features = create_dataset (cst03_right)

cst04__left_features = create_dataset (cst04_left) 
cst04__right_features = create_dataset (cst04_right)

cst05__left_features = create_dataset (cst05_left) 
cst05__right_features = create_dataset (cst05_right)

cst06__left_features = create_dataset (cst06_left) 
cst06__right_features = create_dataset (cst06_right)

cst07__left_features = create_dataset (cst07_left) 
cst07__right_features = create_dataset (cst07_right)

cst08__left_features = create_dataset (cst08_left) 
cst08__right_features = create_dataset (cst08_right)

cst09__left_features = create_dataset (cst09_left) 
cst09__right_features = create_dataset (cst09_right)

cst10__left_features = create_dataset (cst10_left) 
cst10__right_features = create_dataset (cst10_right)

cst11__left_features = create_dataset (cst11_left) 
cst11__right_features = create_dataset (cst11_right)

cst12__left_features = create_dataset (cst12_left) 
cst12__right_features = create_dataset (cst12_right)

cst13__left_features = create_dataset (cst13_left) 
cst13__right_features = create_dataset (cst13_right)

cst14__left_features = create_dataset (cst14_left) 
cst14__right_features = create_dataset (cst14_right)

cst15__left_features = create_dataset (cst15_left) 
cst15__right_features = create_dataset (cst15_right)

cst16__left_features = create_dataset (cst16_left) 
cst16__right_features = create_dataset (cst16_right)

cst17__left_features = create_dataset (cst17_left) 
cst17__right_features = create_dataset (cst17_right)

cst18__left_features = create_dataset (cst18_left) 
cst18__right_features = create_dataset (cst18_right)

cst19__left_features = create_dataset (cst19_left) 
cst19__right_features = create_dataset (cst19_right)

cst20__left_features = create_dataset (cst20_left) 
cst20__right_features = create_dataset (cst20_right)
