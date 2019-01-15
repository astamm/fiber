source("read_utils.R")
source("outliers_utils.R")

###########  Lettura tratto
cst = read_csv("/06007", case = "007", scan = "001") 

###########  Pulizia rispetto a RD e AD
cst$Streamlines = map(cst$Streamlines, remove_point_outliers) 

###########  Divide left and right
cst_divided = divide_cst(cst)

###########  Outlier detection
outliers_left = get_outliers_distance(cst_divided$lhs)
cst_divided$lhs$Streamlines = cst_divided$lhs$Streamlines[-outliers_left]

outliers_right = get_outliers_distance(cst_divided$rhs)
cst_divided$rhs$Streamlines = cst_divided$rhs$Streamlines[-outliers_right]

cst07 = cst_divided

###########  Reparametrization
cst07$lhs = reparametrize_with_interpolation(cst07$lhs)
cst07$rhs = reparametrize_with_interpolation(cst07$rhs)

save(cst07, file = "cst07_rep_no_outliers.RData")



setwd("/Users/ILARIASARTORI/Politecnico di Milano/Luca Torriani - Project StatApp/RData")

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

cst_list = list(cst01, 
                cst02, 
                cst03, 
                cst04, 
                cst05, 
                cst06,
                cst07, 
                cst08, 
                cst09,
                cst10,
                cst11,
                cst12,
                cst13,
                cst14,
                cst15,
                cst16,
                cst17,
                cst18,
                cst19,
                cst20)

save(cst_list, file = "cst_list.RData")