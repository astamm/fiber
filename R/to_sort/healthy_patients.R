library(tidyverse)
library(multidplyr)
library(roahd)
library(fdatractography)
source("R/to_sort/all_patients.R")

# Load database -----------------------------------------------------------

P <- 50L
cl <- create_cluster(cores = 4L) %>%
  cluster_library("fdatractography") %>%
  cluster_library("modelr") %>%
  cluster_library("dplyr") %>%
  cluster_library("purrr") %>%
  cluster_copy(P) %>%
  cluster_copy(approx_tensors)

project_folder <- "~/OneDrive - Politecnico di Milano/fdatractography/"
df <- tibble()
for (i in 1:20) {
  ID <- paste0("06", formatC(i, width = 3, flag = "0"))
  print(ID)
  path <- paste0(project_folder, ID, "/DIFF_30DIR_CST/DIFF_30DIR_CST_0.csv")
  tmp <- load_tract(path, ID, cl)
  df <- bind_rows(df, tmp)
}

save(df, file = "~/OneDrive - Politecnico di Milano/fdatractography/fdakma/healthy_patients_step1.RData")

# Filter out length outliers ----------------------------------------------

out <- 0
df1 <- df

while (length(out) > 0) {
  df1 <- df1 %>%
    mutate(len = map_dbl(data, get_curvilinear_length))

  out <- boxplot.stats(df1$len)$out
  print(length(out))

  df1 <- df1 %>%
    filter(!(len %in% out)) %>%
    select(-len)
}

save(df1, file = "~/OneDrive - Politecnico di Milano/fdatractography/fdakma/healthy_patients_step2.RData")

N <- nrow(df1)

# Define the common grid of curvilinear abscissa --------------------------

s_min <- df1$data %>%
  map(pull, s) %>%
  map_dbl(min) %>%
  max()
s_max <- df1$data %>%
  map(pull, s) %>%
  map_dbl(max) %>%
  min()
grid <- seq(s_min, s_max, length.out = P)

# Reparametrize all streamlines using common grid -------------------------

cl <- create_cluster(cores = 4L) %>%
  cluster_library("fdatractography") %>%
  cluster_library("purrr") %>%
  cluster_library("dplyr") %>%
  cluster_copy(approx_tensors) %>%
  cluster_copy(grid)

df2 <- df1 %>%
  partition(cluster = cl) %>%
  mutate(data = map(data, ~ streamline(
    s = grid,
    x = approx(.$s, .$x, xout = s)$y,
    y = approx(.$s, .$y, xout = s)$y,
    z = approx(.$s, .$z, xout = s)$y,
    t = approx_tensors(.$s, .$t, xout = s)
  ))) %>%
  collect() %>%
  ungroup() %>%
  arrange(StreamlineId) %>%
  select(-PARTITION_ID)

save(df2, file = "~/OneDrive - Politecnico di Milano/fdatractography/fdakma/healthy_patient_step3.RData")

# Set up the arrays required by the fdakma package ------------------------

data <- array(dim = c(N, P, 3L))
for (i in 1:N) {
  data[i, , 1] <- df2$data[[i]]$x
  data[i, , 2] <- df2$data[[i]]$y
  data[i, , 3] <- df2$data[[i]]$z
}

# Mirror x coordinates to bring left CST onto right CST -------------------

lbls <- df2$side
pid <- df2$case
data[, , 1] <- abs(data[, , 1])

# Save the grid and data arrays to disk -----------------------------------

save(
  lbls, pid, grid, data,
  file = "~/OneDrive - Politecnico di Milano/fdatractography/fdakma/healthy_patients_step4.RData"
)

# Plots
data1 <- data[pid == "06001", ,]
lbls1 <- lbls[pid == "06001"]
fda::matplot(grid, t(data1[, , 2]), type = "l", col = c("red", "blue")[(lbls1 == "R") + 1])
fda::matplot(grid, t(data1[, , 1]), type = "l", col = c("red", "blue")[(lbls1 == "R") + 1])

set.seed(1234)

# Functional boxplot ------------------------------------------------------

out <- 0
fd <- fData(grid, data[, , 2])
lbls_red <- lbls
pid_red <- pid
data_red <- data
outliers_amplitude <- NULL

while (length(out) > 0) {
  fb <- fbplot(
    fd,
    adjust = list(
      N_trials = 10,
      trial_size = 5 * fd$N,
      VERBOSE = TRUE
    ),
    display = FALSE
  )
  out <- fb$ID_outliers
  outliers_amplitude[[length(outliers_amplitude) + 1]] <- out
  keep <- setdiff(1:fd$N, out)
  fd <- fd[keep]
  data_red <- data_red[keep, , ]
  lbls_red <- lbls_red[keep]
  pid_red <- pid_red[keep]
}

fda::matplot(grid, t(data_red[, , 2]), type = "l", col = c("red", "blue")[(lbls_red == "R") + 1])
fda::matplot(grid, t(data_red[, , 1]), type = "l", col = c("red", "blue")[(lbls_red == "R") + 1])

save(
  lbls_red, pid_red, grid, data_red, outliers_amplitude,
  file = "~/OneDrive - Politecnico di Milano/fdatractography/fdakma/healthy_patients_step5.RData"
)

data_mir <- data_red
data_mir[lbls_red == "L", , 1] <- -data_mir[lbls_red == "L", , 1]
fda::matplot(grid, t(data_mir[, , 1]), type = "l", col = c("red", "blue")[(lbls_red == "R") + 1])

# Outliergram -------------------------------------------------------------

out <- 0
fd <- fData(grid, data_red[, , 2])
lbls_final <- lbls_red
pid_final <- pid_red
data_final <- data_red
outliers_shape <- NULL

while (length(out) > 0) {
  og <- outliergram(
    fd,
    adjust = list(
      N_trials = 10,
      trial_size = 5 * fd$N,
      VERBOSE = TRUE
    ),
    display = FALSE
  )
  out <- og$ID_outliers
  print(out)
  outliers_shape[[length(outliers_shape) + 1]] <- out
  keep <- setdiff(1:fd$N, out)
  fd <- fd[keep]
  data_final <- data_final[keep, , ]
  lbls_final <- lbls_final[keep]
  pid_final <- pid_final[keep]
}

data1 <- data_final[pid_final == "06001", ,]
lbls1 <- lbls_final[pid_final == "06001"]
fda::matplot(grid, t(data1[, , 2]), type = "l", col = c("red", "blue")[(lbls1 == "R") + 1])
fda::matplot(grid, t(data1[, , 1]), type = "l", col = c("red", "blue")[(lbls1 == "R") + 1])

save(
  lbls_final, pid_final, grid, data_final, outliers_shape,
  file = "~/OneDrive - Politecnico di Milano/fdatractography/fdakma/healthy_patients_step6.RData"
)

data_mir <- data_final
data_mir[lbls_final == "L", , 1] <- -data_mir[lbls_final == "L", , 1]
fda::matplot(grid, t(data_mir[, , 1]), type = "l", col = c("red", "blue")[(lbls_final == "R") + 1])
