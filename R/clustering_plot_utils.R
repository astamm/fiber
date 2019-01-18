# Filtra tenendo solo le streamline del cluster j
get_stream_cl_j = function(features, j){
  return (filter(features, clust == j))
  # return ((features[which(features$clust == j),]))
}

# Filtra tenendo solo le streamline del paziente i
get_stream_pat_i = function(features, i){
  return (filter(features, patient==i))
  # return ((features[which(features$patient == i),]))
}

# Restituisce una lista con 2 elementi:
# - lhs: che contiene tutte le features sinistre
# - rhs: che contiene tutte le features destre
split_left_right = function(features){
  left_right = list (lhs = filter (features, side =="left"),
                     rhs = filter (features, side == "right"))
  # left_right = list(lhs = features[features$side == "left",],
  #                   rhs = features[features$side == "right",])
  return(left_right)
}

# Restituisce una lista di n_cluster elementi. Ogni elemento è a sua volta una lista di 2 elementi:
# -lhs: contiene tutte le features sinistre appartenenti al cluster j
# -rhs: contiene tutte le features destre appartenenti al cluster j
get_streamlines_per_clusters_list = function(features){
  n_cluster = max(as.numeric(as.character(features[[1]]$clust)))
  clusters_list = list()
  # for(j in 1:n_cluster){
  #   assign(paste("cluster_",j,sep=""),map_df(features,get_stream_cl_j,j))
  #   clusters_list [[j]] =  eval(parse(text=paste("cluster_",j,sep ="")))
  # }
  for(j in 1:n_cluster){
    clusters_list [[j]] = purrr::map_df(features,get_stream_cl_j,j)
  }
  clusters_list = purrr::map(clusters_list, split_left_right)
  return (clusters_list)
}

# Restituisce una lista di 20 elementi, che sono a loro volta liste di 2 elementi.
# - lhs: le features sinistre del paziente i
# - rhs: le features destre del paziente j
get_streamlines_per_patients_list = function(features_cl_j, n_patients){
  patients_list = list()

  features_cl_j_sxdx = as.data.frame(rbind(features_cl_j$lhs, features_cl_j$rhs))

  for(i in 1:n_patients){
    # patients_list [[i]] = map_df(features_cl_j,get_stream_pat_i,i)
    patients_list [[i]] = get_stream_pat_i(features_cl_j_sxdx,i)
  }
  patients_list = map(patients_list, split_left_right)

  return (patients_list)
}

# Restituisce una lista di n_cluster*n_patients elementi,  in ordine di cluster e, all'interno del cluster, in ordine di paziente
# Ogni elemento è a sua volta una lista di 2 elementi:
# -lhs: contiene tutte le features sinistre appartenenti al cluster j
# -rhs: contiene tutte le features destre appartenenti al cluster j
get_streamlines_per_clusters_patient_list = function (features, n_patients) {
  clusters_list = get_streamlines_per_clusters_list(features)
  cluster_patients_list = map(clusters_list, get_streamlines_per_patients_list, n_patients = n_patients)

  cluster_patients_list_unpack = list()
  for (j in 1: length(cluster_patients_list)) {
    cluster_patients_list_unpack = c(cluster_patients_list_unpack, cluster_patients_list[[j]])
  }
  return (cluster_patients_list_unpack) # In ordine di cluster e, all'interno del cluster, in ordine di paziente
  # return (cluster_patients_list)
}

# Restituisce una lista di 2 elementi:
# -lhs: streamline media (fittizia) sinistra
# -rhs: streamline media (fittizia) destra
get_clusterMean_streamline = function(left_right_list){
  left_streamline = colMeans(left_right_list$lhs[,1:33])
  right_streamline = colMeans(left_right_list$rhs[,1:33])
  return (list(lhs = left_streamline, rhs =right_streamline))
}


# Restituisce un dataframe con 3 variabili: indexes, side e clust.
# - indexes: indici delle features che stanno nel cluster j (relativamente al cluster e alla side sx-dx)
get_indexes_streamline_in_cluster_j = function (features, j) {
  features_left = filter(features, side == "left")
  idx_left = which( as.numeric(as.character(features_left$clust)) == j)
  left = as.data.frame(idx_left)
  colnames(left) = "indexes"
  left = mutate (left, side =rep("left",length(idx_left)), clust =rep(j, length(idx_left)))

  features_right = filter(features, side == "right")
  idx_right = which( as.numeric(as.character(features_right$clust)) == j)
  right = as.data.frame(idx_right)
  colnames(right) = "indexes"
  right = mutate (right, side =rep("right",length(idx_right)), clust =rep(j, length(idx_right)))

  df = as.data.frame(rbind(left, right))
  return (df)
}


# Restituisce una lista di n_cluster elementi. Ciasuno è una lista di 20 liste di 2 elementi
# -lhs: get_indexes_streamline_in_cluster_j per paziente i, cluster j, side left
# -rhs: get_indexes_streamline_in_cluster_j per paziente i, cluster j, side right
get_indexes_streamlines_per_clusters_list = function(features_list_sxdx){
  n_cluster = max(features_list_sxdx[[1]]$clust)
  clusters_list = list()
  # for(j in 1:n_cluster){
  #   assign(paste("cluster_",j,sep=""),map(features_list_sxdx,get_indexes_streamline_in_cluster_j,j))
  #   clusters_list [[j]] =  eval(parse(text=paste("cluster_",j,sep ="")))
  # }
  for(j in 1:n_cluster){
    clusters_list [[j]] =  map(features_list_sxdx, get_indexes_streamline_in_cluster_j, j=j)
  }
  clusters_list = map(clusters_list, map, split_left_right)
  return (clusters_list)
}


# Data una streamline fittizia (centroide) restituisce una lista con
# - idx: l'indice della streamline reale più vicina (indice rispetto al paziente, la side e il cluster)
# - patient: il paziente a cui questa streamline appartiene
# "streamlines_cluster_j_side" è
#  * in Single_rep_cluster: l'insieme delle streamline di tutti i pazienti che stanno nel cluster j
#  * in Multiple_rep_cluster: l'insieme delle streamline del paziente i che stanno nel cluster j
find_index_nearest_streamline_cluster_j = function(centroid, streamlines_cluster_j_side, single_multiple_rep_cluster) {
  distances = as.numeric(apply(streamlines_cluster_j_side, 1, distance_from_centroid, centroid=centroid))
  idx = which.min(distances)
  patient = streamlines_cluster_j_side[idx,35]

  if (single_multiple_rep_cluster == "Single_rep_cluster") {
    real_idx = which(which(streamlines_cluster_j_side[,35]==patient) == idx)
  }
  else if (single_multiple_rep_cluster== "Multiple_rep_cluster") {
    real_idx= idx
  }
  else {
    stop ("Wrong indication of the requested function find_index_nearest_streamline_cluster_j: \nSingle_rep_cluster or Multiple_rep_cluster")
  }

  index_patient = list(idx = real_idx, patient = patient)
  return (index_patient)
}

# * Single_rep_cluster:
#   Riceve un centroide (fittizio) e le streamline di tutti i pazienti appartenenti a quel cluster
#   Restituisce una lista di 2 elementi:
#    - lhs: una lista con $idx: l'indice della streamline reale più vicina (indice rispetto al paziente, la side e il cluster)
#                         $patient: il paziente a cui questa streamline appartiene
#    - rhs: una lista con $idx: l'indice della streamline reale più vicina (indice rispetto al paziente, la side e il cluster)
#                         $patient: il paziente a cui questa streamline appartiene
# * Multiple_rep_cluster:
#   Riceve un centroide (fittizio) e le streamline di un certo paziente i appartenenti a quel cluster
#   Restituisce una lista di 2 elementi:
#    - lhs: una lista con $idx: l'indice della streamline reale più vicina (indice rispetto al paziente, la side e il cluster)
#                         $patient: il paziente a cui questa streamline appartiene
#    - rhs: una lista con $idx: l'indice della streamline reale più vicina (indice rispetto al paziente, la side e il cluster)
#                         $patient: il paziente a cui questa streamline appartiene
find_indexes_per_cluster_j = function(centroid_cluster_j, streamlines_cluster_j, single_multiple_rep_cluster){
  left = find_index_nearest_streamline_cluster_j(centroid_cluster_j$lhs, streamlines_cluster_j$lhs, single_multiple_rep_cluster)
  right = find_index_nearest_streamline_cluster_j(centroid_cluster_j$rhs, streamlines_cluster_j$rhs, single_multiple_rep_cluster)

  return (list(lhs=left, rhs=right))
}

# Riceve
# - un indice (indice rispetto al paziente, cluster e side)
# - un vettore di indici (indici rispetto al paziente e alla side delle streamline che appartengono al cluster j)
# Restituisce
# - un indice (indice rispett al paziente e alla sie)
find_indexes_base = function (idx, index_patient_i_side){
  result = index_patient_i_side[idx,1]
  return (as.numeric(as.character((result))))
}

# * Single_rep_cluster
# Restituisce una lista di 2 elementi:
# - lhs: lista di 2 elementi - idx: indice (rispetto al paziente e alla side) della streamline rappresentante il cluster j
#                            - patient: paziente a cui appartiene la streamline rappresentante
# - rhs: lista di 2 elementi - idx: indice (rispetto al paziente e alla side) della streamline rappresentante il cluster j
#                            - patient: paziente a cui appartiene la streamline rappresentante
# * Multiple_rep_cluster
# Restituisce una lista di 2 elementi:
# - lhs: lista di 2 elementi - idx: vettore di indici (rispetto al paziente e alla side) della streamline rappresentante il cluster j per il paziente i
#                            - patient: vettore dei pazienti i, a cui si riferisce l'indice sopra
# - rhs: lista di 2 elementi - idx: vettore di indici (rispetto al paziente e alla side) della streamline rappresentante il cluster j per il paziente i
#                            - patient: vettore dei pazienti i, a cui si riferisce l'indice sopra
find_indexes = function(indexes_patients_cluster_j, indexes_cluster_j, single_multiple_rep_cluster){
  patient_left = indexes_patients_cluster_j$lhs$patient
  patient_right = indexes_patients_cluster_j$rhs$patient

  # patient_left e patient_right sono due interi da 1 a 20
  # indexes_patients_cluster_j$lhs$idx e indexes_patients_cluster_j$rhs$idx sono due interi
  if (single_multiple_rep_cluster == "Single_rep_cluster") {
    idx_left = find_indexes_base(indexes_patients_cluster_j$lhs$idx, indexes_cluster_j[[patient_left]]$lhs)
    idx_right = find_indexes_base(indexes_patients_cluster_j$rhs$idx, indexes_cluster_j[[patient_right]]$rhs)

    # left_list = list(idx = idx_left, patient = patient_left, clust = as.numeric(as.character((indexes_cluster_j[[patient_left]]$lhs[1,3]))))
    # right_list = list(idx = idx_right, patient = patient_right, clust = as.numeric(as.character(indexes_cluster_j[[patient_right]]$rhs[1,3])))
    left_list = list(idx = idx_left, patient = patient_left)
    right_list = list(idx = idx_right, patient = patient_right)
  }

  # patient_left e patient_right sono i vettori degli interi da 1 a 20
  # indexes_patients_cluster_j$lhs$idx e indexes_patients_cluster_j$rhs$idx sono due vettori di 20 interi
  else if (single_multiple_rep_cluster == "Multiple_rep_cluster") {
    idx_left = NULL
    idx_right = NULL
    for (i in 1: length(patient_left)) {
      idx_left = c (idx_left, find_indexes_base(indexes_patients_cluster_j$lhs$idx[i], indexes_cluster_j[[patient_left[i]]]$lhs))
      idx_right = c (idx_right, find_indexes_base(indexes_patients_cluster_j$rhs$idx[i], indexes_cluster_j[[patient_right[i]]]$rhs))
    }
    left_list = list(idx = idx_left, patient = patient_left)
    right_list = list(idx = idx_right, patient = patient_right)

  }

  else {
    stop ("Wrong indication of the requested function find_indexes: \nSingle_rep_cluster or Multiple_rep_cluster")
  }

  return (list (lhs = left_list, rhs = right_list))
}


# Data una lista di n_patients*n_clusters elementi la raggruppa in una lista (lunga n_clusters)
# di liste (lunghe n_patients)
collect_in_clusters = function (indexes_pat_clust, n_patients, n_clusters) {
  indexes_clust = list()
  for (j in 1:n_clusters) {
    left_indexes_cluster_j = NULL
    right_indexes_cluster_j = NULL
    left_patients_cluster_j = NULL
    right_patients_cluster_j = NULL

    for (i in 1: n_patients) {
      left_indexes_cluster_j = c(left_indexes_cluster_j, indexes_pat_clust [[(j-1)*n_patients + i]]$lhs$idx)
      right_indexes_cluster_j = c(right_indexes_cluster_j, indexes_pat_clust [[(j-1)*n_patients + i]]$rhs$idx)

      left_patients_cluster_j = c(left_patients_cluster_j, indexes_pat_clust [[(j-1)*n_patients + i]]$lhs$patient)
      right_patients_cluster_j = c(right_patients_cluster_j, indexes_pat_clust [[(j-1)*n_patients + i]]$rhs$patient)
    }

    indexes_clust[[j]] = list(lhs=list(idx = left_indexes_cluster_j, patient =left_patients_cluster_j),
                              rhs=list(idx = right_indexes_cluster_j, patient =right_patients_cluster_j))

  }
  return (indexes_clust)
}


extract_data_sx = function (data) {
  return(data$lhs)
}
extract_data_dx = function (data) {
  return(data$rhs)
}

# evaluate_variance = function(features_clusterj, mean_vec, sd_vec) {
#   std_data = sweep(features_clusterj[,1:33], 2, mean_vec)  # Subtract the mean
#   std_data = sweep(std_data, 2, sd_vec, FUN = "/")   # Divide by the standard deviation
#   return(cov(std_data))
# }


extract_var_sx = function (data) {
  return(data$var_sx)
}

extract_var_dx = function (data) {
  return(data$var_dx)
}

find_max = function(data) {
  return (map_dbl(data, max))
}

find_min = function(data) {
  return (map_dbl(data, min))
}


# Date features_reduced e features_list_sxdx (solo pazienti sani) restituisce
# una lista di 9 (cluster), dove ogni elemento è una lista di 2 componenti:
# - idx: Indice della streamline reale (rispetto a un paziente e alla side) più vicina (distanza L2 delle features)
#        alla streamline fittizia media di quel cluster (media mettendo insieme tutti i pazienti)
# - patient: Paziente a cui appartiene la streamline rappresentante
get_cluster_true_mean_indexes_healthy = function(features_reduced,features_list_sxdx, mean_left, sd_left, mean_right, sd_right){

  representative_clusters_list = get_streamlines_per_clusters_list(features_reduced)
  # Lista di liste: Ogni elemento della list ? un cluster che a sua volta ? una lista
  # contente lhs e rhs (ciascuno dei due ? un dataframe con le streamlines rappresentanti)

  clusterMean_streamlines = map(representative_clusters_list, get_clusterMean_streamline)
  # Lista di 9 cluster, ciascun elemento ? diviso in lhs e rhs e ciascuno contiene la streamline media rappresentante

  streamlines_per_cluster = get_streamlines_per_clusters_list(features_list_sxdx)
  # lista di 9 elementi (ogni cluster) divisi a loroa volta in lhs e rhs contententi tutte le streamline del cluster j unendo tutti
  # i pazienti

  # streamlines_per_cluster_left = map(streamlines_per_cluster, extract_data_sx)
  # streamlines_per_cluster_right = map(streamlines_per_cluster, extract_data_dx)
  # var_sx = map(streamlines_per_cluster_left, evaluate_variance, mean_vec=mean_left, sd_vec=sd_left)
  # var_dx = map(streamlines_per_cluster_right, evaluate_variance, mean_vec=mean_right, sd_vec=sd_right)

  indexes_per_cluster = get_indexes_streamlines_per_clusters_list(features_list_sxdx)
  # Lista di 9 (clusters) di 20 liste (pazienti) ciascuno divisa in lhs e rhs contenenti: indexes  side cluster
  # (rispetto al singolo paziente)

  # max(indexes_per_cluster[[1]][[1]]$rhs$indexes)   # 4460
  # # Attenzione: la colonna degli indici viene letta come un factor
  # healthy_tract_lengths_matrix[1,2]   # 3126

  indexes_patients = map2(clusterMean_streamlines,streamlines_per_cluster,
                          find_indexes_per_cluster_j, single_multiple_rep_cluster = "Single_rep_cluster")

  final_indexes = map2(indexes_patients,indexes_per_cluster,find_indexes, single_multiple_rep_cluster = "Single_rep_cluster")
  return (list(final_indexes=final_indexes))
}


# Date features_reduced e features_list_sxdx (solo pazienti sani) restituisce
# una lista di 9 (cluster), dove ogni elemento e' una lista di due elementi (lhs, rhs),
# ognuno dei quali e' una lista di 2 componenti:
# - idx: Vettore di indici delle streamline reale (rispetto a un paziente e alla side) più vicina (distanza L2 delle features)
#        alla streamline fittizia media di quel cluster per quel paziente
# - patient: Vettore dei pazienti
get_cluster_patient_true_mean_indexes_healthy = function(features_reduced,features_list_sxdx){
  n_patients = length(features_list_sxdx)

  clusterpatientMean_streamlines = get_streamlines_per_clusters_patient_list (features_reduced, n_patients)
  # Lista di 9*20 elementi, corrispondenti a ciascuna combinazione cluster-paziente, ognuno diviso in lhs e rhs.
  # Ciascuno contiene la streamline media rappresentante

  streamlines_per_patient_cluster = get_streamlines_per_clusters_patient_list (features_list_sxdx, n_patients)
  # Lista di 9*20 elementi, corrispondenti a ciascuna combinazione cluster-paziente, ognuno diviso in lhs e rhs.
  # Ciascuno contiene tutte le streamline di quel paziente in quel cluster

  indexes_per_cluster = get_indexes_streamlines_per_clusters_list(features_list_sxdx)
  # Lista di 9 (clusters) di 20 liste (pazienti) ciascuno divisa in lhs e rhs contenenti: indexes  side cluster
  # (rispetto al singolo paziente)
  n_clusters = length(indexes_per_cluster)

  indexes_pat_clust = map2(clusterpatientMean_streamlines,streamlines_per_patient_cluster,
                          find_indexes_per_cluster_j,  single_multiple_rep_cluster = "Multiple_rep_cluster")
  # Lista di 9*20 elementi, corrispondenti a ciascuna combinazione cluster-paziente, ognuno diviso in lhs e rhs.
  # Ciascuno contiene l'indice della streamline media vera (indice nel paziente i cluster j)
  indexes_patients = collect_in_clusters(indexes_pat_clust, n_patients, n_clusters)
  # Lista di 9 liste di 20 elementi, ognuno diviso in lhs e rhs.
  # Ciascuno contiene l'indice della streamline media vera (indice nel paziente i cluster j)

  final_indexes = map2(indexes_patients,indexes_per_cluster,find_indexes, single_multiple_rep_cluster = "Multiple_rep_cluster")
  # Lista di 9 liste di 20 elementi, corrispondenti a ciascuna combinazione cluster-paziente, ognuno diviso in lhs e rhs.
  # Ciascuno contiene l'indice della streamline media vera (indice nel paziente i)
  return (final_indexes)
}

plot_from_indexes = function (cst_list_healthy, indexes_plot) {
  open3d()
  primo=1
  n_cluster = length(indexes_plot)
  for (j in 1:n_cluster){
    for (s in 1:2){   # Left e right
      patient = indexes_plot[[j]][[s]]$patient
      index = indexes_plot[[j]][[s]]$idx
      for (i in 1: length(patient)) {
        x= cst_list_healthy[[patient[i]]][[s]]$Streamlines[[index[i]]]$x
        y= cst_list_healthy[[patient[i]]][[s]]$Streamlines[[index[i]]]$y
        z= cst_list_healthy[[patient[i]]][[s]]$Streamlines[[index[i]]]$z

        X = cbind(x,y,z)

        lines3d(X, asp=1, size=0.5, col=rainbow(n_cluster)[j],axes=FALSE)

        if (primo) {
          axes3d()
          title3d(xlab='x', ylab='y', zlab='z')
          primo=0
        }
      }
    }
  }
}


# NEW : 14/12/18
get_streamline = function(row_df, side, cst_list){
  pat = row_df[2]
  idx = row_df[1]
  if(side == "left")  stream = cst_list[[pat]]$lhs$Streamlines[[idx]]
  else if (side == "right") stream = cst_list[[pat]]$rhs$Streamlines[[idx]]
  else stop("Error side")
  return (stream)
}

# NEW 17/01/18
rescale = function (features_patient, mean_tot, sd_tot){
  features_rescaled = features_patient
  features_rescaled[,1:33] = sweep(features_rescaled[,1:33], 2, mean_tot  )
  # features_rescaled$rhs[,1:33] = sweep(features_patient$rhs[,1:33], 2, mean_tot  )
  features_rescaled[,1:33] = sweep(features_rescaled[,1:33], 2, sd_tot, FUN = "/")
  # features_rescaled$rhs[,1:33] = sweep(features_rescaled$rhs[,1:33], 2, sd_vec, FUN = "/")
  return(features_rescaled)
}

compute_variance = function (features_patient, side){
  features = NULL
  n = NULL
  fac = NULL
  tmp_var = NULL
  var = NULL
  if(side == "left"){
    features = features_patient$lhs
    fac = factor(features$clust)
    n = length(levels(fac))
    tmp_var = split(as.data.frame((features[,1:33])), features$clust)
    var = purrr::map(tmp_var, cov)
  } else if(side == "right"){
    features = features_patient$rhs
    fac = factor(features$clust)
    n = length(levels(fac))
    tmp_var = split(as.data.frame((features[,1:33])), features$clust)
    var = purrr::map(tmp_var, cov)
  }  else stop("Error side")
  return(var)
}

cov_cluster = function(cluster_feature, id_cols){
  return(cov(cluster_feature[,id_cols]))
}
