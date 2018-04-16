tensors <- purrr::map(df$data, "t")
for (i in 1:length(outliers_amplitude)) {
  print(i)
  out <- outliers_amplitude[[i]]
  n <- length(tensors)
  keep <- setdiff(1:n, out)
  tensors <- tensors[keep]
}

save(
  tensors,
  file = "~/OneDrive - Politecnico di Milano/fdatractography/fdakma/healthy_patients_tensors_mda.RData"
)
