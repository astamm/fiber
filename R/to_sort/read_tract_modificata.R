source("utility_functions.R")

read_csv = function(folder, name, case, scan, side = NA_character_){
  setwd(paste(getwd(),folder, sep=""))
  setwd(paste(getwd(),"/DIFF_30DIR_CST", sep=""))

  cst <- read_tract("DIFF_30DIR_CST_0.csv", name = name, case = case, scan = scan, side = side)
  return(cst)
}


read_tract <- function(path, name, case, scan, side = NA_character_) {
  data <- readr::read_csv(path)

  if (ncol(data) < 5L)
    stop("The input CSV file should contain at least 5 columns.")

  required_vars <- c("X", "Y", "Z", "PointId", "StreamlineId")
  if (!all(required_vars %in% names(data)))
    stop("The CSV file should contain at least the variables X, Y, Z, PointId and StreamelineId.")

  if (nrow(data) == 0)
    return(tract(
      name = name,
      case = case,
      scan = scan,
      side = side,
      data = list()
    ))


  data <- data %>%
    dplyr::arrange(StreamlineId, PointId) %>%
    dplyr::mutate(
      Tensors = purrr::pmap(list(`Tensors#0`, `Tensors#1`, `Tensors#2`,
                                 `Tensors#3`,`Tensors#4`, `Tensors#5`), c) %>%
        purrr::map(as_tensor)
      )

  data = data %>%
    dplyr::group_by(StreamlineId) %>%
    dplyr::do(streamlines = streamline(
      x = .$X, y = .$Y, z = .$Z, diffusion = .$Tensors
    )) %>%
    dplyr::ungroup()


  tract = tract(name = name, case = case, scan = scan, side = side, data = data$streamlines)

  AD_list = purrr::map(tract$data, get_axial_diffusivity_streamline)
  RD_list = purrr::map(tract$data, get_radial_diffusivity_streamline)
  for(i in 1:length(tract$data)){
    tract$data[[i]] = tract$data[[i]] %>%
                        dplyr::mutate(
                          ad = AD_list[[i]],
                          rd = RD_list[[i]]
                        )
  }

  lower_rd <- 0.00001
  upper_rd <- 0.001
  lower_ad <- 0.001
  upper_ad <- 0.0024

  data_filtered = purrr::map(tract$data, ~dplyr::filter(., lower_rd < rd, rd < upper_rd, lower_ad < ad, ad < upper_ad) %>% as_streamline())

  for(i in 1:length(tract$data)){
    tract$data[[i]] = dplyr::select(data_filtered[[i]], -ad, -rd)
  }


  return (tract)

}


# # ESEMPIO DI UTILIZZO
# setwd("C:/Users/User/Politecnico di Milano/Aymeric Stamm - fdatractography")
# cst01_prova = read_csv("/06001", name = "patient1", case = "1", scan = "01")



