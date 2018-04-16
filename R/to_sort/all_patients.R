load_tract <- function(path, ID, cl) {
  readr::read_csv(path) %>%
    dplyr::select(-`Fiber weights`) %>%
    dplyr::mutate(Tensors = purrr::pmap(list(`Tensors#0`, `Tensors#1`, `Tensors#2`,
                                             `Tensors#3`, `Tensors#4`, `Tensors#5`), c) %>%
                    purrr::map(as_tensor)) %>%
    dplyr::select(-dplyr::starts_with("Tensors#")) %>%
    dplyr::arrange(StreamlineId, PointId) %>%
    multidplyr::partition(StreamlineId, cluster = cl) %>%
    dplyr::do(data = streamline(x = .$X, y = .$Y, z = .$Z, t = .$Tensors) %>%
                dplyr::do(streamline(
                  s = seq_range(.$s, n = P),
                  x = approx(.$s, .$x, xout = s)$y,
                  y = approx(.$s, .$y, xout = s)$y,
                  z = approx(.$s, .$z, xout = s)$y,
                  t = approx_tensors(.$s, .$t, xout = s)
                ))) %>%
    dplyr::collect() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(StreamlineId) %>%
    dplyr::mutate(
      case = ID,
      side = c("L", "R")[purrr::map_lgl(data, ~ .$x[1] > 0) + 1]
    )
}
