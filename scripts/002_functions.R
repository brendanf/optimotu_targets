# faster than length(intersect(x, y)), but assumes x = unique(x) and y = unique(y)
# this could be even faster since we know they are sorted, but there is no fast R
# function I am aware of.
intersect_length <- function(x, y) sum(!is.na(match(x, y)))

# compiled c++ version gives ~4x speedup
Rcpp::sourceCpp("src/intersect.cpp")

inner_fmeasure <- function(cj, kpartition, nk) {
  nc <- length(cj)
  nc * max(purrr::map_int(kpartition, intersect_length_, cj)/(nc + nk))
}

# Calculate F-measure for delimitation by clustering
f_measure <- function(data, c, k) {
  nseq <- nrow(data)
  cpartition <- split(seq_along(data[[c]]), data[[c]])
  kpartition <- split(seq_along(data[[k]]), data[[k]])
  nk <- purrr::map_int(kpartition, length)
  2 / nseq * sum(
    purrr::map_dbl(
      cpartition,
      inner_fmeasure,
      kpartition,
      nk
    )
  )
}

