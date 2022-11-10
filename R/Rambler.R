#------------------------------------------------
#' @title Rambler
#'
#' @description Malaria infections often contain multiple haplotypes. When
#'   sequenced at several points in time, these haplotypes may be present for a
#'   period beforing cleaing as parasite density decreases. We may also fail to
#'   detect a haplotype at a given point in time due to imperfect sensitivity of
#'   our sequencing. Rambler takes these factors into account and produces
#'   estimates of the time(s) at which each individual in a longitudinal cohort
#'   became infected, and for how long this infection lasted.
#'
#' @docType package
#' @name Rambler
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib Rambler, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("Rambler", libpath)
}
