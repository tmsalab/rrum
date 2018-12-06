#' @useDynLib rrum, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom simcdm sim_rrum_items
#' @details 
#' Implemention of the reduced Reparametrized Unified Model
#' @template rrum-references
#' @aliases rrum-package
"_PACKAGE"

### Import from simcdm. Move to core package

#' @inherit simcdm::sim_alpha_matrix
#' @importFrom simcdm sim_alpha_matrix
#' @export
sim_alpha_matrix = simcdm::sim_alpha_matrix

#' @inherit simcdm::bijectionvector
#' @importFrom simcdm bijectionvector
#' @export
bijectionvector = simcdm::bijectionvector
