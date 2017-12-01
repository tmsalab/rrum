#' @useDynLib rrum, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom simcdm sim_rrum
#' @details 
#' Implemention of the reduced Reparametrized Unified Model
#' @template rrum-references
#' @aliases rrum-package
"_PACKAGE"

### Import from simcdm. Move to core package

#' @inherit simcdm::pi_reference
#' @importFrom simcdm pi_reference
#' @export
pi_reference = simcdm::pi_reference

#' @inherit simcdm::bijectionvector
#' @importFrom simcdm bijectionvector
#' @export
bijectionvector = simcdm::bijectionvector
