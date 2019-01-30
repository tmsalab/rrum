#' Gibbs sampler to estimate the rRUM
#'
#' Obtains samples from posterior distributon for the reduced Reparametrized
#' Unified Model (rRUM).
#' 
#' @param Y            A `matrix` with N rows and J columns, where N reperesnts
#'                     the number of individuals and J the number of items. 
#'                     `Y` indicates the indviduals' responses to each of the
#'                     items.
#' @param Q            A `matrix` with J rows and K columns indicating which
#'                     attributes are required to answer each of the items. 
#'                     An entry of 1 indicates attribute k is required to
#'                     answer item j.  An entry of one indicates attribute `k`
#'                     is not required.
#' @param chain_length A `numeric` indicating the number of iterations of Gibbs
#'                     sampler to be run.  Default is set to 10000.
#' @param as           A `numeric`, parameter for the prior distribution of
#'                     pistar.  High values as encourage higher values of
#'                     pistar and lower values of rstar.
#' @param bs           A `numeric`, parameter for the prior distribution of
#'                     pistar.  High values as encourage lower values of
#'                     pistar and higher values of rstar.
#' @param ag           A `numeric`, parameter for the prior distribution of
#'                     rstar.  High values as encourage higher values of rstar.
#' @param bg           A `numeric`, parameter for the prior distribution of
#'                     pistar.  High values as encourage lower values of rstar.
#' @param delta0       A `vector`, parameters for the Dirichlet prior on pi.
#' 
#' @return 
#' A `list` that contains
#' 
#' - `PISTAR`: A `matrix` where each column represents one draw from the
#'            posterior distribution of pistar.
#' - `RSTAR`: A \eqn{J x K x chain_length} `array` where `J` reperesents the
#'           number of items, and `K` represents the number of attributes.
#'           Each slice represents one draw from the posterior distribution
#'           of `rstar`.
#' - `PI`: A `matrix` where each column reperesents one draw from the posterior
#'        distribution of `pi`.
#' - `ALPHA`: An \eqn{N x K x chain_length} `array` where `N` reperesents the
#'           number of individuals, and `K` represents the number of
#'           attributes. Each slice represents one draw from the posterior
#'           distribution of `alpha`.
#'           
#' @author 
#' Steven Andrew Culpepper, Aaron Hudson, and James Joseph Balamuta
#' 
#' @seealso 
#' [simcdm::sim_rrum_items()]
#' 
#' @export
#' @template rrum-example
#' @template rrum-references
rrum = function(Y, Q, chain_length = 10000L, 
                as = 1, bs = 1, ag = 1, bg = 1, 
                delta0 = rep(1, 2 ^ ncol(Q))) {
  
  if(length(as) != 1 | length(ag) != 1 | length(bs) != 1 | length(bg) != 1) {
    stop("as, ag, bs, and bg must all be numeric and of length 1")
  }
  
  rrum_helper(Y, Q, delta0, chain_length, as, bs, ag, bg)
}
