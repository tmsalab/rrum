#' @examples 
#' ###################################
#' ### Generate random data from rRUM
#' ### and recover parameters with
#' ### Gibbs sampler
#' ###################################
#' set.seed(217)
#' # Define Simulation Parameters
#' N <- 1000 # number of individuals
#' J <- 6 # number of items
#' K <- 2 # number of attributes
#' As <- pi_reference(K) # matrix where rows represent attribute classes
#' pis <- c(.1, .2, .3, .4) # pi
#' Q <- rbind(c(1,0),
#'            c(0,1),
#'            c(1,0),
#'            c(0,1),
#'            c(1,1),
#'            c(1,1))
#' pistar <- rep(.9, J) # pistar
#' rstar <- .5 * Q # rstar
#' alpha <- As[sample(1:(K ^ 2), N, replace = TRUE, pis),]
#' # Simulate data
#' Y <- simrRUM(Q, rstar, pistar, alpha)
#' # Recover simulation parameters with Gibbs Sampler
#' Gibbs.out <- rRUM_Gibbs(Y, Q)
#' burnin <- 1:5000 # iterations to be discarded from chain as burnin
#' rstar.mean <- with(Gibbs.out, apply(RSTAR[,,-burnin], c(1, 2), mean))
#' pistar.mean <- with(Gibbs.out, apply(PISTAR[,-burnin], 1, mean))
#' pis.mean <- with(Gibbs.out, apply(PI[,-burnin], 1 ,mean))
