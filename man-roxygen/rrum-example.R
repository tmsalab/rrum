#' @examples 
#' # Set seed for reproducibility
#' set.seed(217)
#' 
#' ## Define Simulation Parameters
#' 
#' N = 1000 # Number of Individuals
#' J = 6    # Number of Items
#' K = 2    # Number of Attributes
#' 
#' # Matrix where rows represent attribute classes
#' As = sim_alpha_matrix(K) 
#' 
#' # Latent Class probabilities
#' pis = c(.1, .2, .3, .4) 
#' 
#' # Q Matrix
#' Q = rbind(c(1, 0),
#'           c(0, 1),
#'           c(1, 0),
#'           c(0, 1),
#'           c(1, 1),
#'           c(1, 1)
#'     )
#'     
#' # The probabiliies of answering each item correctly for individuals 
#' # who do not lack any required attribute
#' pistar = rep(.9, J)
#' 
#' # Penalties for failing to have each of the required attributes
#' rstar  = .5 * Q
#' 
#' # Randomized alpha profiles
#' alpha  = As[sample(1:(K ^ 2), N, replace = TRUE, pis),]
#' 
#' # Simulate data
#' rrum_items = sim_rrum_items(Q, rstar, pistar, alpha)
#' 
#' \dontrun{
#' # Note: This portion of the code is computationally intensive.
#' 
#' # Recover simulation parameters with Gibbs Sampler
#' Gibbs.out = rrum(rrum_items, Q)
#' 
#' # Iterations to be discarded from chain as burnin
#' burnin = 1:5000 
#' 
#' # Calculate summarizes of posterior distributions
#' rstar.mean  = with(Gibbs.out, apply(RSTAR[,,-burnin], c(1, 2), mean))
#' pistar.mean = with(Gibbs.out, apply(PISTAR[,-burnin], 1, mean))
#' pis.mean    = with(Gibbs.out, apply(PI[,-burnin], 1 ,mean))
#' }

