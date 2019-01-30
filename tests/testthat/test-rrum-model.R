context("test-rrum-model")

# Set a seed for reproducibility
set.seed(888)

# Setup Parameters
N = 15   # Number of Examinees / Subjects
J = 10   # Number of Items
K = 2    # Number of Skills / Attributes

# Simulate identifiable Q matrix
Q = sim_q_matrix(J, K)

# Penalties for failing to have each of the required attributes
rstar  = .5 * Q

# The probabilities of answering each item correctly for individuals 
# who do not lack any required attribute
pistar = rep(.9, J)

# Latent Class Probabilities
pis = c(.1, .2, .3, .4)

# Generate latent attribute profile with custom probability (N subjects by K skills)
subject_alphas = sim_subject_attributes(N, K, prob = pis)

# Simulate rrum items
rrum_items = simcdm::sim_rrum_items(Q, rstar, pistar, subject_alphas)

test_that("Verify Bad Input detected", {

  ## Verify error catches ----
  
  # R-level main
  expect_error(rrum(rrum_items, Q, as = c(1, 2)), info = "Detected Bad `as` data")
  expect_error(rrum(rrum_items, Q, ag = c(3, 2)), info = "Detected Bad `ag` data")
  expect_error(rrum(rrum_items, Q, bs = c(4, 2)), info = "Detected Bad `bs` data")
  expect_error(rrum(rrum_items, Q, bg = c(5, 6)), info = "Detected Bad `bg` data")
  
  # C++-level via helper
  expect_error(rrum(rrum_items, Q, delta0 = c(4, 3, 2)), info = "Detected Bad `delta0` data")
  expect_error(rrum(rrum_items[,-1], Q), info = "Detected Bad `Y` and `Q` data")
  
})

test_that("Verify Model Reproducibility", {
  
  ## Estimate two models ----
  
  # Set a seed for reproducibility
  set.seed(518)
  model1 = rrum(rrum_items, Q, chain_length = 100)
  
  # Set the same seed
  set.seed(518)
  model2 = rrum(rrum_items, Q, chain_length = 100)
  
  expect_equal(model1, model2,
               info = "Models are reproducible under same seed.")
  
})

