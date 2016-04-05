# Data generating process with limited covariate overlap
# We keep assumptions 3.2-3.4 from Abadie-Imbens
# Change the covariate distribution: for mu > 0,
# X ~ N(0,1) for the treated group
# X ~ N(mu,1) for the control group
# Y(0) ~ N(0,1)
# Y(1) = tau
# Limited covariate overlap will mean that few controls, with small X values, will be matched many times.
# In bootstrap samples where these small controls don't appear, matches will be poor.
# In bootstrap samples where they do appear, their weight K will much larger than in the original sample.

DGPLimitedOverlap <- function(n.sam, tau, alpha, mu) {
  
  # Assumption 3.2/3.3
  n0 <- round(n.sam / (1 + alpha))
  W <- rep(1, n.sam)
  W0 <- sample(1:n.sam, n0)
  W[W0] <- 0

  # Generate X
  X <- rnorm(n.sam)
  X[W0] <- X[W0]+mu
  
  # Assumption 3.4
  Y <- rep(tau, n.sam)
  Y[W0] <- rnorm(n0)
#  Y <- X + rnorm(n.sam)
#  Y[W0] <- Y[W0] - tau
  
  return(data.frame(Y = Y, W = W, X = X))
}