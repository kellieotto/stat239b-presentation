# Data generating process used in section 3 of Abadie and Imbens

DGPAbadieImbens <- function(n.sam, tau, alpha) {
  # Assumption 3.1
  X <- runif(n.sam)
  
  # Assumption 3.2/3.3
  n0 <- round(n.sam / (1 + alpha))
  W <- rep(1, n.sam)
  W0 <- sample(1:n.sam, n0)
  W[W0] <- 0
  
  # Assumption 3.4
  Y <- rep(tau, n.sam)
  Y[W0] <- rnorm(n0)
  
  return(data.frame(Y = Y, W = W, X = X))
}