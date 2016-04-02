Bootstrap <- function (x, B, theta, ...) {
  n <- length(x)
  boot.idx <- sample(x, size = n * B, replace = TRUE)
  boot.sample <- matrix(boot.idx, nrow = B)
  thetastar <- apply(boot.sample, 1, theta, ...)

  return(thetastar)
}