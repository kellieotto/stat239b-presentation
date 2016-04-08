SimulateK <- function(dgp, formula = W ~ X, response.name = "Y", ...) {
  
  data <- dgp(...)
  matches <- FindMatches(formula, data)
  CalculateK(matches)
  
}

SimulateBootK <- function(dgp, formula = W ~ X, response.name = "Y", ...) {
  Bootstrap <- function (x, B, theta, ...) {
    n <- length(x)
    boot.idx <- sample(x, size = n * B, replace = TRUE)
    boot.samples <- matrix(boot.idx, nrow = B)
    boot.samples <- split(boot.samples, row(boot.samples))
    
    return(thetastar)
  }
}

CalcDistK <- function(simulations, i = "all", n.sam) {
  if (!is.character(i) && i > 0 && i %% 1 == 0) {
    return(sapply(simulations, function(x) ifelse(i %in% x$c, x$K[which(x$c==i)], 0)))
  } else if (i == "all") {
    return(sapply(simulations, function(x) ifelse(1:n.sam %in% x$c, x$K, 0)))
  } else stop("i incorrectly specified")
  
}