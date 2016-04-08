SimulateK <- function(dgp, formula = W ~ X, response.name = "Y", ...) {
  
  data <- dgp(...)
  matches <- FindMatches(formula, data)
  CalculateK(matches)
  
}
Bootstrap <- function (x, B, theta, ...) {
  n <- length(x)
  boot.idx <- sample(x, size = n * B, replace = TRUE)
  boot.samples <- matrix(boot.idx, nrow = B)
  boot.samples <- split(boot.samples, row(boot.samples))
  
  return(thetastar)
}

BootK <- function(boot.idx, formula, data, response.name) {
  response <- data[boot.idx, response.name]
  data <- data[boot.idx, setdiff(names(data), response.name)]
  matches <- BootFindMatches(boot.idx, formula, data)
  CalculateK(matches)
}

SimulateBootK <- function(dgp, n.boot = 100, 
                          formula = W ~ X, response.name = "Y", ...) {
  data <- dgp(...)
  K.boot <- Bootstrap(1:nrow(data), B = n.boot, theta = BootK, 
                      formula = formula, data = data, response.name = response.name)
  return(K.boot)
}

CalcDistK <- function(simulations, i = "all", n.sam) {
  if (!is.character(i) && i > 0 && i %% 1 == 0) {
    return(sapply(simulations, function(x) ifelse(i %in% x$c, x$K[which(x$c==i)], 0)))
  } else if (i == "all") {
    return(sapply(simulations, function(x) ifelse(1:n.sam %in% x$c, x$K, 0)))
  } else stop("i incorrectly specified")
  
}