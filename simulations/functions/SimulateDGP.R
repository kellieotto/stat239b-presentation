# dgp - function, data generating process that produces a dataframe
# n.boot - integer, number of bootstrap samples to take
# formula - formula, specifies the propensity score W ~ covariates
# response.name - character vector of length 1, denotes name of response in dgp
# ... arguments to be passed to dgp

SimulateDGP <- function(dgp, n.boot = 100, 
                        formula = W ~ X, response.name = "Y", ...) {
  data <- dgp(...)
  matches <- FindMatches(formula, data)
  tau.hat <- CalcMatchedATE(data[[response.name]], matches)
  tau.boot <- Bootstrap(1:nrow(data), n.boot, BootMatchedATE, 
                        formula, data, response.name)
  var.boot <- CalcBootstrapVar(tau.boot, tau.hat)
  return(list(tau.hat = tau.hat, tau.boot = mean(tau.boot), var.boot = var.boot))
}
