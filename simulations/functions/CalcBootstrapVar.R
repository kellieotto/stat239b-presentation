CalcBootstrapVar <- function(tau.boot, tau.hat) {
  mean((tau.boot - tau.hat)^2)
}