library(Rcpp)
sourceCpp('CalcPDist.cpp')

FindMatches <- function (formula, data, distance.method, ...) {
  
  # Save the function call
  mcall <- match.call()
  
  
  # Checks on the function inputs
  
  # Data frame issues
  if (!is.data.frame(data)) {
    stop("Data must be in a dataframe", call. = FALSE)
  }
  if (sum(is.na(data)) > 0) 
    stop("There are missing values in the data")
  # Coerce factors
  cc <- sapply(test, is.character)
  test[cc] <- lapply(test[cc], factor)
  
  
  if (!is.numeric(distance)) {
    fn1 <- paste("Calc", distance, "Dist", sep = "")
    if (!exists(fn1)) 
      stop(distance, "not supported.")
  }
  if (is.numeric(distance)) {
    fn1 <- "distance2user"
  }
  
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  mf <- model.frame(tt, data)
  treat <- model.response(mf)
  X <- model.matrix(tt, data = mf)
  
  # This is expensive
  distance <- as.matrix(dist(X, method = distance.method))
  
    
  out2$call <- mcall
  out2$model <- out1$model
  out2$formula <- formula
  out2$treat <- treat
  if (is.null(out2$X)) {
    out2$X <- X
  }
  out2$distance <- distance

  nn <- matrix(0, ncol = 2, nrow = 4)
  nn[1, ] <- c(sum(out2$treat == 0), sum(out2$treat == 1))
  nn[2, ] <- c(sum(out2$treat == 0 & out2$weights > 0), 
               sum(out2$treat == 1 & out2$weights > 0))
  nn[3, ] <- c(sum(out2$treat == 0 & out2$weights == 0 & out2$discarded == 0), 
               sum(out2$treat == 1 & out2$weights == 0 & out2$discarded == 0))
  dimnames(nn) <- list(c("All", "Matched", "Unmatched"), 
                       c("Control", "Treated"))
  out2$nn <- nn
  return(out2)
}