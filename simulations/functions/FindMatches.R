library(RANN)

FindMatches <- function (formula, data) {
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
  
  terms.formula <- terms(formula)
  mf <- model.frame(terms.formula, data)
  W <- model.response(mf)
  X <- model.matrix(terms.formula, data = mf)
  
  # This is expensive --- need to think of a way to do this without storing dm
  # Only euclidean distance is implemented
  W1 <- W == 1
  nn.idx <- nn2(X[!W1,], X[W1,], k = 1)$nn.idx[,1]
  # Forms list of indices for treatment and control matches
  # Seems like the wrong data structure, but this simplifies for Abadie and Imbens'
  # bootstrap calculation
  matches <- Map(function(x, y) list(t = x, c = y),
                 which(W1), which(!W1)[nn.idx])
  
  return(matches)
}