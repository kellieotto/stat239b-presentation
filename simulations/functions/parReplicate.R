# Wrapper function for parSapply that emulates vanilla replicate

# cl - cluster object created by parallel package
# n - integer, the number of replications
# expr - the expression (a language object, usually a call) to evaluate repeatedly
# simplify - logical should the result be simplified to a matrix if possible?

parReplicate <- function(cl, n, expr, simplify = T) {
  ans <- parLapply(cl, integer(n), 
    eval.parent(substitute(function(...) expr)))
  if (simplify) do.call(rbind, lapply(ans, do.call, what = cbind))
  else ans
}