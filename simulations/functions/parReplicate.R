# Wrapper function for parSapply that emulates vanilla replicate

# cl - cluster object created by parallel package
# n - integer, the number of replications
# expr - the expression (a language object, usually a call) to evaluate repeatedly
# simplify - logical or character string; should the result be simplified to a
#   vector, matrix or higher dimensional array if possible?

parReplicate <- function(cl, n, expr, simplify = "array") {
  parSapply(cl, integer(n), 
            eval.parent(substitute(function(...) expr)), 
            simplify = simplify)
}