# Utility function that finds additional matches during bootstrapping
# Exploits the fact that units are duplicates in a bootstrap sample

BootFindMatches <- function(boot.idx, formula, data) {
  matches <- FindMatches(formula, data)
  dup <- lapply(boot.idx, function(x) which(boot.idx == x))
  lapply(matches, function(x) list(t = x$t, c = dup[[x$c]]))
}
