# Utility function that finds additional matches during bootstrapping
# Exploits the fact that units are duplicates in a bootstrap sample

DupMatchInBootstrap <- function(matches, boot.idx) {
  dup <- lapply(boot.idx, function(x) which(boot.idx == x))
  lapply(matches, function(x) list(t = x$t, c = dup[[x$c]]))
}
