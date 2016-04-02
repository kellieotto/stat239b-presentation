CalcMatchedATE <- function(response, matches) {
  diffs <- vapply(matches,
                  function(x) response[x$t] - mean(response[x$c]), numeric(1))
  mean(diffs)
}