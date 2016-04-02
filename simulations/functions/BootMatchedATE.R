BootMatchedATE <- function(boot.idx, formula, data, response.name) {
  response <- data[boot.idx, response.name]
  data <- data[boot.idx, setdiff(names(data), response.name)]
  matches <- BootFindMatches(boot.idx, formula, data)
  CalcMatchedATE(response, matches)
}