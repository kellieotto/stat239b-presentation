BootMatchedATE <- function(formula, data, response.name, boot.idx) {
  response <- data[response.name]
  data[response.name] <- NULL
  matches <- FindMatches(formula, data)
  matches <- DupMatchInBootstrap(matches, boot.idx)
  CalcMatchedATE(response, matches)
}