if (!require(dplyr)) install.packages("dplyr")
library(dplyr)

CalculateK <- function(matches) {
  do.call(rbind, lapply(matches, do.call, what = data.frame)) %>%
    group_by(t) %>%
    mutate(size.J = n()) %>%
    group_by(c) %>%
    summarise(K = sum(1/size.J)) %>%
    as.data.frame
}