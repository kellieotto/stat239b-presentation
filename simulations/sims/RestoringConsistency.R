rm(list = ls())
library(parallel)
library(rlecuyer)
library(Matching)
library(boot)
library(dplyr)
# Load the functions
#path.to.functions <- "~/Dropbox/berkeley/stat239B/new/stat239b-presentation/simulations/functions/"
#invisible(sapply(list.files(path.to.functions, full.names = TRUE), source))


mcreplicate <- function(n, expr,  ncores = 35, simplify = T) {
  ans <- mclapply(integer(n), 
                  eval.parent(substitute(function(...) expr)), mc.preschedule = T,
                  mc.cores = ncores)
  if (simplify) simplify2array(ans)
  else ans
}

DGPAbadieImbens <- function(n.sam, tau, alpha) {
  # Assumption 3.1
  X <- runif(n.sam)
  
  # Assumption 3.2/3.3
  n0 <- round(n.sam / (1 + alpha))
  W <- rep(1, n.sam)
  W0 <- sample(1:n.sam, n0)
  W[W0] <- 0
  
  # Assumption 3.4
  Y <- rep(tau, n.sam)
  Y[W0] <- rnorm(n0)
  
  return(data.frame(Y = Y, W = W, X = X))
}

#Implementing multiple matches by modify ATD functions
require(RANN)
FindMatches <- function (formula, data, k=1) {
  # Checks on the function inputs
  
  # Data frame issues
  if (!is.data.frame(data)) {
    stop("Data must be in a dataframe", call. = FALSE)
  }
  if (sum(is.na(data)) > 0) 
    stop("There are missing values in the data")
  # Coerce factors
  cc <- sapply(data, is.character)
  data[cc] <- lapply(data[cc], factor)
  
  terms.formula <- terms(formula)
  mf <- model.frame(terms.formula, data)
  W <- model.response(mf)
  X <- model.matrix(terms.formula, data = mf)
  
  # This is expensive --- need to think of a way to do this without storing dm
  # Only euclidean distance is implemented
  W1 <- W == 1
  nn.idx <- nn2(X[!W1,], X[W1,], k = k)$nn.idx
  # Forms list of indices for treatment and control matches
  # Seems like the wrong data structure, but this simplifies for Abadie and Imbens'
  # bootstrap calculation
  
  match_locations <- apply(nn.idx, 1, function(z) list(which(!W1)[z]))
  
  matches <- Map(function(x, y) list(t = x, c = unname(unlist(y))),
                 which(W1), match_locations)
  
  
  return(matches)
}

CalcMatchedATE <- function(response, matches, weights) {
  diffs <- vapply(matches,
                  function(x) response[x$t] - mean(mean(response[x$c])), numeric(1))
  sum(diffs * weights)/sum(weights)
}


BootMatchedATE <- function(boot.idx, formula, data, response.name, k = 1) {
  
  bootdata <- data[boot.idx,]
  bootdata <- bootdata[order(bootdata$X),]
  reps <- group_size(group_by(bootdata, X))
  boot_unique <- as.data.frame(unique(as.matrix(bootdata), byrow = TRUE))
  
  response <- boot_unique[, response.name]
  bootx <- boot_unique[, setdiff(names(data), response.name)]
  
  matches <- FindMatches(formula, boot_unique, k = k)
  
  weights <- reps[boot_unique$W == 1]
  
  boot.var <- (Match(Y = boot_unique$Y, Tr = boot_unique$W, X = boot_unique$X, M = k)$se)^2
  boot.est <- CalcMatchedATE(response, matches, weights = weights)
  c(boot.est, boot.var)
}

BootMatchedATEWrapper <- function(data, idx, k) {
  boot.est <- BootMatchedATE(idx, W ~ X, data, "Y", k = k)
  bootdd <- data[idx,]
  boot.var <- (Match(Y = bootdd$Y, Tr = bootdd$W, X = bootdd$X, M = k)$se)^2
  
  return(c(boot.est, boot.var))
}


StudentCovers <- function(n, type = "perc") {
  tau <- 1  
  dd <- DGPAbadieImbens(n, tau, 1)
  k <- ceiling(log(n)/3)
  
  bootdd <- boot(dd, BootMatchedATEWrapper, R = 1000, stype = "i", k =k)
  require(Matching)
  var.t0 <- (Match(Y = dd$Y, Tr = dd$W, X = dd$X, M = k)$se)^2
  if(type == "student") ci <- boot.ci(bootdd, var.t0 = var.t0, type = "stud")$student[4:5]
  if(type == "perc")   ci <- boot.ci(bootdd, var.t0 = var.t0, type = "perc")$perc[4:5]
  covers <- ci[1] < tau & tau < ci[2]
  return(covers)
}

system.time(tmp <- mean(mcreplicate(70, StudentCovers(100, type = "perc"))))
tmp

#res <- mclapply(c(50, 75, 100, 200), function(n) mean(replicate(2000, Efroncovers(n))),
#a         mc.preschedule = FALSE, mc.cores = parallel::detectCores()-1)
#mean(reps)



# #Residual bootstrap
# ResidualCover <- function(n){
#   tau <- 1  
#   dd <- DGPAbadieImbens(n, tau, 0.5)
#   matches <- FindMatches(W~X, data = dd)
#   response <- dd$Y
#   diffs <- vapply(matches,
#                   function(x) response[x$t] - mean(response[x$c]), numeric(1))
#   
# }



#Abadie and Imbens Normal approximation
AIcovers <- function(n, alpha= 0.05) {
  tau <- 1  
  dd <- DGPAbadieImbens(n, tau, 1)
  k <- ceiling(log(n)/3)
  
  
  sekhon <- Match(Y = dd$Y, Tr = dd$W, X = dd$X, M = k)
  ci <- c(sekhon$est +  qnorm(alpha/2)*sekhon$se, 
          sekhon$est +  qnorm(1-alpha/2)*sekhon$se)
  covers <- ci[1] < tau & tau < ci[2]
  return(covers)
}

system.time(tmp <- mean(mcreplicate(35*20, AIcovers(100))))

# ASKanalyticcovers <- function(n, alpha = 0.5) {
#   tau <- 1  
#   dd <- DGPAbadieImbens(n, tau, 0.5)
#   matches <- FindMatches(W~X, data = dd)
#   response <- dd$Y
#   diffs <- vapply(matches,
#                   function(x) response[x$t] - mean(response[x$c]), numeric(1))
#   
#   se <- sd(diffs)/sqrt(n)
#   est <- mean(diffs)
#   ci <- c(est +  qnorm(alpha/2)*se, 
#           est +  qnorm(1-alpha/2)*se)
#   covers <- ci[1] < tau & tau < ci[2]
# }

simulate <- function(n) {
  cat(paste("n = ", n, ":", "Time is "))
  cat(paste(as.character(Sys.time()), "\n"))
  n_ai <- 2000
  n_efron <- 210
  ai <- mean(mcreplicate(n_ai, AIcovers(n)))
  efron <- mean(mcreplicate(n_efron, StudentCovers(n, type = "perc")))
  sd_ai <- sqrt(ai * (1-ai)/n_ai)
  sd_efron <- sqrt(efron * (1-efron)/n_efron)
  return(data.frame(est = c(ai,  efron), sd = c(sd_ai, sd_efron), n = n, type = c("ai", "efron")))
}  
ns <- seq(1000, 2000, 50)


coverage_ai <-
  sapply(c(seq(100, 400, 100), 1000),
         function(n_ai) {
           cat(paste("n = ", n_ai, ":", "Time is "))
           cat(paste(as.character(Sys.time()), "\n"))
           (mean(mcreplicate(9975*5, AIcovers(n_ai))))
         }
  )


coverage_ai2 <-
  sapply(c(seq(500, 900, 100)),
         function(n_ai) {
           cat(paste("n = ", n_ai, ":", "Time is "))
           cat(paste(as.character(Sys.time()), "\n"))
           (mean(mcreplicate(9975*5, AIcovers(n_ai))))
         }
  )



coverage_efron3 <-
  sapply(c(seq(100, 400, 100), 1000),
         function(n) {
           cat(paste("n = ", n, ":", "Time is "))
           cat(paste(as.character(Sys.time()), "\n"))
           (mean(mcreplicate(280, StudentCovers(n,type = "perc"))))
         }
  )


coverage_efron_alt <-
  sapply(c(seq(500, 900, 100)),
         function(n) {
           cat(paste("n = ", n, ":", "Time is "))
           cat(paste(as.character(Sys.time()), "\n"))
           (mean(mcreplicate(980, StudentCovers(n,type = "perc"))))
         }
  )

coverage_efron_alt2 <-
  sapply(c(seq(100, 1000, 100)),
         function(n) {
           cat(paste("n = ", n, ":", "Time is "))
           cat(paste(as.character(Sys.time()), "\n"))
           (mean(mcreplicate(980, StudentCovers(n,type = "perc"))))
         }
  )
#n 280, 700, 280

coverage_efron_agg <- mapply(
  function(x,y,z) (280*x + 700 * y + 280 * z)/(280+700+280),
  coverage_efron, coverage_efron2, coverage_efron3)

#system.time(coverage_efron <- c(coverage_efron,  mean(mcreplicate(280, StudentCovers(1000,type = "perc")))))




dd <- data.frame(coverage = c(coverage_ai, coverage_ai2, coverage_efron_agg[1:4], coverage_efron_alt, coverage_efron_agg[5]),
                 method = rep(c("Asymptotic", "Bootstrap"), times = c(10,10)),
                 n = c(seq(100, 1000, 100),seq(100,1000,100)), 
                 n_sim = rep(c(9975*5, 1260), times = c(10,10))
)

dd$sd <- sqrt(dd$coverage * (1 - dd$coverage)/dd$n_sim)

require(ggplot2)                 
ggplot(dd, mapping = aes(x = n, y = coverage, ymin = coverage - 2 * dd$sd,
                         ymax = coverage + 2*sd, col = method, fill = method)) +
  geom_line() + theme_bw() + geom_ribbon(alpha = 0.2) + 
  geom_hline(yintercept = 0.95, linetype = 2) + ylab("True coverage") +
  ggtitle("True coverage for asymptotic and bootstrap CIs with nominal 0.95 level")
ggsave("MainFigure.PDF")


sim_res <- sapply(ns, simulate, simplify = FALSE)
sim_res <- do.call(rbind, sim_res)
save(sim_res, file = "simresults.Rda")
require(ggplot2)
require(reshape2)

ggplot(data = as.data.frame(sim_res), 
       aes(x = n, y = est,
           ymin = est - 2 * sd, ymax = est + 2*sd, fill = type)) +
  geom_line(aes(col= type)) + geom_ribbon(alpha = 0.2) +
  theme_bw() + geom_hline(yintercept = 0.95) + xlab("n") + ylab("Coverage") + 
  ggtitle("Comparing the Bootstrap and AI (2006) for multiple matches")
ggsave("ComparingAI.pdf")





system.time(tmp <- simulate(100))

