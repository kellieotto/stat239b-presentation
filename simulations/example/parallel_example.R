library(parallel)
library(rlecuyer)

# Load the functions.  Change path as necessary
path.to.functions <- "../functions/"
invisible(sapply(list.files(path.to.functions), source))

# Parallel Backend Setup
RNGkind("L'Ecuyer-CMRG") # So we actually have random data
n.cores <- detectCores() - 2
cl <- makeCluster(n.cores, type="SOCK")
clusterExport(cl=cl, 
              varlist = c("BootFindMatches", 
                          "BootMatchedATE", 
                          "Bootstrap", 
                          "CalcBootstrapVar", 
                          "CalcMatchedATE", 
                          "DGPAbadieImbens", 
                          "FindMatches", 
                          "SimulateDGP"))
clusterEvalQ(cl, {
  library(RANN)
  n.sam <- 100 # sample size in each iteration
  alpha <- 2
  tau <- 5})

n.sim <- 600 # number of times to simulate
simulation <- parReplicate(cl, n.sim, 
  SimulateDGP(DGPAbadieImbens, n.sam = n.sam, alpha = alpha, tau = tau))

# Free up resources on your computer again
stopCluster(cl)