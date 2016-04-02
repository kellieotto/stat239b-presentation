library(parallel)
library(rlecuyer)

# Load the functions
path.to.functions <- "../functions/"
invisible(sapply(list.files(path.to.functions), source))

# Parameters to tweak
n.sam <- 100
alpha <- 2
tau <- 5
n.sim <- 6000

# <------------------------- Parallel Backend Setup ------------------------->#

# No need to modify this
RNGkind("L'Ecuyer-CMRG") # So we actually have random data
n.cores <- detectCores() - 2
cl <- makeCluster(n.cores, type="SOCK")

# This passes global variables to each cluster
# Most likely you'll have to add parameters from your dgp to
# varlist (unless you only use n.sam, alpha, and tau)
clusterExport(cl=cl, 
              varlist = c("BootFindMatches", 
                          "BootMatchedATE", 
                          "Bootstrap", 
                          "CalcBootstrapVar", 
                          "CalcMatchedATE", 
                          "DGPAbadieImbens", 
                          "FindMatches", 
                          "SimulateDGP",
                          "nn2", # nn2 from RANN package
                          "n.sam",
                          "alpha",
                          "tau"
              )
)

# <------------------------- Parallel Backend Setup ------------------------->#

simulation <- parReplicate(cl, n.sim, 
  SimulateDGP(DGPAbadieImbens, n.sam = n.sam, alpha = alpha, tau = tau))

# Free up resources on your computer again
stopCluster(cl)


