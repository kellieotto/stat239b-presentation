library(parallel)
library(rlecuyer)

# Load the functions
path.to.functions <- "../functions/"
functions <- list.files(path.to.functions)
invisible(sapply(paste0(path.to.functions, functions), source))

# parameters to tweak
n.sam <- 1000
alpha <- 0.1
tau <- 3
n.sim <- 10000

# Parallel Backend Setup
RNGkind("L'Ecuyer-CMRG") # So we actually have random data
n.cores <- detectCores()
core.message <- sprintf("Using %i cores", n.cores)
print(core.message)
cl <- makeCluster(n.cores, type="SOCK")
clusterExport(cl=cl, 
              varlist = c("DGPAbadieImbens", 
                          "FindMatches", 
                          "nn2", # nn2 from RANN package
                          "tau",
                          "alpha",
                          "n.sam",
                          "SimulateK",
                          "CalculateK"
              )
)
clusterEvalQ(cl, {library(dplyr)})

simulations <- parReplicate(cl, n.sim, simplify = FALSE,
  expr = SimulateK(DGPAbadieImbens, n.sam = n.sam, tau = tau, alpha = alpha))

pn.k <- CalcDistK(simulations, i = 888, n.sam)
