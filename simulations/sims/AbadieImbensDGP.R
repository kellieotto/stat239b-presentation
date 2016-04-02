library(parallel)
library(rlecuyer)

# Load the functions
path.to.functions <- "../functions/"
invisible(sapply(list.files(path.to.functions), source))

# parameters to tweak
n.sam <- 100
alpha <- 2
tau <- 5
n.sim <- 6000

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
                          "SimulateDGP",
                          "nn2", # nn2 from RANN package
                          "n.sam",
                          "alpha",
                          "tau"
                          )
              )
                          

# Simulation Call
simulation <- parReplicate(cl, n.sim, 
                           SimulateDGP(DGPAbadieImbens, n.sam = n.sam, alpha = alpha, tau = tau))

# Free up resources on your computer again
stopCluster(cl)

# <--- Verifying Lemma 3.1 --->#

var(simulation[ ,"tau.hat"])
n0 <- round(n.sam/(1+alpha))
n1 <- n.sam - n0
exact.var.tau.hat <- 1/n1 + 3/2 * (n1 -1) * (n0 + 8/3) / n1 / (n0 + 1) / (n0 + 2)

# 3.1.iii
scaled <- sqrt(n1) * (simulation[ ,"tau.hat"] - tau)
mean(scaled)
var(scaled)
