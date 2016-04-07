library(parallel)
library(rlecuyer)
library(plyr)

# Load the functions
wd <- Sys.getenv("R_WD")
setwd(wd)
path.to.functions <- "../functions/"
functions <- list.files(path.to.functions)
invisible(sapply(paste0(path.to.functions, functions), source))

# parameters to tweak
n.sam <- seq(500, 10000, 500)
alpha <- 1:10
tau <- 3
n.sim <- 1000

# Parallel Backend Setup
RNGkind("L'Ecuyer-CMRG") # So we actually have random data
n.cores <- detectCores()
core.message <- sprintf("Using %i cores", n.cores)
print(core.message)
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
                          "tau"
                          )
              )
                          

# Simulation Call

simulation <- alply(expand.grid(n = n.sam, a = alpha), 1,
  function(x) {
    env <- environment()
    n <- x$n; a <- x$a
    clusterExport(cl = cl, varlist = c("n", "a"), envir = env)
    print(sprintf("Now working on n = %i and a = %i", n, a))
    parReplicate(cl, n.sim, SimulateDGP(DGPAbadieImbens, 
                                        n.sam = n, 
                                        alpha = a, 
                                        tau = tau))
    }
)

# Save Results
save(simulation, file = "AbadieImbensDGP_Results.RData")

# Free up resources on your computer again
stopCluster(cl)

# <--- Verifying Lemma 3.1 --->#
if (FALSE) {
var(simulation[ ,"tau.hat"])
n0 <- round(n.sam/(1+alpha))
n1 <- n.sam - n0
exact.var.tau.hat <- 1/n1 + 3/2 * (n1 -1) * (n0 + 8/3) / n1 / (n0 + 1) / (n0 + 2)

# 3.1.iii
scaled <- sqrt(n1) * (simulation[ ,"tau.hat"] - tau)
mean(scaled)
var(scaled)
}
