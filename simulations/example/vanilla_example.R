# Load the functions.  Change path as necessary
path.to.functions <- "../functions/"
invisible(sapply(list.files(path.to.functions), source))

n.sim <- 60 # number of times to simulate
n.sam <- 100 # sample size in each iteration
alpha <- 2
tau <- 5

simulations <- replicate(n.sim, 
  SimulateDGP(DGPAbadieImbens, n.sam = n.sam, alpha = alpha, tau = tau))