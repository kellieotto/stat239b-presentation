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
n.sim <- 100000

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
pn.k <- data.frame(k = CalcDistK(simulations, i = 1, n.sam), subtitle = "Empirical Average")

boot.simulations <- SimulateBootK(DGPAbadieImbens, n.boot = 10000, n.sam = n.sam, tau = tau, alpha = alpha)
pn.k2 <- data.frame(k = CalcDistK(boot.simulations, i = 1, n.sam), subtitle = "Bootstrap Average")

p.k <- rbind(pn.k, pn.k2) %>% filter(k != 0)
ggplot(p.k) +
  geom_histogram(aes(x = k, y = ..density..), binwidth = .25) +
  facet_grid(~subtitle) +
  labs(title = "Distribution of K, W = 0, tau = 3, alpha = 0.25, 100000 Simulations",
       y = "Density") +
  theme(
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave("../../fig/KDistribution.png")


