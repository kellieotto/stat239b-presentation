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


save(sims, results, file = "AbadieImbensDGP_Results.Rdata")

library(dplyr)
library(tidyr)
library(ggplot2)

n.sam <- seq(1000, 10000, 1000)
alpha <- c(1, 10)
n.sim <- 1000

exact.tau <- expand.grid(sample.size = n.sam, alpha = alpha) %>%
  mutate(n0 = round(n.sam/(1+alpha)),
         n1 = sample.size - n0,
         exact.var.tau.hat = 1/n1 + 3/2 * (n1 -1) * (n0 + 8/3) / n1 / (n0 + 1) / (n0 + 2))

results <- sims %>%
  group_by(sample.size, alpha) %>%
  summarise(mean.var.boot = mean(var.boot),
            var.tau.hat = var(tau.hat),
            d.varboot.vartauhat = mean.var.boot - var.tau.hat) %>%
  left_join(exact.tau) %>%
  mutate(error.varboot = mean.var.boot / exact.var.tau.hat,
         error.var.tau.hat = var.tau.hat / exact.var.tau.hat) %>%
  gather(estimator, error.ratio, error.varboot, error.var.tau.hat) %>%
  mutate(alpha = factor(alpha),
         Estimator = ifelse(estimator == "error.varboot", 
                            "Average Bootstrap Variance Estimate",
                            "Variance of Tau Estimates"))

ggplot(results) +
  geom_line(aes(x = sample.size, y = error.ratio, color = alpha)) +
  facet_grid(~Estimator) +
  labs(x = "Sample Size",
       y = "Estimate of Variance / Limit of True Variance",
       title = "Ratio of Variance Estimate and Limit Variance \n across 1000 simulations, tau = 3",
       subtitle = "1000 simulations") +
  theme(
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave("../../fig/AndrewSimulationBias.png")

# <--- Verifying Lemma 3.1 --->#
if (FALSE) {
var(simulation[ ,"tau.hat"])


# 3.1.iii
scaled <- sqrt(n1) * (simulation[ ,"tau.hat"] - tau)
mean(scaled)
var(scaled)
}
