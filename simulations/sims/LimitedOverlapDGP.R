library(parallel)
library(rlecuyer)
library(dplyr)
library(reshape2)
library(ggplot2)


system.time({
  
# Load the functions
setwd("../functions/")
invisible(sapply(list.files(), source))
setwd("../sims/")

# parameters to tweak
n.sam <- 1000
alpha <- 2
tau <- 5
n.sim <- 1000
mu <- seq(0, 5, by = 0.5)

# Parallel Backend Setup
RNGkind("L'Ecuyer-CMRG") # So we actually have random data
n.cores <- detectCores()-2
cl <- makeCluster(n.cores, type="SOCK")
clusterExport(cl=cl, 
              varlist = c("BootFindMatches", 
                          "BootMatchedATE", 
                          "Bootstrap", 
                          "CalcBootstrapVar", 
                          "CalcMatchedATE", 
                          "DGPLimitedOverlap", 
                          "FindMatches", 
                          "SimulateDGP",
                          "nn2", # nn2 from RANN package
                          "n.sam",
                          "alpha",
                          "tau",
                          "mu"
              )
)


# Simulation Call
simulation <- lapply(mu, function(mm){
                          parReplicate(cl, n.sim, 
                           SimulateDGP(DGPLimitedOverlap, n.sam = n.sam, alpha = alpha, tau = tau, mu = mm))
              })

# Free up resources on your computer again
save(simulation, file = "LimitedOverlapDGP_Results.RData")
stopCluster(cl)

})

##################################################################################################################
###### Make simulation results usable

res <- data.frame(do.call(rbind, simulation))
res <- res %>% mutate(mu = rep(mu, each = n.sim))

###### RMSE of the bootstrap variances

# suppose that the variance of the n.sim tau_hats is the "true" variance
# take RMSE relative to that.

res2 <- res %>% group_by(mu) %>%
  summarise(tau.hat.overall = mean(tau.hat),
            tau.hat.var.obs = var(tau.hat),
            tau.boot.overall = mean(tau.boot),
            tau.boot.var.obs = var(tau.boot),
            var.boot.overall = mean(var.boot),
            var.boot.var     = var(var.boot),
            var.boot.median  = median(var.boot),
            var.boot.rmse    = sqrt(mean((var.boot - var(tau.hat))^2))
  )
res2 <- res2 %>% mutate("Empirical Variance" = tau.hat.var.obs, "Bootstrap Variance" = var.boot.overall)

###### Define plotting theme

report_theme <- theme(
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 16),
  title = element_text(size = 20),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14)
)


###### Plot the bias
res2_melt <- melt(res2, id.vars = "mu")
png("../fig/KellieSimulationBias.png", width = 1000)
res2_melt %>% 
  filter(variable %in% c("Empirical Variance", "Bootstrap Variance")) %>%
  ggplot(aes(x = mu, y = value, color = variable)) + 
  geom_line() +
  xlab("mu") +
  ylab("Variance") +
  ggtitle("True and Bootstrap Variance") +
  report_theme +
  theme(legend.title = element_blank()) 
dev.off()

###### Plot the variance
res_melt <- melt(res, id.vars = "mu")

additional_points <- data.frame("mu" = mu, "variable" = "truth", "value" = res2$tau.hat.var.obs)
res_melt <- rbind(res_melt, additional_points)

png("../fig/KellieSimulationBias.png", width = 1000)
res_melt %>% filter(variable == "var.boot") %>%
  ggplot(aes(x = factor(mu), y = value)) +
  geom_boxplot() +
  geom_point(data = additional_points, aes(x = factor(mu), y = value), color = "red", size = 5) +
  xlab("mu") +
  ylab("Bootstrap Variance") +
  ggtitle("Distribution of the Bootstraped Variance across 1000 Simulations") +
  report_theme
dev.off()
