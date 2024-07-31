# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(rjags)
library(rstan)
library(bayesplot)
library(coda)


# Read the data
library(readr)
data_MZ_Oncho_iu <- read_csv("data-MZ-Oncho-iu.csv")
data_MZ_Oncho_iu = data.frame(data_MZ_Oncho_iu)
head(data_MZ_Oncho_iu)

#data_MZ_LF_iu$PopTot = as.numeric(data_MZ_LF_iu$PopTot)
#data_MZ_LF_iu$PopTreat = as.numeric(data_MZ_LF_iu$PopTreat)

#summary(data_MZ_LF_iu$PopTot)

#data_MZ_SCH_iu <- read_csv("C:/Users/User/Desktop/Debashis 2024/Tropical Disease PLOS/Datasets/data-MZ-SCH-iu.csv")
#data_MZ_STH_iu <- read_csv("C:/Users/User/Desktop/Debashis 2024/Tropical Disease PLOS/Datasets/data-MZ-STH-iu.csv")


data <- data_MZ_Oncho_iu
#data <- data_MZ_SCH_iu
#data <- data_MZ_STH_iu


# Display the first few rows of the data
head(data)
summary(data)
# Check the column names
colnames(data)

# Prepare data by filtering out rows with NA values in PopTot and PopTreat
data <- data %>%
  filter(!is.na(PopTot) & !is.na(PopTreat) & PopTot > 0 & PopTreat >= 0 & PopTreat <= PopTot) %>%
  mutate(prevalence = PopTreat / PopTot)

# Display summary statistics
summary(data)
sum(is.na(data))

#data$log_PopTot = log(data$PopTot)
#data$log_Trear= log(data$PopTreat)

# JAGS model with Jeffrey's prior
jeffreys_model <- "
model {
  for (i in 1:N) {
    cases[i] ~ dbin(theta[i], population[i])
    theta[i] ~ dbeta(0.5, 0.5)
  }
}
"

summary(data$PopTreat)
#plot(density(log(data$PopTreat)))

#log(data$PopTreat)
# Prepare data for JAGS
jags_data <- list(
  cases = data$PopTreat,
  population = data$PopTot,
  N = nrow(data)
)

# Parameters to monitor
parameters <- c("theta")

# Initial values
inits <- function() list(theta = runif(nrow(data)))

# Run the JAGS model
jags_fit <- jags.model(
  textConnection(jeffreys_model),
  data = jags_data,
  inits = inits,
  n.chains = 3,
  n.adapt = 1000
)

# Update and sample
update(jags_fit, 1000)
samples <- coda.samples(jags_fit, parameters, n.iter = 5000)

# Combine chains and convert to dataframe
combined_samples <- as.mcmc(do.call(rbind, samples))
samples_df <- as.data.frame(combined_samples)


# Summary of the samples
jags_summary <- summary(samples)
print(jags_summary)
# Extract parameter names and print them
parameter_names <- varnames(samples)
print(parameter_names)
ss_orch = samples[,c(1:4)]
# Select the first few parameters for plotting
selected_parameters <- parameter_names[1:4]

# Plot posterior distributions for the first few parameters
jpeg("jags_posteriors.jpg")
mcmc_areas(samples, pars = selected_parameters, prob = 0.95) +
  ggtitle("Posterior distributions with 95% credible intervals for first few parameters")
  
dev.off()

# Trace plots for the first few parameters
jpeg("jags_trace.jpg")
mcmc_trace(samples, pars = selected_parameters) +
  ggtitle("Trace plots for first few parameters")
dev.off()

head(as.data.frame(jags_summary$statistics))


# Diagnostics for JAGS model
jpeg("jags_gelman.jpg")
gelman_diag <- gelman.diag(ss_orch)
print(gelman_diag)
dev.off()

gelman.plot(ss_orch)

# Save Gelman-Rubin diagnostics for JAGS to a CSV file
write.csv(as.data.frame(gelman_diag$psrf), "jags_gelman.csv")

# Save effective sample size for JAGS and Stan
jpeg("jags_ess.jpg")
jags_ess <- effectiveSize(samples)
print(jags_ess)
dev.off()
