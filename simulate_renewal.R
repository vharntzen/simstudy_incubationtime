

library(dplyr)

setwd( "/Users/arntzenvh/Documents/Incubation and latent period Vietnam/incubation_and_latent_period/simulation_renewal")

# Load functions.
source("fun_simulate_data_renewal.R")
source("fun_estimate_renewal_new.R")



# Number of Monte Carlo replications
n.mc <- 1000



# Percentiles (vector of length six)
percentiles <- c(0.25, 0.5, 0.9, 0.95, 0.975, 0.99) # Change quantiles for Burr if you change this: does not run on server.



# Define scenarios.



### Data generation

n <- c(500, 1200)
pi <- c(0, 0.1, 0.2)
distribution_T <- c("Weibull", "heavytail",
                    "lognormal"
                    )
generation_method <- c("DengQin", "Leitzinger"
                       )



### Estimation
assumed_distribution <- c("lognormal", "gamma", "Weibull")
estimation_method <- c("pi", "nopi")



### Get all combinations
scenarios <- expand.grid(generation_method, n, pi, distribution_T, 
                         assumed_distribution, estimation_method)
colnames(scenarios) <- c("generation_method", "n", "pi", "distribution_T", 
                         "assumed_distribution", "estimation_method")
scenarios <- as.data.frame(scenarios)



### Remove those that make no sense.

#scenarios <- scenarios[which( (scenarios$n %in% c(500) & scenarios$generation_method == "Leitzinger") == FALSE) , ]
scenarios <- scenarios[which( (scenarios$generation_method == "Leitzinger" & scenarios$pi %in% c(0.1, 0.2)) == FALSE) , ]

### Add parameters.
# Parameterization based on appendix table 2 of Lauer (2020)
scenarios$par1 <- scenarios$par2 <- NA
scenarios$par1[which(scenarios$distribution_T == "lognormal")] <- 1.621 
scenarios$par2[which(scenarios$distribution_T == "lognormal")] <- 0.418
scenarios$par1[which(scenarios$distribution_T == "heavytail")] <- 8.5 # m
scenarios$par2[which(scenarios$distribution_T == "heavytail")] <- 2 # s, and f = 2
scenarios$par1[which(scenarios$distribution_T == "Weibull")] <- 2.453 
scenarios$par2[which(scenarios$distribution_T == "Weibull")] <- 6.258

scenarios_fordatagen <- scenarios %>% select(generation_method, n, pi, distribution_T, par1, par2) %>% unique()

# Select n = 1200 for data generation only, and randomly select 500 for each data set of 500.
scenarios_fordatagen <- scenarios_fordatagen %>% filter(n == 1200)

################################# Data generation ################################

for(run in 1:nrow(scenarios_fordatagen)){ # Start runs.

cat("Starting iteration", run, "of ", nrow(scenarios_fordatagen), ".\n")

# Select the settings.
n <- scenarios_fordatagen$n[run]
pi <- scenarios_fordatagen$pi[run]
generation_method <- scenarios_fordatagen$generation_method[run]
distribution_T <- scenarios_fordatagen$distribution_T[run]
par1 <- scenarios_fordatagen$par1[run]
par2 <- scenarios_fordatagen$par2[run]


# Generate data.
set.seed(123)
if(generation_method == "DengQin"){
mat_dat500 <- matrix(ncol = 1000, nrow = 500)
mat_dat <- replicate(n.mc, Simulate_data_DengQin(n = n,
                                               pi = pi,
                                               distribution = distribution_T,
                                               par1 = par1, par2 = par2,
                                               output = "forward") )

for(i in 1:n.mc){
mat_dat500[,i] <- mat_dat[,i][sample(1:length(mat_dat[,i]), 500, replace = FALSE)]
}

} else if(generation_method == "Leitzinger"){
mat_dat500 <- vector(mode = "list", length = n.mc)
mat_dat <- replicate(n.mc, Simulate_data_Leitzinger_new(daily_new_cases = 125, population_size = 11000000, leaveenter_per_day_prior19th = 245000,
                                                  leave_per_day_from19thonwards = 490000,
                                                  growth_rate = 0.14,
                                                  distribution = distribution_T,
                                                  par1 = par1, par2 = par2, output = "forward"), simplify = TRUE)
for(i in 1:n.mc){
mat_dat500[[i]] <- mat_dat[[i]][sample(1:length(mat_dat[[i]]), 500, replace = FALSE)]
}

}

print(head(mat_dat))

save(mat_dat, file = paste("data/", "n", n, "pi", pi, "generation_method", generation_method, "dist", distribution_T, "par1", par1, "par2", par2, ".RData", sep = ""))
mat_dat <- mat_dat500
save(mat_dat, file = paste("data/", "n", 500, "pi", pi, "generation_method", generation_method, "dist", distribution_T, "par1", par1, "par2", par2, ".RData", sep = ""))


} # End of runs.


#################################### Estimation ####################################

for(run in 1:nrow(scenarios)){ # Start of loop over runs.


cat("Running scenario", run, "of", nrow(scenarios), ".\n")
  
  
# Select the settings.
n <- scenarios$n[run]
pi <- scenarios$pi[run]
generation_method <- scenarios$generation_method[run]
distribution_T <- scenarios$distribution_T[run]
par1 <- scenarios$par1[run]
par2 <- scenarios$par2[run]
assumed_distribution <- scenarios$assumed_distribution[run]
estimation_method <- scenarios$estimation_method[run]

#If n = 500, still select data sets with 1200:
#n <- 1200

# Load data.
load(file = paste("data/", "n", n, "pi", pi, "generation_method", generation_method, "dist", distribution_T, 
                  "par1", par1, "par2", par2, ".RData", sep = "")) # Name is mat_dat

# Estimate.
# if(generation_method == "DengQin"){
#   
#   if(estimation_method == "nopi"){
#     out <- t( sapply(1:n.mc, function(i){Estimate_renewal_noPi(mat_dat[,i][sample(1:length(mat_dat[,i]), 500, replace = FALSE)], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
#   } else {
#     out <- t( sapply(1:n.mc, function(i){Estimate_renewal_Pi(mat_dat[,i][sample(1:length(mat_dat[,i]), 500, replace = FALSE)], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
#   }
#   
# } else if (generation_method == "Leitzinger"){
#   
#   # Note that mat_dat is a list because the sample size varies a bit: always around 1200, but may be a bit less or more.
#   
#   if(estimation_method == "nopi"){
#     out <- t( sapply(1:n.mc, function(i){Estimate_renewal_noPi(mat_dat[[i]][sample(1:length(mat_dat[[i]]), 500, replace = FALSE)], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
#   } else {
#     out <- t( sapply(1:n.mc, function(i){Estimate_renewal_Pi(mat_dat[[i]][sample(1:length(mat_dat[[i]]), 500, replace = FALSE)], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
#   }
#   
# }

# For n approximately 1200:

if(generation_method == "DengQin"){

        if(estimation_method == "nopi"){
out <- t( sapply(1:n.mc, function(i){Estimate_renewal_noPi(mat_dat[,i], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
        } else {
out <- t( sapply(1:n.mc, function(i){Estimate_renewal_Pi(mat_dat[,i], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
        }

} else if (generation_method == "Leitzinger"){

# Note that mat_dat is a list because the sample size varies a bit: always around 1200, but may be a bit less or more.

        if(estimation_method == "nopi"){
out <- t( sapply(1:n.mc, function(i){Estimate_renewal_noPi(mat_dat[[i]], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
        } else {
out <- t( sapply(1:n.mc, function(i){Estimate_renewal_Pi(mat_dat[[i]], assumed_distribution = assumed_distribution, percentiles = percentiles)}) )
        }

}

save(out, file = paste("results_per_run/", paste(colnames(scenarios), scenarios[run,], collapse = "_"), "_n500.RData", sep = ""))

print(head(out))

# Summarize.
out <- as.data.frame(out)

if(run == 1){
dat_results <- Summarize_results(results = out, scenario = scenarios[run,])
} else {dat_results <- rbind(dat_results, 
                               Summarize_results(results = out, scenario = scenarios[run,])
                             )
}

save(dat_results, file = paste("results/Results_after_run_", run, ".RData", sep = ""))

} # End of loop over runs.


