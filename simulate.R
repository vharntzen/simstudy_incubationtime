#######################################
### R script to run on ALICE server ###
#######################################

# Loading libraries (make sure everything is installed)
# For installation: https://wiki.alice.universiteitleiden.nl/index.php?title=R_on_ALICE#Loading_and_Installing_R_packages
library(doParallel, quiet = TRUE)  # installed
library(foreach) # installed
library(dplyr) # installed
library(tidyr) # installed
library(doSNOW)
library(boot)
# For function Reich:
### coarseDataTools
### dplyr # installed
### magrittr # installed
# For data simulation:
### dplyr # installed
### tidyr # installed

##### ALICE specific vs. try-out locally #####

# Needed folders: results_run, results_each

# Unhash for try-out:
setwd("/Users/arntzenvh/Documents/Incubation and latent period Vietnam/incubation_and_latent_period/simulation")
#setwd("/Users/arntzenvh/To_ALICE/Run_files")
cores <- 3 # Has 4, uses 3

# For try-out:
max_cpu <- 10000 
#max_cpu <- 100

# Hash for try-out:
# Get the number of cores
#cores <- Sys.getenv(paste("SLURM_CPUS_PER_TASK"))
#jobname <- Sys.getenv(paste("SLURM_JOB_ID"))

####################

# Loading functions
source("fun_estimate_single_interval_censored.R")
#source("fun_estimate_doubly_interval_censored.R")
source("fun_simulate_data_new.R")

load("exposure_window_lengths_based_on_real_data.RData") # lengths

# Get the starting time of the script
start_time <- proc.time()

###################################
##### SETTINGS ####################
###################################

# VARY: CI METHODS FOR MLE.
#       

##### Simulated data ##############
n = c(
100#, 
#500
  )

type = c("single"#,
      #  "doubly"
         )

distribution_T = c(#"lognormal", 
 "Weibull"#,
# "heavytail"
  ) 

windowlength <- c("standard"#, 
 #"double", "squared"
                  )

exposure_dist = c("unif"#, "exp_growth", 
#"exp_household"
  )

##### Method #####################
approach = c(#"Reich x lognormal", 
#"Reich x Weibull", "Reich x gamma", 
#  "Groeneboom"#,
#"MLE_surv x lognormal",  
#"MLE_surv x Weibull",
"MLE_flexsurv x lognormal"#,  
#"MLE_flexsurv x Weibull"#,
#"MLE_flexsurv x gamma",
#"NPMLE", #"real_quantiles"
#"PGM"
  )

# https://cran.r-project.org/web/packages/flexsurv/flexsurv.pdf
# The generalized gamma distribution simplifies to the gamma, log-normal and Weibull distributions
# with the following parameterisations:
# p.40
# dgengamma(x, mu, sigma, Q=sigma) = dgamma(x, shape=1/sigma^2, rate=exp(-mu) / sigma^2)

##### Monte Carlo simulations #####
N <- 1000

##### Percentiles #####
# NB: needs to be length 5
#percentiles <- c(0.025, 0.05, 0.5, 0.95, 0.975)
percentiles <- c(0.5, 0.9, 0.95, 0.975, 0.99)
row_names_percentiles <- c("q50", "q90", "q95", "q97.5", "q99")

###################################
###################################
###################################

# Combine all options and select rows.
# Make a loop from here to do multiple scenarios in once.

scenarios <- expand_grid(n, type, exposure_dist, distribution_T, approach, N, windowlength)
##### SUBSET SCENARIOS
#scenarios <- scenarios[19,]  # For try-out
n_scenarios <- nrow(scenarios)

save(scenarios, file = "scenarios.RData")

# Build in a check for availability of all data sets for the scenarios.

# ##### Generate data #######
# 
# 
#  library(rmutil) # N.B. Server cannot run this because of different version of R.
# 
#  for (run in 1:n_scenarios){
# 
# set.seed(123)
# 
# # Select settings
# n <- scenarios$n[run]
# type <- scenarios$type[run]
# exposure_dist <- scenarios$exposure_dist[run]
# distribution_T <- scenarios$distribution_T[run]
# windowlength <- scenarios$windowlength[run]
# 
# # Select pool of windows to sample from
# if(windowlength == "standard"){
#   window_lengths = lengths # sampled with equal probability
# } else if(windowlength == "double"){
#   window_lengths = 2*lengths
# } else if(windowlength == "squared"){
#   window_lengths = lengths^2
# }
# 
# # Parameterization based on appendix table 2 of Lauer (2020)
# if(distribution_T == "lognormal"){
#   par1_distribution_T = 1.621
#   par2_distribution_T = 0.418
# } else if(distribution_T == "Weibull"){
#   par1_distribution_T = 2.453
#   par2_distribution_T = 6.258
# } else if(distribution_T == "gamma"){
#   par1_distribution_T = 5.807
#   par2_distribution_T = 0.948
# } else if(distribution_T == "heavytail"){
#   par1_distribution_T = 8.5 # m
#   par2_distribution_T = 2 # s
# }
# 
# list_dat <- replicate(n = N, Simulate_data_windows(n = n, window_lengths = window_lengths, type = type,
#                              distribution_T = distribution_T,
#                              par_1 = par1_distribution_T,
#                              par_2 = par2_distribution_T,
#                              exposure_dist = exposure_dist), simplify = FALSE) # unif or exp_growth
# 
# approach_nr <- which(colnames(scenarios) == "approach")
# 
# filename_temp = paste(colnames(scenarios[run, -approach_nr]), c(scenarios[run, -approach_nr], ".RData"), sep = "_")
# 
# # Save data sets
# save(list_dat, file = paste("data/", paste(filename_temp, collapse = "_"), sep = ""))
# 
# cat("Starting data generation of scenario ", run, ".\n")
# print(scenarios[run,])
# print(head(list_dat[[1]]))
# 
# }

##### Function #####
# The function that does the actual work
myProc <- function(iteration){

# Load the data
approach_nr <- which(colnames(scenarios) == "approach")
filename_temp = paste(colnames(scenarios[run, -approach_nr]), c(scenarios[run, -approach_nr], ".RData"), sep = "_")
load(file = paste("data/", paste(filename_temp, collapse = "_"), sep = ""))
dat <- list_dat[[iteration]]
    
# Estimate interval

#try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf)
#{
#  try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE) 
#}

out <- try({setTimeLimit(cpu = max_cpu, elapsed = Inf);
  if(method == "Reich" & type == "single"){  
    out <- Estimate_single_Reich(dat, distribution = assumed_dist_T) %>% as.data.frame()
  }else if(method == "Groeneboom" & type == "single"){
    out <- Estimate_single_Groeneboom(dat) %>% as.data.frame()
  }else if(method == "Reich" & type == "doubly"){
    out <- Estimate_doubly_Reich(dat, distribution = assumed_dist_T) %>% as.data.frame()
  }else if(method == "Groeneboom" & type == "doubly"){
    out <- Estimate_doubly_Groeneboom(dat) %>% as.data.frame()
  }else if(method == "NPMLE"){
   out <- Estimate_single_npmle_survival(dat, conf.int = F) %>% as.data.frame()
    #out <- Estimate_single_npmle_survival(dat, conf.int = "survfit") %>% as.data.frame()
   #out <- Estimate_single_npmle_survival(dat, conf.int = "bootstrap") %>% as.data.frame()
  }else if(method == "PGM"){
    #out <- Estimate_single_npmle_survival(dat, conf.int = F) %>% as.data.frame()
   out <- Estimate_single_PGM(dat, conf.int = "boot", n.boots = N) %>% as.data.frame()
#    list_out <- Estimate_single_PGM(dat, conf.int = FALSE)
#    out <- list_out[[1]] %>% as.data.frame()
#    save(list_out, file = paste("results_par/Scen_", run,  "_Iteration_", iteration, ".RData", sep = ""))
  }else if(method == "MLE_surv"){
    out <- Estimate_single_mle_survival(dat, distribution = assumed_dist_T) %>% as.data.frame() # Not available for gamma.
    #list_out <- Estimate_single_mle_survival(dat, distribution = assumed_dist_T)
    #out <- list_out[[1]] %>% as.data.frame()
    #save(list_out, file = paste("results_par/Scen_", run,  "_Iteration_", iteration, ".RData", sep = ""))
  }else if(method == "MLE_flexsurv"){
    out <- Estimate_single_mle_flexsurv(dat, distribution = assumed_dist_T#, 
                                        #init.vals = c(par1_distribution_T, par2_distribution_T) 
                                        ) %>% as.data.frame() 
  }else if(method == "real_quantiles"){
    out <- data.frame(est = quantile(dat$Ti, probs = percentiles)) %>% as.data.frame()
  }
}#, silent = TRUE
)  # Try() function for error handling.
  
# Discriminate between error output and valid output
if( (class(out) == "try-error") == TRUE){
  cat(paste(out, "\n", sep = ""))
  estimates <- data.frame(est = rep(NA, length(percentiles)) )
  
  save(dat, file = paste("Data_of_scenario_", run, "_that_timed_out_or_gives_error.RData", sep = ""))
  
  # Indicate if CI present
  CI_present <- "no"
  
} else {
  estimates <- out # data.frame with column est.
  
  # Indicate if CI present
  CI_present <- ifelse( ("lower_CI" %in% colnames(out)) == T, "yes", "no")
  
}

  # Determine real values
  
  if(distribution_T == "gamma"){
    real_p <- qgamma(percentiles, shape = par1_distribution_T, scale = par2_distribution_T)
  }else if(distribution_T == "Weibull"){
    real_p <- qweibull(percentiles, shape = par1_distribution_T, scale = par2_distribution_T)
  }else if(distribution_T == "lognormal"){
    real_p <- qlnorm(percentiles, meanlog = par1_distribution_T, sdlog = par2_distribution_T)
  }else if(distribution_T == "heavytail"){
    #library(rmutil)
    #par1_distribution_T = 8.5 # m
    #par2_distribution_T = 2 # s
    real_p <- c(5.470551, 12.498982, 15.838618, 19.613748, 25.500000)
      #qburr(percentiles, m = par1_distribution_T, s = par2_distribution_T, f = 2) # rmutil
  }
  
  
  # Combine
  if(CI_present == "no"){
  ests <- data.frame(real_p = real_p, estimates = estimates$est)
  } else {ests <- data.frame(real_p = real_p, estimates = estimates$est, 
                             lower_CI = estimates$lower_CI, upper_CI = estimates$upper_CI)}
  
  # Evaluate
  
    if(CI_present == "yes"){
      results_tab <- ests %>% mutate( # Add necessary compounds to calculate performance measures later on.
        coverage = ifelse(ests$real_p > ests$lower_CI & ests$real_p < ests$upper_CI, 1, 0),  
        # If true parameter value is in estimated CI 1, otherwise 2.
        bias = ests$estimates - ests$real_p,
        CI_length = ests$upper_CI - ests$lower_CI
        # Difference between estimate and true parameter value.
        
      ) %>% select(coverage, 
                   bias, CI_length)
    } else {
      results_tab <- ests %>% mutate( # Add necessary compounds to calculate performance measures later on.
      coverage = NA,
      bias = ests$estimates - ests$real_p,
      CI_length = NA
      # Difference between estimate and true parameter value.
      
    ) %>% select(coverage, 
                 bias, CI_length)
    }
    
    save(results_tab, file = paste("results_each/RESULT_Scen_", run,
                                   "_Iter_", iteration, "_",
                                   paste(
                                     paste(colnames(scenarios[run, ]), c(scenarios[run, ], ".RData"), sep = "_"),
                                     collapse = "_"),
                                   sep = ""))

    results_tab
    
} # End of function that does actual work

###################### Start loop for multiple scenarios ######################
for(run in 1:n_scenarios){
#for(run in c(1,4)){
  
  # Select settings
  n <- scenarios$n[run]
  type <- scenarios$type[run]
  exposure_dist <- scenarios$exposure_dist[run]
  distribution_T <- scenarios$distribution_T[run]
  approach <- scenarios$approach[run]
  
  # Parameterization based on appendix table 2 of Lauer (2020)
  if(distribution_T == "lognormal"){
  par1_distribution_T = 1.621 
  par2_distribution_T = 0.418
  } else if(distribution_T == "Weibull"){
  par1_distribution_T = 2.453 
  par2_distribution_T = 6.258
  } else if(distribution_T == "gamma"){
  par1_distribution_T = 5.807
  par2_distribution_T = 0.948
  } else if(distribution_T == "heavytail"){ # Burr distribution
  par1_distribution_T = 8.5 # m
  par2_distribution_T = 2 # s
    # f = 2
  }
  
  if(approach == "Reich x lognormal"){
    method = "Reich"; assumed_dist_T = "lognormal"
  }else if(approach == "Reich x Weibull"){
    method = "Reich"; assumed_dist_T = "Weibull"
  }else if(approach == "Reich x gamma"){
    method = "Reich"; assumed_dist_T = "gamma"
  }else if(approach == "Groeneboom"){
    method = "Groeneboom"; assumed_dist_T = "none"
  }else if(approach == "MLE_surv x lognormal"){
    method = "MLE_surv"; assumed_dist_T = "lognormal"
  }else if(approach == "MLE_surv x Weibull"){
    method = "MLE_surv"; assumed_dist_T = "Weibull"
  }else if(approach == "MLE_flexsurv x lognormal"){
    method = "MLE_flexsurv"; assumed_dist_T = "lognormal"
  }else if(approach == "MLE_flexsurv x Weibull"){
    method = "MLE_flexsurv"; assumed_dist_T = "Weibull"
  }else if(approach == "MLE_flexsurv x gamma"){
    method = "MLE_flexsurv"; assumed_dist_T = "gamma"
  }else if(approach == "real_quantiles"){
    method = "real_quantiles"#; assumed_dist_T = "none"
  }else if(approach == "NPMLE"){
    method = "NPMLE"#; assumed_dist_T = "none"
  }else if(approach == "PGM"){
    method = "PGM"#; assumed_dist_T = "none"
  }

# Go through the simulation runs in parallel
#set.seed(123) # Make sure that generated data is always the same (but different per run)

cat("Now starting scenario ", run, ". \n Errors occurring: \n", sep = "")

# Initiate compute environment for doparallel
cl <- makeCluster(as.numeric(cores)-1, type = "SOCK")
#registerDoParallel(cl)

#pb <- txtProgressBar() # utils

#https://stackoverflow.com/questions/5423760/how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r

registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

result <- foreach(i = 1:N,
                  .packages = c("magrittr", #"coarseDataTools",
                               "tidyr", "dplyr"
                                ), # Solves the bug: https://stackoverflow.com/questions/41117127/r-parallel-programming-error-in-task-1-failed-could-not-find-function
                  #.combine = 'list',
                  .multicombine = TRUE, # Returns a list of data.frames using default (list).
                  .inorder = TRUE,
                  .options.snow = opts
                  ) %dopar% { 
       try(myProc(i), silent = TRUE)
                  }
close(pb)

# Remove compute environment
invisible(stopCluster(cl)) 



# Error in if (update.init && fitA$fail == 0 && fitA$degree.smooth$df >  : 
# missing value where TRUE/FALSE needed
# Error in t.default(ccoef) : 
# 'rho' must be an environment not pairlist: detected in C-level applyClosure

# tab <- vector("list", N)
# for(i in 933:N){
#               tab[[i]] <- try(myProc(i))
#               print(i)
#               }
# tab
# doSnow .options


if( (class(result) == "try-error") == TRUE){result <- replicate(
                                                                 N, data.frame(coverage = rep(NA, length(percentiles)), 
                                                                 bias = rep(NA, length(percentiles))), simplify = F
                                                                  )
  }

# For Reich, MLEs can sometimes not be estimated for all bootstrap replications such that CI is based on less.

# Here, SAVE result as a list/table
save(result, file = paste("results_each/List_Scen_", run, 
                               paste(
                                 paste(colnames(scenarios[run, ]), c(scenarios[run, ], ".RData"), sep = "_"), 
                                 collapse = "_"), 
                               sep = ""))

# Assemble summary
coverage <- bias <- bias_q25 <- bias_q75 <- mse <- mab <- mean_CI_length <- CI_q25 <- CI_q75 <- invalid <- as.numeric( rep(NA, length(percentiles)) )

for(j in 1:length(percentiles)){
  coverage[j] <- mean( sapply(result, FUN = function(tab){tab[j,1]}), na.rm = TRUE ) # coverage
  bias[j] <- mean( sapply(result, FUN = function(tab){tab[j,2]}) , na.rm = TRUE) # average bias
  
  if(N == 1000){ # Variation in the bias across runs.
  bias_q25[j] <- sort( sapply(result, FUN = function(tab){tab[j,2]}))[250] # average bias
  bias_q75[j] <- sort( sapply(result, FUN = function(tab){tab[j,2]}))[750] # average bias
  }
  
  mse[j] <- mean( sapply(result, FUN = function(tab){tab[j,2]})^2 , na.rm = TRUE) # MSE
  mab[j] <- mean( sapply(result, FUN = function(tab){abs(tab[j,2])}) , na.rm = TRUE) # MAB
  mean_CI_length[j] <- mean( sapply(result, FUN = function(tab){abs(tab[j,3])}) , na.rm = TRUE) # mean_CI_length
  
  if(N == 1000){ # Variation in the CI length across runs.
    CI_q25[j] <- sort( sapply(result, FUN = function(tab){tab[j,3]}))[250]
    CI_q75[j] <- sort( sapply(result, FUN = function(tab){tab[j,3]}))[750] 
  }
  
  invalid[j] <- 100 * sum(is.na(sapply(result, FUN = function(tab){abs(tab[j,2])})))/N
  }

result_tab <- data.frame(coverage = coverage, bias = bias, mse = mse, 
                         mab = mab, bias_q25 = bias_q25, bias_q75 = bias_q75, 
                         mean_CI_length = mean_CI_length, CI_q25 = CI_q25, CI_q75 = CI_q75, 
                         invalid = invalid) %>% round(3)

rownames(result_tab) <- paste("q", 100*percentiles, sep = "")

result_tab

print(scenarios[run, ])
print(result_tab)

# Save summary as .RData object
#save(result_tab, file = paste("Scenario_nr_", run, ".RData", sep = ""))
save(result_tab, file = paste("results_run/", paste(
  paste(colnames(scenarios[run, ]), c(scenarios[run, ], ".RData"), sep = "_"),
  collapse = "_"), sep = ""))


} ### End the loop over n scenarios.

# Summary table
results <- data.frame()

for(run in 1:nrow(scenarios)){
load(file = paste("results_run/", paste(
  paste(colnames(scenarios[run, ]), c(scenarios[run, ], ".RData"), sep = "_"),
  collapse = "_"), sep = ""))

temp <- cbind(percentiles, scenarios[run, ], result_tab)
results <- rbind(results, temp)
}

save(results, file = paste("Results_run_", date(), ".RData", sep = ""))

#copy_to_clipboard(results)


# Get the running time of the script and print it
running_time <- proc.time() - start_time
cat(paste("Total running time was ", running_time, ".\n", sep = ""))


