


#### Function to simulate data as in Deng/Qin et al. (2020).
# Returns backward or forward time.
# -------------------------------------------------------------------------
Simulate_data_DengQin <- function(n, pi, distribution, par1, par2, output){
  
# Necessary for generating from a heavytail (Burr) distribution.
require(rmutil)

# Risk of getting infected while travelling is pi.
  pp <- 1 - pi
  
# Create container.
  if(output == "both"){
    obs_time = as.data.frame(matrix(NA, ncol = 2, nrow = n))
    colnames(obs_time) <- c("forward", "backward")
  } else{
    obs_time = rep(NA, n)
  }
  
# For every observation, do
    
for (i in 1:n){
  
      # Define status: travel on (A) or after (B)infection day.
      select = rbinom(1,1,pp)
      
      # A. Travelling at same day as infection.
      if (select==0){
        
              if(output == "both"){    
          if(distribution == "gamma"){
          obs_time[i, 1] = rgamma(1,par1,par2)
        }
        else if(distribution == "Weibull"){
          obs_time[i, 1] = rweibull(1,par1,par2) 
        }
        else if(distribution == "lognormal"){
          obs_time[i, 1] = exp(rnorm(1,par1,par2)) 
        } 
        else if(distribution == "heavytail"){
          obs_time[i, 1] = rburr(1, m = par1, s = par2, f = 2) 
        }
        
        obs_time[i, 2] <- 0 # Backward time is zero.
        
              } else {   
          if(distribution == "gamma"){
          obs_time[i] = rgamma(1,par1,par2)
        }
        else if(distribution == "Weibull"){
          obs_time[i] = rweibull(1,par1,par2) 
        }
        else if(distribution == "lognormal"){
          obs_time[i] = exp(rnorm(1,par1,par2)) 
        } 
        else if(distribution == "heavytail"){
          obs_time[i] = rburr(1, m = par1, s = par2, f = 2) 
        }
      
        }
  

      
      # B. Travelling after infection.
      } else{ # If select is non-zero.
      
        
        while (TRUE){
          C = runif(1,0,30)		#departure time
          if(distribution == "gamma"){
            Y = rgamma(1,par1,par2)
          }
          else if(distribution == "Weibull"){
            Y = rweibull(1,par1,par2) 
          }
          else if(distribution == "lognormal"){
            Y = exp(rnorm(1,par1,par2))
          }
          else if(distribution == "heavytail"){
            Y = rburr(1, m = par1, s = par2, f = 2) 
          }
          
          if (Y>C) break
        }
        
        # Collect output.
        if(output == "forward"){ 
          obs_time[i] <- Y-C #observed forward time
        }
        else if(output == "backward"){  
          obs_time[i] = Y-(Y-C)		#observed backward time
        }  else if(output == "both"){
          obs_time[i,1] <- Y-C #observed forward time
          obs_time[i,2] = Y-(Y-C)	#backward time
          
        }
      
    } # End of part B.
      
      
} # End of loop over n.
  
  return(round(obs_time,0))
    
}  

#Simulate_data_DengQin(1200, 0.2, "gamma", 5, 0.8, "forward")

  
# -----------------------------------------------------------------------------

#### Function to simulate data based on work by Leitzinger (master thesis, 2021)
# Returns backward or forward time.
# ------------------------------------------------------------------------------

Simulate_data_Leitzinger <- function(daily_new_cases, 
                                     leaveenter_per_day_prior19th, 
                                     leave_per_day_from19thonwards, 
                                     growth_rate, 
                                     distribution, 
                                     par1, par2, output){
  
  require(rmutil)
  
  infected_list <- list()
  incubation_periods <- list()
  infected_travel <- list()
  time_obs <- list()
  left_wuhan <- list()
  population <- c(1:11000000) # population at risk
  susceptibles <- c(1:11000000) # susceptibles
  #same as population_total at the start - to ensure that not > 1 infections per person
  
  day_of_inf <- vector(mode = "list", length = 19)
  
  ### Start the loop over the first 19 days.
  
      for (i in 1:19){
        
              # 1. Infection: draw a sample from people who are still in Wuhan and have not been infected before
            infected_list[[i]] <- sample(x = susceptibles, size = daily_new_cases, replace = FALSE)
            
              ### Vector with IDs of infecteds
            
          if (i < 15){
              
              # Incoming travelers: add to population (day 1 to 15 only)
              population <- c(population, (tail(population, 1) + 1):(tail(population, 1) + leaveenter_per_day_prior19th))
              
              ### Vector with population_at_risk: add IDs of incoming travelers and remove the exiting ones.
              ### Update the vector with susceptibles based on both.
              
              # Exiting travelers: extract from population
              left_wuhan[[i]] <- sample(x = population, size = leaveenter_per_day_prior19th, replace = FALSE)
              population <- population[!(population %in% left_wuhan[[i]])]
              susceptibles <- population[!(population %in% unlist(infected_list[1:i]))] 
              
          } else {
              
              # Exiting travelers: subset those that are infected - those are the ones of interest.
              left_wuhan[[i]] <- sample(x = population, size = leave_per_day_from19thonwards, replace = FALSE)
              population <- population[!(population %in% left_wuhan[[i]])]
              susceptibles <- population[!(population %in% unlist(infected_list[1:i]))] 
              infected_travel[[i]] <- intersect(left_wuhan[[i]], unlist(infected_list))
              
              # Store the time of infection for each infected traveler
              for (k in 1:length(infected_travel[[i]])){
                day_of_inf[[i]][k] <- which(sapply(X = infected_list, FUN = function(X) infected_travel[[i]][k] %in% X))[1]
              }
              # and generate incubation periods
              if(distribution == "gamma"){
                incubation_periods[[i]] <- rgamma(n = length(infected_travel[[i]]), shape = par1, rate = par2)
              } else if(distribution == "Weibull"){
                incubation_periods[[i]] <- rweibull(n = length(infected_travel[[i]]), shape = par1, scale = par2)
              } else if(distribution == "lognormal"){
                incubation_periods[[i]] <- rlnorm(n = length(infected_travel[[i]]), meanlog = par1, sdlog = par2)
              } else if(distribution == "heavytail"){
                incubation_periods[[i]] <- rburr(length(infected_travel[[i]]), m = par1, s = par2, f = 2) 
              }
              
              if(output == "forward"){ # Forward time
                # the forward times when symptom onset is not before day of travel are:
                time_obs[[i]] <- (day_of_inf[[i]][day_of_inf[[i]] + incubation_periods[[i]] > i] +
                                    incubation_periods[[i]][day_of_inf[[i]] + incubation_periods[[i]] > i]) - i
              }
              else if(output == "backward"){ # Backward time
                time_obs[[i]] <- i - (day_of_inf[[i]][day_of_inf[[i]] + incubation_periods[[i]] > i])
              }
            }
            
            # Incidence in number of cases per day
            daily_new_cases <- daily_new_cases + round(daily_new_cases * (exp(growth_rate) - 1)) # mistake
            
      } ### End the loop over the first 19 days.
          
  return( as.numeric(unlist(time_obs)) ) 
}
# -----------------------------------------------------------------------------


# To try
# daily_new_cases = 125
# leaveenter_per_day_prior19th = 245000
# leave_per_day_from19thonwards = 490000
# growth_rate = 0.14
# distribution = "gamma"
# par1 = 5; par2 = 0.8; output = "forward"
#Simulate_data_Leitzinger(daily_new_cases = 125, 245000, 490000, 0.14, "gamma", 5, 0.8, "forward")
# 
# 


#### Function to simulate data based on work by Leitzinger (master thesis, 2021)
# Returns backward or forward time.
# ------------------------------------------------------------------------------

Simulate_data_Leitzinger_new <- function(daily_new_cases, 
                                     population_size = 11000000,
                                     leaveenter_per_day_prior19th, 
                                     leave_per_day_from19thonwards, 
                                     growth_rate, 
                                     distribution, 
                                     par1, par2, output, print_daily_mat = FALSE){
  
  require(rmutil) # For Burr distribution

  # Create container.
  mat <- matrix(ncol = 3, nrow = population_size)
  
  # Column 1: infection
  # Column 2: travel from Wuhan
  # Column 3: symptom onset
  
  ### Start the loop over the first 19 days.
  
  for (i in 1:19){
    
    if(print_daily_mat == TRUE){cat("Matrix of day", i, "(col1 infection; col2 travel; col3 symptom):\n") 
                                print(head(mat) )}
    
    # 1. Infection: draw a sample from people who are still in Wuhan and have not been infected before
    
    ID_infecteds <- sample(x = which(is.na(mat[ , 1]) & is.na(mat[ , 2])), size = daily_new_cases, replace = FALSE) 
    # By doing is.na[mat,1] we only sample from susceptibles and by (is.na[ , 2]) from those did not leave Wuhan yet.
    
    mat[ID_infecteds, 1] <- i # Store the infection day.
    
    if (i < 15){
    
      # Incoming travelers: add to population (day 1 to 15 only)
      incoming_travelers <- matrix(ncol = 3, nrow = leaveenter_per_day_prior19th)
      mat <- rbind(mat, incoming_travelers) # NB: potentially slow
      
      # Exiting travelers: store travel date (and later on exclude these)
      ID_left_wuhan <- sample(x = which(is.na(mat[ , 2])), size = leaveenter_per_day_prior19th, replace = FALSE)
      # is.na(mat[ , 2]) to sample from those that did not travel from Wuhan (yet) only.
      mat[ID_left_wuhan, 2] <- i
      
    } else {
      
      # No incoming travelers.
      
      # Exiting travelers: sample from those still there (i.e. that did not leave Wuhan yet, so is.na(travel_day)).
      ID_left_wuhan <- sample(x = which(is.na(mat[ , 2])), size = leave_per_day_from19thonwards, replace = FALSE)
      mat[ID_left_wuhan, 2] <- i # Store the day of travel.
      
    }
  
  # Incidence in number of cases per day.
  daily_new_cases <- round(daily_new_cases*exp(growth_rate)) # This is the same as I_0 * exp(g*t)
  
  } # End of loop.
    
  ### mat <- mat[which(mat[,3] >= mat[,2]) , ]
  ### if(print_daily_mat) # NB TRUE hoeft niet
  ### 
  
      # Subset those that were infected AND traveled on the 15th, 16th, 17th, 18th or 19th.
  
      # Column 1: infection
      # Column 2: travel from Wuhan
      # Column 3: symptom onset
  
      mat <- mat[which(!is.na(mat[ , 1]) & mat[ , 2]>=15) , ] # mat[,1] non-missing for infecteds and mat[,2>=15] for selecting late travelers only
      
      # Generate incubation periods and define symptom onset date as infection day + incubation period.
      if(distribution == "gamma"){
        incubation_periods <- rgamma(n = nrow(mat), shape = par1, rate = par2)
      } else if(distribution == "Weibull"){
        incubation_periods <- rweibull(n = nrow(mat), shape = par1, scale = par2)
      } else if(distribution == "lognormal"){
        incubation_periods <- rlnorm(n = nrow(mat), meanlog = par1, sdlog = par2)
      } else if(distribution == "heavytail"){
        incubation_periods <- rburr(n = nrow(mat), m = par1, s = par2, f = 2) 
      }
      
      mat[, 3] <- mat[, 1] + incubation_periods
      
      # Subset those with symptom onset on or after travel.
      mat <- mat[which(mat[,3] >= mat[,2]) , ] # mat[,3] symptom onset day; mat[,2] travel day
      
      # Collect output.
      
        # Forward time
      if(output == "forward"){ 
        out <- mat[, 3] - mat[, 2] # Symptom onset day minus travel day
        # Backward time
      } else if(output == "backward"){ 
        out <- mat[, 2] - mat[, 1] # Travel day minus infection day
      } else if(output == "both"){
        
        out <- data.frame(forward = mat[, 3] - mat[, 2], backward = mat[, 2] - mat[, 1])
        
      }
      
    round(out, 0)
    
    }
    

    
  

# -----------------------------------------------------------------------------


#test <- Simulate_data_Leitzinger_new(daily_new_cases = 125, population_size = 11000000, 245000, 490000, 0.14, "gamma", 5, 0.8, "forward")
#length(test) 



#### Function to simulate data as in Deng/Qin et al. (2020).
# Returns backward or forward time. UNROUNDED
# -------------------------------------------------------------------------
Simulate_data_DengQin_unrounded <- function(n, pi, distribution, par1, par2, output){
  
  # Necessary for generating from a heavytail (Burr) distribution.
  require(rmutil)
  
  # Risk of getting infected while travelling is pi.
  pp <- 1 - pi
  
  # Create container.
  if(output == "both"){
    obs_time = as.data.frame(matrix(NA, ncol = 2, nrow = n))
    colnames(obs_time) <- c("forward", "backward")
  } else{
    obs_time = rep(NA, n)
  }
  
  # For every observation, do
  
  for (i in 1:n){
    
    # Define status: travel on (A) or after (B)infection day.
    select = rbinom(1,1,pp)
    
    # A. Travelling at same day as infection.
    if (select==0){
      
      if(output == "both"){    
        if(distribution == "gamma"){
          obs_time[i, 1] = rgamma(1,par1,par2)
        }
        else if(distribution == "Weibull"){
          obs_time[i, 1] = rweibull(1,par1,par2) 
        }
        else if(distribution == "lognormal"){
          obs_time[i, 1] = exp(rnorm(1,par1,par2)) 
        } 
        else if(distribution == "heavytail"){
          obs_time[i, 1] = rburr(1, m = par1, s = par2, f = 2) 
        }
        
        obs_time[i, 2] <- 0 # Backward time is zero.
        
      } else {   
        if(distribution == "gamma"){
          obs_time[i] = rgamma(1,par1,par2)
        }
        else if(distribution == "Weibull"){
          obs_time[i] = rweibull(1,par1,par2) 
        }
        else if(distribution == "lognormal"){
          obs_time[i] = exp(rnorm(1,par1,par2)) 
        } 
        else if(distribution == "heavytail"){
          obs_time[i] = rburr(1, m = par1, s = par2, f = 2) 
        }
        
      }
      
      
      
      # B. Travelling after infection.
    } else{ # If select is non-zero.
      
      
      while (TRUE){
        C = runif(1,0,30)		#departure time
        if(distribution == "gamma"){
          Y = rgamma(1,par1,par2)
        }
        else if(distribution == "Weibull"){
          Y = rweibull(1,par1,par2) 
        }
        else if(distribution == "lognormal"){
          Y = exp(rnorm(1,par1,par2))
        }
        else if(distribution == "heavytail"){
          Y = rburr(1, m = par1, s = par2, f = 2) 
        }
        
        if (Y>C) break
      }
      
      # Collect output.
      if(output == "forward"){ 
        obs_time[i] <- Y-C #observed forward time
      }
      else if(output == "backward"){  
        obs_time[i] = Y-(Y-C)		#observed backward time
      }  else if(output == "both"){
        obs_time[i,1] <- Y-C #observed forward time
        obs_time[i,2] = Y-(Y-C)	#backward time
        
      }
      
    } # End of part B.
    
    
  } # End of loop over n.
  
  return(obs_time)
  
}  

#Simulate_data_DengQin(1200, 0.2, "gamma", 5, 0.8, "forward")



#### Function to simulate data based on work by Leitzinger (master thesis, 2021)
# Returns backward or forward time. UNROUNDED
# ------------------------------------------------------------------------------

Simulate_data_Leitzinger_new_unrounded <- function(daily_new_cases, 
                                         population_size = 11000000,
                                         leaveenter_per_day_prior19th, 
                                         leave_per_day_from19thonwards, 
                                         growth_rate, 
                                         distribution, 
                                         par1, par2, output, print_daily_mat = FALSE){
  
  require(rmutil) # For Burr distribution
  
  # Create container.
  mat <- matrix(ncol = 3, nrow = population_size)
  
  # Column 1: infection
  # Column 2: travel from Wuhan
  # Column 3: symptom onset
  
  ### Start the loop over the first 19 days.
  
  for (i in 1:19){
    
    if(print_daily_mat == TRUE){cat("Matrix of day", i, "(col1 infection; col2 travel; col3 symptom):\n") 
      print(head(mat) )}
    
    # 1. Infection: draw a sample from people who are still in Wuhan and have not been infected before
    
    ID_infecteds <- sample(x = which(is.na(mat[ , 1]) & is.na(mat[ , 2])), size = daily_new_cases, replace = FALSE) 
    # By doing is.na[mat,1] we only sample from susceptibles and by (is.na[ , 2]) from those did not leave Wuhan yet.
    
    mat[ID_infecteds, 1] <- i # Store the infection day.
    
    if (i < 15){
      
      # Incoming travelers: add to population (day 1 to 15 only)
      incoming_travelers <- matrix(ncol = 3, nrow = leaveenter_per_day_prior19th)
      mat <- rbind(mat, incoming_travelers) # NB: potentially slow
      
      # Exiting travelers: store travel date (and later on exclude these)
      ID_left_wuhan <- sample(x = which(is.na(mat[ , 2])), size = leaveenter_per_day_prior19th, replace = FALSE)
      # is.na(mat[ , 2]) to sample from those that did not travel from Wuhan (yet) only.
      mat[ID_left_wuhan, 2] <- i
      
    } else {
      
      # No incoming travelers.
      
      # Exiting travelers: sample from those still there (i.e. that did not leave Wuhan yet, so is.na(travel_day)).
      ID_left_wuhan <- sample(x = which(is.na(mat[ , 2])), size = leave_per_day_from19thonwards, replace = FALSE)
      mat[ID_left_wuhan, 2] <- i # Store the day of travel.
      
    }
    
    # Incidence in number of cases per day.
    daily_new_cases <- round(daily_new_cases*exp(growth_rate)) # This is the same as I_0 * exp(g*t)
    
  } # End of loop.
  
  ### mat <- mat[which(mat[,3] >= mat[,2]) , ]
  ### if(print_daily_mat) # NB TRUE hoeft niet
  ### 
  
  # Subset those that were infected AND traveled on the 15th, 16th, 17th, 18th or 19th.
  
  # Column 1: infection
  # Column 2: travel from Wuhan
  # Column 3: symptom onset
  
  mat <- mat[which(!is.na(mat[ , 1]) & mat[ , 2]>=15) , ] # mat[,1] non-missing for infecteds and mat[,2>=15] for selecting late travelers only
  
  # Generate incubation periods and define symptom onset date as infection day + incubation period.
  if(distribution == "gamma"){
    incubation_periods <- rgamma(n = nrow(mat), shape = par1, rate = par2)
  } else if(distribution == "Weibull"){
    incubation_periods <- rweibull(n = nrow(mat), shape = par1, scale = par2)
  } else if(distribution == "lognormal"){
    incubation_periods <- rlnorm(n = nrow(mat), meanlog = par1, sdlog = par2)
  } else if(distribution == "heavytail"){
    incubation_periods <- rburr(n = nrow(mat), m = par1, s = par2, f = 2) 
  }
  
  mat[, 3] <- mat[, 1] + incubation_periods
  
  # Subset those with symptom onset on or after travel.
  mat <- mat[which(mat[,3] >= mat[,2]) , ] # mat[,3] symptom onset day; mat[,2] travel day
  
  # Collect output.
  
  # Forward time
  if(output == "forward"){ 
    out <- mat[, 3] - mat[, 2] # Symptom onset day minus travel day
    # Backward time
  } else if(output == "backward"){ 
    out <- mat[, 2] - mat[, 1] # Travel day minus infection day
  } else if(output == "both"){
    
    out <- data.frame(forward = mat[, 3] - mat[, 2], backward = mat[, 2] - mat[, 1])
    
  }
  
out
  
}





# -----------------------------------------------------------------------------


#test <- Simulate_data_Leitzinger_new(daily_new_cases = 125, population_size = 11000000, 245000, 490000, 0.14, "gamma", 5, 0.8, "forward")
#length(test) 



