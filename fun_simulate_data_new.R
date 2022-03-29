
# _____________________________________________________________

# Function to generate single- or doubly interval-censored data
# using window length.

Simulate_data_windows <- function(n = 10, type = "doubly", window_lengths = 1:5, 
                          distribution_T = "Weibull", 
                          par_1 = 1, par_2 = 7, exposure_dist = "unif"){
  
  # Inputs
  
  # Output
  
  # Packages
  require(dplyr)
  require(tidyr)
 # require(poweRlaw)
  require(rmutil)  # Burr
  
  mat <- matrix(nrow = n, ncol = 7)
  colnames(mat) <- c("Li", "Ri", "Ti", "L0", "L1", "R0", "R1")
  
  for (i in 1:n){ # For every observation
    
    # Generate left window
    if(type == "single"){
    L0 <- 0
    } else if(type == "doubly"){
      
    # Daily incidence (derivative of y = a*exp(0.15t) ). 
    # plot(x = seq(0, 60, 0.01), y = sapply(X = seq(0, 60, 0.01), 
    #                                        FUN = function(t){0.15*exp(0.15*t)}), type = "l")
      
    # CalcIncidence <- function(t){0.15*1*exp(0.15*t)}
    # t <- seq(0, 30, 1)
    # probs <- CalcIncidence(seq(0, 30, 1))/sum(CalcIncidence(seq(0, 30, 1)))
  
    L0 <- runif(1, min = 0, max = 90) %>% round()
    }
  
    L1 <- L0 + sample(window_lengths, 1)  # interval L0 to L1 is between length 1 and 5.
    
    mat[i,4:5] <- c(L0, L1)
    
    # Generate L, R and T
    
    if(exposure_dist == "exp_growth"){
    #Li <- runif(1, min = L0, max = L1) #%>% round(0)
    # plot(x = seq(0, 60, 0.01), y = sapply(X = seq(0, 60, 0.01), 
    #                   FUN = function(t){exp(0.15*t)}), type = "l") # cumulative incidence
    a <- 1; g = 0.14  # a starting value; g growth factor
    #x <- runif(1, min = 0, max = L1 - L0)  # draw from uniform
    x <- runif(1, min = exp(0.14*(L0)), max = exp(0.14*(L1)))  # draw from uniform, 
    # for truncation: https://blogs.sas.com/content/iml/2013/07/22/the-inverse-cdf-method.html
    Li <- (1/g)*log(x/a) # inverse of cumulative incidence
    #y = -log(1-x)  # inverse exponential
    #Li <- L1 - y
   
    } else if(exposure_dist == "unif"){
    Li <- runif(1, min = L0, max = L1)  
    } else if(exposure_dist == "exp_household"){
    
    # p(1-p)^t
    p <- 0.2
    x <- runif(1, min = (p*(1-p)^L0) / log(1-p), max = (p*(1-p)^L1) / log(1-p)) # Truncation by including CDF
    Li <- log(x*log(1-p)/p) / log(1-p) # inverse cdf
    
    # Old
    #x <- runif(1, min = pexp(L0), max = pexp(L1)) # Truncation by including CDF
    #Li <- -log(1-x) # inverse cdf
   
     # ## Verify:
    #  y <- vector(length = 1000)
    #  for(i in 1:1000){
    #    x <- runif(1, min = pexp(L0), max = pexp(L1)) # Truncation by including CDF
    #  y[i] = -log(1-x)
    #  max(y)
    #  
    #  }
    # #
    # hist(y)
    }
    
    if(distribution_T == "Weibull"){
      Ti <- rweibull(1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "gamma"){
      Ti <- rgamma(1, shape = par_1, scale = par_2) #%>% round(0)
    } else if(distribution_T == "lognormal"){
      Ti <- rlnorm(1, meanlog = par_1, sdlog = par_2) #%>% round(0)
    } else if(distribution_T == "density"){
      Ti <- sample(x = values_T, size = 1, replace = T, prob = prob_T) #%>% round(0)
    } else if(distribution_T == "mixture"){
      choose <- rbinom(prob = 0.4, n = 1, size = 1)
      Ti <- ifelse(choose == 1, rgamma(1, shape = 0.5, scale = 5), rgamma(1, 0.1, 2)) #%>% round(0)
      
    } else if(distribution_T == "heavytail"){
    
    Ti <- rburr(n = 1, m = par_1, s = par_2, f = 2) # From rmutil
      
  
    #dlevy  
    # Ti <- rlevy(1000000, m=6, s=0.05)   #location parameter equal to m and dispersion equal to s (scale)
    # plot(x = seq(5.1,20,0.00001), dlevy(y = seq(5.1,20,0.00001), m=5, s=0.05), type = "l", xlim = c(4,14))
    # plot(x = seq(0.0001,30,0.01), dlevy(y = seq(0.0001,30,0.01), m=0, s=5), type = "l", xlim = c(0,10))
    # 
    # m <- 0:10; s <- seq(0.05, 0.2, 0.01)
    # tab <- expand.grid(m, s)
    # tab$median <- tab$ninetyfive <- NA
    # tab$median <- qlevy(0.5, m = tab$Var1, s = tab$Var2)
    # tab$ninetyfive <- qlevy(0.95, m = tab$Var1, s = tab$Var2)
    # tab
  
      
# qlevy(0.5, m = 5, s = 0.05)
# qlevy(0.95, m = 5, s = 0.05) #127
    
    # m = 6, s = 0.05 
    # 95% 18.71572  50% 6.1099055
    
    #https://en.wikipedia.org/wiki/L%C3%A9vy_distribution
      
    # Ti <- rep(NA, 1000)
    # for(i in 1:1000){
    #   
    # Ti[i] <- ifelse(rbinom(n = 1, size = 1, prob = 0.9) == 1, # This way, there is still a jump early on.
    #          rweibull(n = 1, shape = 2.453, scale = 6.258),
    #          rplcon(n = 1, xmin = 0, alpha = 1.8)
    #   )
    # }
    # 
    # 
    # Ti <- rcauchy(1000, location = 5, scale = 1)
    #   
    # plot(x = seq(-5,20,0.001), dcauchy(x = seq(-5,20,0.001), location = 0, scale = 1), 
    #      type = "l", xlim = c(-3,2))
    # 
    # library(rmutil)
    # plot(x = seq(0.001,20,0.001), dburr(y = seq(0.001,20,0.001), m = 8.5, s = 2, f = 2), 
    #      type = "l", xlim = c(0,20), ylim = c(0,0.2), col = "blue", xlab = "Incubation time T (days since infection)", 
    #      ylab = "P.D.F. f(T)", sub = "50% and 95% percentiles: Burr 5.47 and 15.84; Weibull 5.39 and 9.79")
    # lines(x = seq(0.001,20,0.001), dweibull(x = seq(0.001,20,0.001), shape = 2.453, scale = 6.258), col = "red")
    # legend("topright", legend = c("Weibull (shape 2.453, scale 6.258)", "Burr (m = 8.5, s = 2, f = 2)"), 
    #        col = c("red", "blue"), lty = c(1,1), bty = "n")
    # 
    # qburr(c(0.5, 0.95), m = 8.5, s = 2, f = 2) # m location (), s dispersion, f is family par (2)
    # qweibull(c(0.5, 0.95), shape = 2.453, scale = 6.258)
    # 
    
    
    
    # plot(x = seq(0.01,20,0.001), dburr(y = seq(0.01,20,0.001), m = 5, s = 3, f = 2), 
    #      type = "l", xlim = c(0,20))
    # 
    # qburr(0.5, m = 5, s = 3, f = 2)
    # qburr(0.95, m = 5, s = 1.5, f = 2)
    
    # Burr
    # Frechet
    #install.packages("evd")
    #library(evd)
    
    #dfrechet()
    
  #  x, loc=0, scale=1, shape=1, log = FALSE) 
#pfrechet(q, loc=0, scale=1, shape=1, lower.tail = TRUE) 
#qfrechet(p, loc=0, scale=1, shape=1, lower.tail = TRUE)
#rfrechet(n, loc=0, scale=1, shape=1)
    
      # Can be improved by saving the distribution. Not resampling.
      #install.packages("poweRlaw")
      #library("poweRlaw") # doi: 10.18637/jss.v000.i00
   #   xmin = 7.5; alpha = 1.8; x = seq(0,100,0.01)
  #    pdf_heavytailed <- dplcon(x, xmin, alpha) # xmin = 7.5, alpha = 1.8
  #    pdf_Weibull <- dweibull(x = seq(0,100,0.01), shape = 2.453, scale = 6.258)
      
      # plot(seq(0,100,0.01), type = "l", ylim = c(0,0.2), y = c(pdf_Weibull[1:750], 
      #   #  0.5*pdf_Weibull[72:401] + 0.5* 
      #   pdf_heavytailed[751:10001]), col = "blue")
      # lines(seq(0,100,0.01), pdf_heavytailed, lty = 3)
      # lines(seq(0,100,0.01), pdf_Weibull, lty = 2)
    
      # Cut-off at 100.
      # T <- sample(x = seq(0,100,0.01), size = 100000,  prob = c(pdf_Weibull[1:750],
      #                   pdf_heavytailed[751:10001]), replace = TRUE )
      # plot(density(T) )
      # quantile(T, c(0.5, 0.95) )
      
   #   Ti <- sample(x = seq(0,100,0.01), size = 1,  prob = c(pdf_Weibull[1:750],
    #                                                            pdf_heavytailed[751:10001]), replace = TRUE )
    }
    
    Ri <- Li + Ti
    mat[i,c(1:3)]  <- c(Li, Ri, Ti)
    
} # End of loop per individual.

# Generate/specify R0 and R1    
    if(type == "single"){
      # Make sure that R cannot be smaller than L1: set to R in that case.
      # This is a commonly made assumption, see Backer 2020, Tindale 2020.
      mat <- mat %>% as.data.frame()
      
     if(length(mat$L1[which(mat$Ri<mat$L1)])>0){
     mat$L1[which(mat$Ri<mat$L1)] <- mat$Ri[which(mat$Ri<mat$L1)]
     }
      
      mat[, 6:7] <-  mat[, 2] # Ri = R0 = R1
      
   } else if(type == "doubly"){
      
     for(i in 1:n){
    ## Old way:  
    #   # first_test <- min(60/L1*2, 30)  # Linear decrease in days between last exposure and first test.
    #   t <- seq(0, 30, 1)
    #   CalcIncidence <- function(t){0.15*1*exp(0.15*t)}
    #   probs <- CalcIncidence(seq(0, 30, 1))/sum(CalcIncidence(seq(0, 30, 1)))
    #   cutpoints <- c(0, sample(t, prob = probs, size = 10), max(mat[, "Ri"], na.rm = TRUE))
    #   first_test <- max(cutpoints[which(cutpoints <= Ri)])
    #   second_test <- min(cutpoints[which(cutpoints >= Ri)])
    #   
    #   # determine place by Ri
    #   R0 <- ifelse( rbinom(1, 1, prob = 0.2) == 1, # About 80% tests positive at first test.
    #         L0, # Set start of window infectiousness to start exposure.
    #         first_test  # Otherwise first test date.
    #   )
    #   R1 <- (second_test)  # Implicit assumption that windows get shorter in time.
      Ri <- mat[i, 2]  
       
     cutpoints <- c(0, runif(n = 15, min = 0, max = max(mat[, "Ri"]) ), max(mat[, "Ri"]) )
     R0 <- max(cutpoints[which(cutpoints <= Ri)])
     R1 <- min(cutpoints[which(cutpoints >= Ri)])
     
      # window_right <- sample(window_lengths, 1)
      # R0_to_R <- runif(1, min = 0, max = window_right)
      # R_to_R1 <- window_right - R0_to_R
      # 
      # R0 <- Ri - R0_to_R
      # R1 <- Ri + R_to_R1
      
      mat[i, 6:7] <- c(R0, R1)
      
      L0 <- mat[i, 4]
      L1 <- mat[i, 5]
      mat[i, 5] <- ifelse(R1 < L1, R1, L1)
      mat[i, 6] <- ifelse(R0 < L0, L0, R0)
    } }
    


# Problem: using the following statement observations with only zeros may occur.  
 # mat[, 4:7] <- floor(mat[, 4:7])
  
# New approach:
#  mat[, c(4,6)] <- floor(mat[, c(4,6)])
#  if(type == "doubly"){
#  mat[, c(5,7)] <- ceiling(mat[, c(5,7)])
#  } else{ mat[, 5] <- ceiling(mat[, 5])
#          mat[, 7] <- floor(mat[,7]) + 1 # + 0.5 to prevent 0s 
#          mat[, 6] <- mat[, 6] + 0.65} # R0 = R1
  
  return(mat %>% as.data.frame)
  
  # End of function   
  
}

#____________________________________________________________________________
# # Function to generate single- or doubly interval-censored data
# # using cutpoints
# 
# Simulate_data_cutpoints <- function(n = 10, type = "single", window_start = 5, 
#                           n_cutpoints = 3, distribution_T = "Weibull", 
#                           par_1 = 1, par_2 = 7, values_T = 3:9, 
#                           prob_T = c(0.0463850922, 0.2466837048, 0.0024858945, 
#                                      0.1126655228, 0.1347501680, 0.2058210187, 0.2512085991)){
#   
#   # Inputs
#   
#   # Output
#   
#   # Packages
#   require(dplyr)
#   require(tidyr)
#   
#   mat <- matrix(nrow = n, ncol = 7)
#   colnames(mat) <- c("Li", "Ri", "Ti", "L0", "L1", "R0", "R1")
#   
#   for (i in 1:n){
#     
#     # Generate L, R and T
#     Li <- runif(1, min = 0, max = window_start) %>%
#       floor()
#     
#     if(distribution_T == "Weibull"){
#       Ti <- rweibull(1, shape = par_1, scale = par_2) %>% round(0)
#     } else if(distribution_T == "gamma"){
#       Ti <- rgamma(1, shape = par_1, scale = par_2) %>% round(0)
#     } else if(distribution_T == "lognormal"){
#       Ti <- rlnorm(1, meanlog = par_1, sdlog = par_2) %>% round(0)
#     } else if(distribution_T == "density"){
#       Ti <- sample(x = values_T, size = 1, replace = T, prob = prob_T) %>% round(0)
#     } else if(distribution_T == "mixture"){
#       choose <- rbinom(prob = 0.4, n = 1, size = 1)
#       Ti <- ifelse(choose == 1, rgamma(1, shape = 0.5, scale = 5), rgamma(1, 0.1, 2)) %>% round(0)
#     }
#     
#     Ri <- Li + Ti
#     
#     mat[i,1:3]  <- c(Li, Ri, Ti)
#     
#   }
#   
#   # Define upper bounds for the interval to sample cut points from
#   Lmax <- max(mat[,"Li"])
#   Tmax <- max(mat[,"Ti"])
#   
#   for (i in 1:n){
#     
#     # Generate cut points C
#     Ci <- runif(n_cutpoints, 0, Lmax) %>%
#       round(0)
#     
#     Li <- mat[i, "Li"]
#     
#     Ci_and_Li <- c(Li, Ci) %>% sort() %>% unique()
#     
#     mat[i, 4] <- Ci_and_Li[
#       max(1, which(Ci_and_Li == Li) - 1) # L0
#     ]
#     
#     mat[i, 5] <- Ci_and_Li[
#       which(Ci_and_Li == Li) + 1 # L1
#     ] %>%  # This can give NA, hence:
#       replace_na(Lmax)
#     
#     if(type == "single"){
#       mat[i, 6:7] <-  mat[i, 2] # Ri = R0 = R1
#     }
#     
#     # Make sure that R cannot be smaller than L1: set to R in that case.
#     # This is a commonly made assumption, see Backer 2020, Tindale 2020.
#     mat <- mat %>% as.data.frame()
#     if(length(mat$L1[which(mat$R0<mat$L1)])>0){
#       mat$L1[which(mat$R0<mat$L1)] <- mat$R0[which(mat$R0<mat$L1)]
#     }
#     
#   }
#   
#   return(mat %>% as.data.frame())
#   
#   # End of function   
#   
# }

