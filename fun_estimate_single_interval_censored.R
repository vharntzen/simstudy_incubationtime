
# Estimate_single_Groeneboom
# Estimate_single_Reich
# Estimate_single_npmle_survival
# Estimate_single_mle_survival

percentiles <- c(0.5, 0.9, 0.95, 0.975, 0.99)
row_names_percentiles <- c("q50", "q90", "q95", "q97.5", "q99")


#percentiles <- c(0.025, 0.050, 0.500, 0.950, 0.975)

# Penalized Gaussian mixture

# Function to estimate incubation time from 
# single interval censored 
# observations 
# using penalized Gaussian mixture.
Estimate_single_PGM <- function(dat, pdf_cdf = FALSE, conf.int = FALSE, n.boots = 1000, lambda = NA){  
  
  # Inputs:
  # dat     data frame with columns L0, L1, R0, R1
  
  require('smoothSurv')
  require('tidyr')
  require('magrittr')
  require('dplyr')
  require('survival')
  require('boot')
  
  dat$shortest <- dat$R0 - dat$L1
  dat$longest <- dat$R0 - dat$L0
  
  # Each observation as a time interval with (-infinity, t) for left censored, (t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
  dat$shortest[dat$shortest == 0] <- dat$shortest[dat$shortest == 0] + 0.000001 # To make the function work.
  
  #dat_new$t.left <- (dat_new$R0 - dat_new$L1) #+ 0.5
  #dat_new$t.right <- (dat_new$R0 - dat_new$L0) #+ 0.5
  
  # From: https://stats.stackexchange.com/questions/176376/invalid-survival-times-for-this-distribution
  #The survreg function in R does not allow time = 0. This is because for several of the distributions, including the lognormal distribution, having events occur at time = 0 will result in an undefined estimator.
 
invisible(capture.output(
  
  if(!is.na(lambda)){
    
  fit <- try(smoothSurvReg(Surv(shortest, longest, type = "interval2") ~ 1, aic = TRUE, lambda = lambda,
                       # the optimal value of λ is found using a grid search, based on AIC.
                       update.init = TRUE, init.dist = "best",
                       data = dat))
  } else {
  fit <- try(smoothSurvReg(Surv(shortest, longest, type = "interval2") ~ 1, aic = TRUE,
                         # the optimal value of λ is found using a grid search, based on AIC.
                         update.init = TRUE, init.dist = "best",
                         data = dat)
    )
  }
                     #  smoothSurvReg.control(difforder = 2))  # m (the order of the finite difference used in
                                    # the penalty term in smoothSurvReg.control) is 3 by default
  ))
  
  # t.right may be NA for right interval censored
  
  # Reading:
  # p.100 of Survival analysis with interval-censored data: a practical approach with examples in R (...)
  
  invisible(capture.output(
  S <- survfit(fit, by = 0.0001, xlim = c(0, 30))
  ))
  
  pdf <- fdensity(fit, xlim = c(0, 30))
  colnames(pdf) <- c("t", "f(t)")
  
  estimates <- sapply(percentiles, FUN = function(p){S$x[which(S$y1 >= (1 - p)) %>% tail(1)]}) # Percentiles
  
  colnames(S) <- c("t", "S(t)")
  
  if(conf.int == FALSE){
    tab <- data.frame(est = estimates)
    
    rownames(tab) <- row_names_percentiles
    
    
  } else if(conf.int %in% c("boot", TRUE)){ # Using package 'boot'
    
  
   
    tab <- data.frame(est = estimates)
    
    rownames(tab) <- row_names_percentiles
    
    ###### Function for bootstrap #####
    fun_boot  <- function(data, indices){
      
      boot_dat <- data[indices, ]
      
        # Model fit
        boot_dat$shortest <- boot_dat$R0 - boot_dat$L1
        boot_dat$longest <- boot_dat$R0 - boot_dat$L0
        
        # Each observation as a time interval with (-infinity, t) for left censored, (t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
        boot_dat$shortest[boot_dat$shortest == 0] <- boot_dat$shortest[boot_dat$shortest == 0] + 0.000001 # To make the function work.
        
        invisible(capture.output(
          if(!is.na(lambda[1])){
           
            fit <- try(smoothSurvReg(Surv(shortest, longest, type = "interval2") ~ 1, aic = TRUE, lambda = lambda,
                                 # the optimal value of λ is found using a grid search, based on AIC.
                                 update.init = TRUE, init.dist = "best",
                                 data = boot_dat))
            
          } else {
            
            fit <- try(smoothSurvReg(Surv(shortest, longest, type = "interval2") ~ 1, aic = TRUE,
                                 # the optimal value of λ is found using a grid search, based on AIC.
                                 update.init = TRUE, init.dist = "best",
                                 data = boot_dat))
            
            }
          #  smoothSurvReg.control(difforder = 2))  # m (the order of the finite difference used in
          # the penalty term in smoothSurvReg.control) is 3 by default
        ))
        
       # End repeat. This is to make sure that a single error in one bootstrap run does not stop the whole system.
      
      
       if( (class(fit) == "try-error") == TRUE ){rep(NA, 5)} else {
       S <- survfit(fit, by = 0.0001, xlim = c(0, 30), plot = FALSE)
       sapply(percentiles, FUN = function(p){S$x[which(S$y1 >= (1 - p)) %>% tail(1)]}) 
        }
       # ) Percentiles
    
    }
    
    # bootstrapping with 1000 replications
    results <- boot(data = dat, statistic = fun_boot, R = n.boots)
    
    # Remove all instances of NA, see https://stat.ethz.ch/pipermail/r-help/2007-September/141868.html
    results.naomit <- results
    cat("NA estimates in bootstrap runs per percentile:", colSums(is.na(results$t)))
    results.naomit$t0 = na.omit(results.naomit$t0)
    results.naomit$t = na.omit(results.naomit$t)
    results.naomit$R = length(results.naomit$t[,1])
    results.naomit$call[4] = results.naomit$R
    cat("\n So estimates are based on ", results.naomit$R, " bootstrap samples.")
    
    # get 95% confidence interval
    one <- boot.ci(results.naomit, type="perc", index = 1)
    two <- boot.ci(results.naomit, type="perc", index = 2)
    three <- boot.ci(results.naomit, type="perc", index = 3)
    four <- boot.ci(results.naomit, type="perc", index = 4)
    five <- boot.ci(results.naomit, type="perc", index = 5)
    
    boot_CIs <- matrix(NA, nrow = 5, ncol = 2)
    boot_CIs[1,] <- one$percent[c(4,5)]
    boot_CIs[2,] <- two$percent[c(4,5)]
    boot_CIs[3,] <- three$percent[c(4,5)]
    boot_CIs[4,] <- four$percent[c(4,5)]
    boot_CIs[5,] <- five$percent[c(4,5)]
    
    colnames(boot_CIs) <- c("lower_CI", "upper_CI")
    
    tab <- cbind(tab, boot_CIs)
    
 } ##### End of CIs, method 2
  
  
  
  # Output
  if(pdf_cdf == TRUE){list(tab, pdf, S)} else {tab}
  
}

# Survival package NPMLE

# Function to estimate incubation time from 
# single interval censored 
# observations 
# using nonparametric approach in survival package.
Estimate_single_npmle_survival <- function(dat, cdf = FALSE, conf.int = FALSE){  # L-BFGS-B, BFGS, SANN, Nelder-Mead, CG, Brent

  require('survival')
  require('tidyr')
  require('magrittr')
  require('dplyr')
  require('Icens')

  # Inputs:
  # dat     data frame with columns L0, L1, R0, R1

  shortest <- dat$R0 - dat$L1
  longest <- dat$R0 - dat$L0
  
  # NPMLE
  # Each observation as a time interval with (-infinity, t) for left censored, (t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
  
  # Using survival package:
  fit <- survfit(Surv(shortest, longest, type='interval2') ~ 1, data = dat, conf.type = "plain", 
                 conf.int = 0.95) # NPMLE, with pointwise CI: https://www2.karlin.mff.cuni.cz/~vavraj/cda/exercise_04.html
  
  # Confidence limits for the values are based on the intersection of the horizontal line at 1-k with the upper and lower limits for the survival curve. Hence confidence limits use the same p-value as was in effect when the curve was created, and will differ depending on the conf.type option of survfit. If the survival curves have no confidence bands, confidence limits for the quantiles are not available.
  # When a horizontal segment of the survival curve exactly matches one of the requested quantiles the returned value will be the midpoint of the horizontal segment; this agrees with the usual definition of a median for uncensored data. 
  
  if(conf.int == FALSE){
  estimates <- quantile(fit, probs = percentiles, conf.int = F) # https://stat.ethz.ch/R-manual/R-devel/library/survival/html/quantile.survfit.html
  tab <- data.frame(est = estimates)

  rownames(tab) <- row_names_percentiles

  } else if(conf.int == "survfit"){
    estimates <- quantile(fit, probs = percentiles, conf.int = T) # https://stat.ethz.ch/R-manual/R-devel/library/survival/html/quantile.survfit.html
    tab <- data.frame(est = estimates$quantile, lower_CI = estimates$lower, upper_CI = estimates$upper)

    rownames(tab) <- row_names_percentiles

  } else if(conf.int == "bootstrap"){
    estimates <- quantile(fit, probs = percentiles, conf.int = F) # https://stat.ethz.ch/R-manual/R-devel/library/survival/html/quantile.survfit.html
        tab <- data.frame(est = estimates$quantile)

    rownames(tab) <- row_names_percentiles

    n.boots <- 1000

    for_CI <- matrix(nrow = 5, ncol = n.boots)

    for(boot in 1:n.boots){
    boot_dat <- dat[sample(1:round(sqrt(nrow(dat))), replace = T),] # m out of n bootstrap

    # Model fit
    shortest <- boot_dat$R0 - boot_dat$L1
    longest <- boot_dat$R0 - boot_dat$L0

    # NPMLE
    # Each observation as a time interval with (-infinity, t) for left censored, (t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
    fit <- survfit(Surv(shortest, longest, type='interval2') ~ 1, data = boot_dat, conf.type = "plain",
                   conf.int = 0.95) # NPMLE, with pointwise CI: https://www2.karlin.mff.cuni.cz/~vavraj/cda/exercise_04.html

    for_CI[ , boot] <- quantile(fit, probs = percentiles, conf.int = F) # https://stat.ethz.ch/R-manual/R-devel/library/survival/html/quantile.survfit.html
    }

    boot_CIs <- t(apply(for_CI, 1, FUN = sort))[ , c(25,975)]
    colnames(boot_CIs) <- c("lower_CI", "upper_CI")

    tab <- cbind(tab, boot_CIs)

  }
  #### End of CIs
  
# Output
 if(cdf == TRUE){
  
   cdf <- data.frame(t = fit$time, S = fit$surv)
   list(tab, cdf)
   
 } else{
  tab
 }

}


# Flexsurv package MLE (for gamma)

# Function to estimate incubation time from 
# single interval censored 
# observations 
# using nonparametric approach in survival package.
Estimate_single_mle_flexsurv <- function(dat, distribution = "gamma", par_estimates = FALSE, 
                                         init.vals = NA){  
  
  # Inputs:
  # dat     data frame with columns L0, L1, R0, R1

  tab <- as.data.frame(matrix(nrow = length(percentiles), ncol = 3))
  
  shortest <- dat$R0 - dat$L1
  longest <- dat$R0 - dat$L0
  
  require('survival')
  require('flexsurv')
  require('tidyr')
  require('magrittr')
  require('dplyr')

  tab <- as.data.frame(matrix(nrow = length(percentiles), ncol = 3))

  if(distribution == "lognormal"){distribution <- "lnorm"}
  if(distribution == "Weibull"){distribution <- "weibullPH"}
    
  # MLE
  # Each observation as a time interval with (-infinity, t) for left censored, (t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
  shortest[shortest == 0] <- shortest[shortest == 0] + 0.000001 # To make the function work.
  
  if(length(init.vals) == 1){
  fit <- flexsurvreg(data = dat, Surv(shortest, longest, type='interval2') ~ 1, dist = distribution)
  } else {
    fit <- flexsurvreg(data = dat, Surv(shortest, longest, type='interval2') ~ 1, dist = distribution, 
                       inits = init.vals)
  }
  
  if(par_estimates == TRUE){
  
    params <-  fit$res[,1]
    names(params) <- c("shape", "rate")  
    
  params} else{

# p.20 possible distributions: https://cran.r-project.org/web/packages/flexsurv/flexsurv.pdf     

 # From: https://stats.stackexchange.com/questions/497735/formula-for-predicting-median-survival-using-a-royston-parmar-model
  tab <- summary(fit, type = "quantile", quantiles = percentiles) %>% as.data.frame()
  
  rownames(tab) <- row_names_percentiles
  colnames(tab) <- c("quantile", "est", "lower_CI", "upper_CI")
  tab <- tab %>% select(est, lower_CI, upper_CI)
  
tab
#  list(tab, params)
}
  
}


# using nonparametric approach in survival package.
Estimate_single_mle_survival <- function(dat, distribution = "lognormal", par_estimates = FALSE){  
  
  # Inputs:
  # dat     data frame with columns L0, L1, R0, R1
  
  shortest <- dat$R0 - dat$L1
  longest <- dat$R0 - dat$L0
  
  require('survival')
  require('tidyr')
  require('magrittr')
  require('dplyr')
  
  tab <- as.data.frame(matrix(nrow = length(percentiles), ncol = 3))
  
  # MLE
  # Each observation as a time interval with (-infinity, t) for left censored, (t, infinity) for right censored, (t,t) for exact and (t1, t2) for an interval. This is the approach used for type = interval2.
  shortest[shortest == 0] <- shortest[shortest == 0] + 0.000001 # To make the function work.
  
  if(distribution == "Weibull"){
    fit <- survreg(data = dat, Surv(shortest, longest, type='interval2') ~ 1, dist = "weibull")
    summary(fit)
    
  } else if(distribution == "lognormal"){ # gamma not available.
    fit <- survreg(data = dat, Surv(shortest, longest, type='interval2') ~ 1, dist = "lognormal")
    # meanlog <- fit$coefficients[[1]]
    #  sdlog <- fit$scale
    #  est_vec <- qlnorm(p = percentiles, meanlog = meanlog, sdlog = sdlog)
    
    # from ??survreg.distributions: survreg scale parameter maps to 1/shape, linear predictor to log(scale)
    # and https://stat.ethz.ch/pipermail/r-help/2008-February/155127.html
    
    # for(perc in 1:length(percentiles)){ # Based on https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf page 85.
    # q1 <- predict(fit, type = 'quantile',
    #               p = percentiles[perc], se.fit = T)
    # ci1 <- cbind(q1$fit, q1$fit - 1.96*q1$se.fit, q1$fit + 1.96*q1$se.fit)
    # tab[perc,] <- ci1[1,]
    # }
    
    #  params <- c(fit$icoef[1], exp(fit$icoef[2]))
    #  names(params) <- c("shape", "scale")
  }
  
  if(par_estimates == TRUE){
    
    if(distribution == "lognormal"){
      params <- c(fit$coefficients, fit$scale)
      names(params) <- c("meanlog", "sdlog")  
      
    } else if(distribution == "Weibull"){
      params <- c(1/fit$scale, exp(fit$coefficients))
      names(params) <- c("shape", "scale")
    }
    
    params} else{
      
      for(perc in 1:length(percentiles)){ # Based on https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf page 85.
        q2 <- predict(fit, type = 'uquantile',
                      p = percentiles[perc], se.fit = T)
        ci2 <- cbind(q2$fit, q2$fit - 1.96*q2$se.fit, q2$fit + 1.96*q2$se.fit)
        tab[perc,] <- exp(ci2[1,]) # convert from log scale to original y
      }
      
      rownames(tab) <- row_names_percentiles
      colnames(tab) <- c("est", "lower_CI", "upper_CI")
      
      tab
      #  list(tab, params)
    }
  
}



# Reich

# Function to estimate incubation time from 
# single interval censored 
# observations 
# using parametric approach.
# Using Reich's package coarseDataTools.
Estimate_single_Reich <- function(dat, distribution = "lognormal", 
                                  opt_method = "BFGS"){  # L-BFGS-B, BFGS, SANN, Nelder-Mead, CG, Brent
  
  require('tidyr')
  require('magrittr')
  require('dplyr')
  require('coarseDataTools')
  
  # Inputs:
  # dat     data frame with columns L0, L1, R0, R1
  
  #dat <- danang_n92

  # Renaming columns
  dat <- dat %>% transmute(
    ER = L1,
    SL = R0,
    SR = R1,
    EL = L0
  ) %>% 
    
    # Specifying censoring types
    mutate(type = 1) %>%
    as.matrix()
  
  
  
  # Fit models
  if(distribution == "lognormal"){
    fit <- dic.fit(dat = dat, dist = "L", opt.method = opt_method, ptiles = c(0.025, 0.05, 0.5, 0.95, 0.975))  # Lognormal distribution
  } else if(distribution == "gamma"){
    fit <- dic.fit(dat = dat, dist = "G", n.boots = 1000, opt.method = opt_method, ptiles = c(0.025, 0.05, 0.5, 0.95, 0.975))  # Gamma "
  } else if(distribution == "Weibull"){
    fit <- dic.fit(dat = dat, dist = "W", n.boots = 1000, 
                   ptiles = c(0.025, 0.05, 0.5, 0.95, 0.975), opt.method = opt_method)  # Weibull 
  }
  
  # Print table with censoring types
  # cat(" # doubly interval-censored: ", dat %>% filter(type == 0) %>% nrow(), "\n",
  #     "# single interval-censored: ", dat %>% filter(type == 1) %>% nrow(), "\n",
  #     "# exact observation: ", dat %>% filter(type == 2) %>% nrow(), "\n" )
  
  #fit@samples # If bootstrap is used this contains the bootstrapped samples of parameters.
  par1 <- fit@ests[1]
  par2 <- fit@ests[2]
  # One can get the confidence bounds also.
  
  xvals <- seq(0, 30, by = .1)
  
  ##### Plot
  if(distribution == "lognormal"){
    df_p <- data.frame(
      xvals = seq(0, 30, by = .1),
      pdf = dlnorm(xvals, meanlog = par1, sdlog = par2)
    )
    
  } else if(distribution == "gamma"){
    df_p <- data.frame(
      xvals = seq(0, 30, by = .1),
      pdf = dgamma(xvals, shape = par1, scale = par2)
    )
    
  } else if(distribution == "Weibull"){
    df_p <- data.frame(
      xvals = seq(0, 30, by = .1),
      pdf = dweibull(xvals, shape = par1, scale = par2)
    )
  }
  
  #### Collect output
  tab <- fit@ests
  # tab[1:2,1] <- exp(tab[1:2,1])
  tab <- tab[-1:-2,-4]
  rownames(tab) <- row_names_percentiles
  colnames(tab) <- c("est", "lower_CI", "upper_CI")
  ####
  
  # inc_or_lat_time <- setClass("incubation or latency time", 
  #                             slots = c(ests = "data.frame", 
  #                                       dist = "character", 
  #                                       par = "vector", 
  #                                       pdf = "data.frame", 
  #                                       CI_method = "character")
  # )
  
  pdf <- df_p %>% as.data.frame()
  colnames(pdf) <- c("t", "p")
  
  plot(pdf$t, pdf$p, type = "l")
  
  # out <- inc_or_lat_time(ests = tab %>% as.data.frame(), 
  #                        dist = distribution, 
  #                        par = c(par1, par2),
  #                        pdf = pdf,
  #                        CI_method = "")
  
  tab
  
}

# Groeneboom

# Function to estimate incubation time from 
# single interval censored observations using
# nonparametric approach with ICM algorithm.

# Based on Groeneboom, 2020 (Github pietg).
Estimate_single_Groeneboom <- function(dat){
  
  # Inputs:
  # dat     data.frame with columns L0, L1, R0, R1

  require('survival')
  require('tidyr')
  require('magrittr')
  require('dplyr')
  require('Rcpp')
  
  sourceCpp("NPMLE_ICM.cpp")
  
  if(sum(dat$R0 == dat$R1) != nrow(dat)){
    warning("Symptom onset is not exact for at least one individual. 
                 Use other method.")
  }
  
  dat <- dat %>% mutate(R = R0)
  
  A <- data.frame(V1 = dat$R - dat$L1, V2 = dat$R - dat$L0)
  # Yields data frame A with column 1 is S-E1 and column 2 is S-E0.
  
  output1 <- NPMLE(A)
  
  C1 <- output1$SMLE # smoothed nonparametric likelihood estimate, cdf
  D1 <- output1$dens  # pdf
  
  t<-D1[,1]
  u<-D1[,2]
  plot(t, u, type = "l")
  
  # Based on smoothed nonparametric cdf
  tab <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975), 
                function(p) max(C1[which(C1[,2]<=p),1]) )
  
  #### Collect output
  tab <- matrix(tab, nrow = 5, ncol = 1)
  rownames(tab) <- row_names_percentiles
  colnames(tab) <- c("est")
  ####
  
  # BOOTSTRAPPED CIs
  # 
  # mat <- matrix(nrow = 5, ncol = 1000)
  # rownames(mat) <- c("q2.5", "q5", "q50", "q95", "q975")
  # 
  # for (i in 1:1000){
  #   set.seed(i)
  #   samp = sample(nrow(dat),replace=T)
  #   dat_boot = dat[samp,]
  #   
  #   A <- data.frame(V1 = dat_boot$R - dat_boot$L1, V2 = dat_boot$R - dat_boot$L0)
  #   # Yields data frame A with column 1 is S-E1 and column 2 is S-E0.
  #   
  #   output1 <- NPMLE(A)
  #   C1 <- output1$SMLE # smoothed nonparametric likelihood estimate, cdf
  #   
  #   # Based on smoothed nonparametric cdf
  #   mat[,i] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975), 
  #                     function(p) max(C1[which(C1[,2]<=p),1]) )
  #   
  # } ## END BOOTSTRAPPED CIs
  # 
  # CIs <- t( apply(mat, 1, FUN = function(x) sort(x)[c(26,975)]) ) %>% as.data.frame()
  # colnames(CIs) <- c("lower_CI", "upper_CI")
  # 
  # 
  
  
  # inc_or_lat_time <- setClass("incubation or latency time", 
  #                             slots = c(ests = "data.frame", 
  #                                       dist = "character", 
  #                                       par = "vector", 
  #                                       pdf = "data.frame", 
  #                                       CI_method = "character")
  # )
  # 
  # pdf <- cbind(t, u) %>% as.data.frame()
  # colnames(pdf) <- c("t", "p")
  # 
  # out <- inc_or_lat_time(ests = cbind(tab, CIs), 
  #                        dist = "nonparametric", 
  #                        par = NA,
  #                        pdf = pdf,
  #                        CI_method = "bootstrap")
  # 
  # out
  
  tab
}


