# Calculates estimates based on Deng's and Qin's data generation process.



Estimate_renewal_Pi <- function(dat, assumed_distribution = "gamma", percentiles = c(0.25, 0.5, 0.9, 0.95, 0.975, 0.99)){
 
# percentiles is a vector of length 6
 # aa <- par1
 # bb <- par2
 # pp <- 1-pi
  x <- dat
  
  if(assumed_distribution == "gamma"){
    ##########################################gamma incubation
    f <- function(y,a,b){
      return(dgamma(y,shape=a,rate=b))
    }
    
    #pdf of forward time
    h <- function(y,a,b){
      RE = b/a*(1-pgamma(y,shape=a,rate=b))
      return(RE)
    }
    FI <- function(y,a,b) return(pgamma(y,shape=a,rate=b))
    
    F <- function(k,a,b,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REh = pgamma(k[posit],shape=a+1,rate=b)+k[posit]*b/a*(1-pgamma(k[posit],shape=a,rate=b))
      if (p>0.999) RE[posit] = REh
      else{
        REf = pgamma(k[posit],shape=a,rate=b)
        RE[posit] = REf*(1-p)+REh*p
      }
      return(RE)
    }
    
    
    
  }else if(assumed_distribution == "Weibull"){
    ##########################################weibull incubation
    f <- function(y,a,b){
      return(dweibull(y,shape=a,scale=b))
    }
    h <- function(y,a,b){
      RE = (1-pweibull(y,shape=a,scale=b))/gamma(1+1/a)/b
      return(RE)
    }
    FI <- function(y,a,b) return(pweibull(y,shape=a,scale=b))
    F <- function(k,a,b,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REh = pgamma((k[posit]/b)^a,shape=1/a,rate=1)
      if (p>0.999) RE[posit] = REh
      else{
        REf = pweibull(k[posit],shape=a,scale=b)
        RE[posit] = REf*(1-p)+REh*p
      }
      return(RE)
    }
    
  }else if(assumed_distribution == "lognormal"){
    ##########################################lognormal incubation
    f <- function(y,u,s){
      RE = rep(0,length(y))
      posit = which(y>0)
      RE[posit]=dnorm(log(y[posit]),u,s)/y[posit]
      return(RE)
    }
    h <- function(y,u,s){
      RE = exp(-u-s^2/2)*(1-pnorm(log(y),u,s))
      return(RE)
    }
    FI <- function(y,u,s){
      RE = rep(0,length(y))
      posit = which(y>0)
      RE[posit]=pnorm(log(y[posit]),u,s)
      return(RE)
    }
    F <- function(k,u,s,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REf = pnorm(log(k[posit]),u,s)
      REh = pnorm(log(k[posit]),u+s^2,s)+exp(-u-s^2/2)*k[posit]*(1-pnorm(log(k[posit]),u,s))
      RE[posit] = REf*(1-p)+REh*p  # equation 2 of paper Deng
      return(RE)
    }
    
  }
  
  loglik_qin <- function(pa){
    a = pa[1]
    b = pa[2]
    p = pa[3]
    if (p>0.999) P=h(x,a,b)
    else P = p*h(x,a,b)+(1-p)*f(x,a,b)
    P[P<0.00001]=0.00001
    RE = - sum( log(P) )
    if (RE>10000) RE=10000
    return(RE)
  }
  
  
  loglik_deng <- function(pa){
    a = pa[1]
    b = pa[2]
    p = pa[3]
    P = F(x+0.5,a,b,p)-F(x-0.5,a,b,p)
    P[P<0.00001]=0.00001
    RE = - sum( log(P) )
    if (RE>10000) RE=10000
    return(RE)
  }
  
      
      # ESTIMATION Deng
      
      if(assumed_distribution == "gamma"){
        par_1 <- optim(par=c(4,0.5,0.8),loglik_deng, method='L-BFGS-B',
                       lower=c(1.1,0.1,0),upper=c(10,2,1))$par
      } else if(assumed_distribution == "Weibull"){
        par_1 <- optim(par=c(2,10,0.8), loglik_deng,method='L-BFGS-B',
                       lower=c(1.1,2,0),upper=c(5,15,1))$par
        
      } else if(assumed_distribution == "lognormal"){
        par_1 <- optim(par=c(2,0.4,0.8), loglik_deng,method='L-BFGS-B',
                       lower=c(1,0.1,0),upper=c(5,1,1))$par
      }
      

      a <- par_1[1]
      b <- par_1[2]
      p <- par_1[3]
 
      #mean, median, quartiles
      if(assumed_distribution == "gamma"){
        #########################gamma
        est <- c(
          c(a,b,1-p) ,  # estimates for parameters alpha and beta, and pi
          a/b,         # mean
          qgamma(percentiles,a,b))  # estimated quantiles
        # -round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "Weibull"){
        #########################weibull
        est <- c(
          c(a,b,1-p) ,   # estimates for parameters k and lambda, and pi
          b*gamma(1+1/a) , # mean
          qweibull(percentiles,shape=a,scale=b))  # estimated quantiles
        #-round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "lognormal"){
        ########################lognormal
        est <- c(
          c(a,b,1-p),  # estimates for parameters mu and sigma, and pi
          exp(a+b^2/2), # mean
         exp(qnorm(percentiles,a,b)))  # estimated quantiles
        #-round(loglik(c(a,b)),2)
        
      }
      
      #mean, median, quartiles
vec_ests <- est[c(1, 2, 3, 4:10)]
names(vec_ests) <- c("par1_deng", "par2_deng", "pi_deng", "mean_deng", "Q1_deng", "Q2_deng", "Q3_deng", "Q4_deng", "Q5_deng", "Q6_deng") 
      
      # ESTIMATION Qin
      
      if(assumed_distribution == "gamma"){
        par_2 <- optim(par=c(4,0.5,0.8),loglik_qin,method='L-BFGS-B',
                       lower=c(1.1,0.1,0),upper=c(10,2,1))$par
      }
      else if(assumed_distribution == "Weibull"){
        par_2 <- optim(par=c(2,10,0.8),loglik_qin,method='L-BFGS-B',
                       lower=c(1.1,2,0),upper=c(5,15,1))$par
      }
      else if(assumed_distribution == "lognormal"){
        par_2 <- optim(par=c(2,0.4,0.8),loglik_qin,method='L-BFGS-B',
                       lower=c(1,0.1,0),upper=c(5,1,1))$par
      }
      
      #model parameter and pi
      
         #mean, median, quartiles
      if(assumed_distribution == "gamma"){
        #########################gamma
        est <- c(
          c(a,b,1-p),  # estimates for parameters alpha and beta, and pi
          a/b,         # mean
          qgamma(percentiles,a,b))  # estimated quantiles
        # -round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "Weibull"){
        #########################weibull
        est <- c(
          c(a,b,1-p),   # estimates for parameters k and lambda, and pi
          b*gamma(1+1/a), # mean
          qweibull(percentiles,shape=a,scale=b))  # estimated quantiles
        #-round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "lognormal"){
        ########################lognormal
        est <- c(
          c(a,b,1-p),  # estimates for parameters mu and sigma, and pi
          exp(a+b^2/2), # mean
         exp(qnorm(percentiles,a,b)))  # estimated quantiles
        #-round(loglik(c(a,b)),2)
        
      }
      
      
      #mean, median, quartiles, bias
      
vec_ests_qin <- est[c(1, 2, 3, 4:10)]
names(vec_ests_qin) <- c("par1_qin", "par2_qin", "pi_qin","mean_qin", "Q1_qin", "Q2_qin", "Q3_qin", "Q4_qin", "Q5_qin", "Q6_qin") 

return(c(vec_ests, vec_ests_qin))
    
   
}

#out <- t(sapply(1:10, function(i){Estimate_renewal_Pi(dat = list_dat[[i]], assumed_distribution = "gamma")}))

# ----------------


# Calculates estimates based on exp_data_generation.R data generation process, when dropping pi (i.e., extra infections at day of travel).
##   
#for (w in 1:3){
#  if (w%%3==1) m=600
#  if (w%%3==2) m=1200
#  if (w%%3==0) m=1800

Estimate_renewal_noPi  <- function(dat, assumed_distribution = "gamma", percentiles = c(0.25, 0.5, 0.9, 0.95, 0.975, 0.99)){

  
  x <- dat
  
  if(assumed_distribution == "gamma"){
    ##########################################gamma incubation
    f <- function(y,a,b){
      return(dgamma(y,shape=a,rate=b))
    }
    h <- function(y,a,b){
      RE = b/a*(1-pgamma(y,shape=a,rate=b))
      return(RE)
    }
    FI <- function(y,a,b) return(pgamma(y,shape=a,rate=b))
    F <- function(k,a,b,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REh = pgamma(k[posit],shape=a+1,rate=b)+k[posit]*b/a*(1-pgamma(k[posit],shape=a,rate=b))
      if (p>0.999) RE[posit] = REh
      else{
        REf = pgamma(k[posit],shape=a,rate=b)
        RE[posit] = REf*(1-p)+REh*p
      }
      return(RE)
    }
    
    
    
  }else if(assumed_distribution == "Weibull"){
    ##########################################weibull incubation
    f <- function(y,a,b){
      return(dweibull(y,shape=a,scale=b))
    }
    h <- function(y,a,b){
      RE = (1-pweibull(y,shape=a,scale=b))/gamma(1+1/a)/b
      return(RE)
    }
    FI <- function(y,a,b) return(pweibull(y,shape=a,scale=b))
    F <- function(k,a,b,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REh = pgamma((k[posit]/b)^a,shape=1/a,rate=1)
      if (p>0.999) RE[posit] = REh
      else{
        REf = pweibull(k[posit],shape=a,scale=b)
        RE[posit] = REf*(1-p)+REh*p
      }
      return(RE)
    }
    
  }else if(assumed_distribution == "lognormal"){
    ##########################################lognormal incubation
    f <- function(y,u,s){
      RE = rep(0,length(y))
      posit = which(y>0)
      RE[posit]=dnorm(log(y[posit]),u,s)/y[posit]
      return(RE)
    }
    h <- function(y,u,s){
      RE = exp(-u-s^2/2)*(1-pnorm(log(y),u,s))
      return(RE)
    }
    FI <- function(y,u,s){
      RE = rep(0,length(y))
      posit = which(y>0)
      RE[posit]=pnorm(log(y[posit]),u,s)
      return(RE)
    }
    F <- function(k,u,s,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REh = pnorm(log(k[posit]),u+s^2,s)+exp(-u-s^2/2)*k[posit]*(1-pnorm(log(k[posit]),u,s))
      if (p>0.999) RE[posit]=REh
      else {
        REf = pnorm(log(k[posit]),u,s)
        RE[posit] = REf*(1-p)+REh*p
      }
      return(RE)
    }
  }
  
  
  loglik0 <- function(pa){
    a = pa[1]
    b = pa[2]
    P = h(x,a,b)
    P[P<0.00001]=0.00001
    RE = - sum( log(P) )
    if (RE>10000) RE=10000
    return(RE)
  }
  
  
  Loglik0 <- function(pa){
    a = pa[1]
    b = pa[2]
    P = F(x+0.5,a,b,1)-F(x-0.5,a,b,1)
    P[P<0.00001]=0.00001
    RE = - sum( log(P) )
    if (RE>10000) RE=10000
    return(RE)
  }
  
        if(assumed_distribution == "gamma"){
        par_1 <- optim(par=c(4,0.5),Loglik0,method='L-BFGS-B',
                       lower=c(1.1,0.1),upper=c(10,2))$par
      }
      else if(assumed_distribution == "Weibull"){
        par_1 <- optim(par=c(2,10),Loglik0,method='L-BFGS-B',
                       lower=c(1.1,2),upper=c(5,15))$par
        
      }
      else if(assumed_distribution == "lognormal"){
        par_1 <- optim(par=c(2,0.4),Loglik0,method='L-BFGS-B',
                       lower=c(1,0.1),upper=c(5,1))$par
      }
      
      #model parameter
      
      a = par_1[1]
      b = par_1[2]
      
      if(assumed_distribution == "gamma"){
        #########################gamma
        est <- c(
          c(a,b) ,  # estimates for parameters alpha and beta, and pi
          a/b,         # mean
          qgamma(percentiles,a,b))  # estimated quantiles
        # -round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "Weibull"){
        #########################weibull
        est <- c(
          c(a,b) ,   # estimates for parameters k and lambda, and pi
          b*gamma(1+1/a) , # mean
          qweibull(percentiles,shape=a,scale=b))  # estimated quantiles
        #-round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "lognormal"){
        ########################lognormal
        est <- c(
          c(a,b),  # estimates for parameters mu and sigma, and pi
          exp(a+b^2/2), # mean
         exp(qnorm(percentiles,a,b)))  # estimated quantiles
        #-round(loglik(c(a,b)),2)
      }  
      
      vec_ests <- est[c(1, 2, 3:9)]
      names(vec_ests) <- c("par1_deng", "par2_deng", "mean_deng", "Q1_deng", "Q2_deng", "Q3_deng", "Q4_deng", "Q5_deng", "Q6_deng") 
      
      
      if(assumed_distribution == "gamma"){
        par_2 <- optim(par=c(4,0.5,0.8),loglik0,method='L-BFGS-B',
                       lower=c(1.1,0.1,0),upper=c(10,2,1))$par
      }
      else if(assumed_distribution == "Weibull"){
        par_2 <- optim(par=c(2,10,0.8),loglik0,method='L-BFGS-B',
                       lower=c(1.1,2,0),upper=c(5,15,1))$par
      }
      else if(assumed_distribution == "lognormal"){
        par_2 <- optim(par=c(2,0.4,0.8),loglik0,method='L-BFGS-B',
                       lower=c(1,0.1,0),upper=c(5,1,1))$par
      }
      
      #model parameter and pi
      
      a = par_2[1]
      b = par_2[2]
      
      #mean, median, quartiles
      if(assumed_distribution == "gamma"){
        #########################gamma
        est <- c(
          c(a,b),  # estimates for parameters alpha and beta, and pi
          a/b,         # mean
          qgamma(percentiles,a,b))  # estimated quantiles
        # -round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "Weibull"){
        #########################weibull
        est <- c(
          c(a,b),   # estimates for parameters k and lambda, and pi
          b*gamma(1+1/a), # mean
         qweibull(percentiles,shape=a,scale=b))  # estimated quantiles
        #-round(loglik(c(a,b)),2),
        
      } else if(assumed_distribution == "lognormal"){
        ########################lognormal
        est <- c(
          c(a,b),  # estimates for parameters mu and sigma, and pi
          exp(a+b^2/2), # mean
          exp(qnorm(percentiles,a,b)))  # estimated quantiles
        #-round(loglik(c(a,b)),2)
      }
      

      vec_ests_qin <- est[c(1, 2, 3:9)]
      names(vec_ests_qin) <- c("par1_qin", "par2_qin", "mean_qin", "Q1_qin", "Q2_qin", "Q3_qin", "Q4_qin", "Q5_qin", "Q6_qin") 
      
      return(c(vec_ests, vec_ests_qin))
      

    }  


Summarize_results <- function(results, percentiles = c(0.25, 0.5, 0.9, 0.95, 0.975, 0.99),
                              scenario){

  distribution_T <- scenario$distribution_T
  par1 <- scenario$par1
  par2 <- scenario$par2
  estimation_method <- scenario$estimation_method
  
  out <- as.data.frame(results)
  
  ### Real quantiles
  if(distribution_T == "lognormal"){
    q_real <- qlnorm(percentiles, meanlog = par1, sdlog = par2)} else if(distribution_T == "Weibull"){
      q_real <- qweibull(percentiles, shape = par1, scale = par2)} else if(distribution_T == "heavytail"){
        #library(rmutil)
        #q_real <- qburr(percentiles, m = 8.5, s = 2, f = 2)
        q_real <- c(3.343219, 5.470551, 12.498982, 15.838618, 19.613748, 25.500000)}
  
  # Bias
  out$bias_Q1_deng <- (out$Q1_deng - q_real[1]); out$bias_Q2_deng <- (out$Q2_deng - q_real[2]); 
  out$bias_Q3_deng <- (out$Q3_deng - q_real[3]); out$bias_Q4_deng <- (out$Q4_deng - q_real[4]); 
  out$bias_Q5_deng <- (out$Q5_deng - q_real[5]); out$bias_Q6_deng <- (out$Q6_deng - q_real[6])
  
  out$bias_Q1_qin <- (out$Q1_qin - q_real[1]); out$bias_Q2_qin <- (out$Q2_qin - q_real[2]); 
  out$bias_Q3_qin <- (out$Q3_qin - q_real[3]); out$bias_Q4_qin <- (out$Q4_qin - q_real[4]); 
  out$bias_Q5_qin <- (out$Q5_qin - q_real[5]); out$bias_Q6_qin <- (out$Q6_qin - q_real[6])
  
  # Mean bias and q2.5 and q97.5
  mean_bias_Q1_deng <- mean(out$bias_Q1_deng)
  mean_bias_Q2_deng <- mean(out$bias_Q2_deng)
  mean_bias_Q3_deng <- mean(out$bias_Q3_deng)
  mean_bias_Q4_deng <- mean(out$bias_Q4_deng)
  mean_bias_Q5_deng <- mean(out$bias_Q5_deng)
  mean_bias_Q6_deng <- mean(out$bias_Q6_deng)
  
  mean_bias_Q1_deng_q25 <- quantile(out$bias_Q1_deng, probs = 0.025); mean_bias_Q2_deng_q25 <- quantile(out$bias_Q2_deng, probs = 0.025)
  mean_bias_Q3_deng_q25 <- quantile(out$bias_Q3_deng, probs = 0.025); mean_bias_Q4_deng_q25 <- quantile(out$bias_Q4_deng, probs = 0.025)
  mean_bias_Q5_deng_q25 <- quantile(out$bias_Q5_deng, probs = 0.025); mean_bias_Q6_deng_q25 <- quantile(out$bias_Q6_deng, probs = 0.025)
  mean_bias_Q1_deng_q975 <- quantile(out$bias_Q1_deng, probs = 0.975); mean_bias_Q2_deng_q975 <- quantile(out$bias_Q2_deng, probs = 0.975)
  mean_bias_Q3_deng_q975 <- quantile(out$bias_Q3_deng, probs = 0.975); mean_bias_Q4_deng_q975 <- quantile(out$bias_Q4_deng, probs = 0.975)
  mean_bias_Q5_deng_q975 <- quantile(out$bias_Q5_deng, probs = 0.975); mean_bias_Q6_deng_q975 <- quantile(out$bias_Q6_deng, probs = 0.975)
  
  mean_bias_Q1_qin <- mean(out$bias_Q1_qin)
  mean_bias_Q2_qin <- mean(out$bias_Q2_qin)
  mean_bias_Q3_qin <- mean(out$bias_Q3_qin)
  mean_bias_Q4_qin <- mean(out$bias_Q4_qin)
  mean_bias_Q5_qin <- mean(out$bias_Q5_qin)
  mean_bias_Q6_qin <- mean(out$bias_Q6_qin)
  
  mean_bias_Q1_qin_q25 <- quantile(out$bias_Q1_qin, probs = 0.025); mean_bias_Q2_qin_q25 <- quantile(out$bias_Q2_qin, probs = 0.025)
  mean_bias_Q3_qin_q25 <- quantile(out$bias_Q3_qin, probs = 0.025); mean_bias_Q4_qin_q25 <- quantile(out$bias_Q4_qin, probs = 0.025)
  mean_bias_Q5_qin_q25 <- quantile(out$bias_Q5_qin, probs = 0.025); mean_bias_Q6_qin_q25 <- quantile(out$bias_Q6_qin, probs = 0.025)
  mean_bias_Q1_qin_q975 <- quantile(out$bias_Q1_qin, probs = 0.975); mean_bias_Q2_qin_q975 <- quantile(out$bias_Q2_qin, probs = 0.975)
  mean_bias_Q3_qin_q975 <- quantile(out$bias_Q3_qin, probs = 0.975); mean_bias_Q4_qin_q975 <- quantile(out$bias_Q4_qin, probs = 0.975)
  mean_bias_Q5_qin_q975 <- quantile(out$bias_Q5_qin, probs = 0.975); mean_bias_Q6_qin_q975 <- quantile(out$bias_Q6_qin, probs = 0.975)
  
  
  # Mean pi and confidence interval
  if(estimation_method == "nopi"){
    mean_pi_deng <- mean_pi_qin <- mean_pi_deng_lowerCI <- mean_pi_deng_upperCI <- mean_pi_qin_lowerCI <- mean_pi_qin_upperCI <- NA
  } else{
    mean_pi_deng <- mean(out$pi_deng)
    mean_pi_qin <- mean(out$pi_qin)
    
    CI_meanpi_deng <- t.test(out$pi_deng)$conf.int # Normal approximation for CI.
    CI_meanpi_qin <- t.test(out$pi_qin)$conf.int
    
    mean_pi_deng_lowerCI <- CI_meanpi_deng[1]
    mean_pi_deng_upperCI <- CI_meanpi_deng[2]
    mean_pi_qin_lowerCI <- CI_meanpi_qin[1]
    mean_pi_qin_upperCI <- CI_meanpi_qin[2]
  }
  
  vec_results <- c(mean_pi_deng, mean_pi_deng_lowerCI, mean_pi_deng_upperCI,
                   mean_bias_Q1_deng, mean_bias_Q1_deng_q25, mean_bias_Q1_deng_q975,
                   mean_bias_Q2_deng, mean_bias_Q2_deng_q25, mean_bias_Q2_deng_q975,
                   mean_bias_Q3_deng, mean_bias_Q3_deng_q25, mean_bias_Q3_deng_q975,
                   mean_bias_Q4_deng, mean_bias_Q4_deng_q25, mean_bias_Q4_deng_q975,
                   mean_bias_Q5_deng, mean_bias_Q5_deng_q25, mean_bias_Q5_deng_q975,
                   mean_bias_Q6_deng, mean_bias_Q6_deng_q25, mean_bias_Q6_deng_q975,
                   
                   mean_pi_qin, mean_pi_qin_lowerCI, mean_pi_qin_upperCI,
                   mean_bias_Q1_qin, mean_bias_Q1_qin_q25, mean_bias_Q1_qin_q975,
                   mean_bias_Q2_qin, mean_bias_Q2_qin_q25, mean_bias_Q2_qin_q975,
                   mean_bias_Q3_qin, mean_bias_Q3_qin_q25, mean_bias_Q3_qin_q975,
                   mean_bias_Q4_qin, mean_bias_Q4_qin_q25, mean_bias_Q4_qin_q975,
                   mean_bias_Q5_qin, mean_bias_Q5_qin_q25, mean_bias_Q5_qin_q975,
                   mean_bias_Q6_qin, mean_bias_Q6_qin_q25, mean_bias_Q6_qin_q975
  )
  
  vec_results <- c(scenario, vec_results)
  
  names(vec_results) <- c(colnames(scenario),
    
                          "mean_pi_deng", "mean_pi_deng_lowerCI", "mean_pi_deng_upperCI",
                          "mean_bias_Q1_deng", "mean_bias_Q1_deng_q25", "mean_bias_Q1_deng_q975",
                          "mean_bias_Q2_deng", "mean_bias_Q2_deng_q25", "mean_bias_Q2_deng_q975",
                          "mean_bias_Q3_deng", "mean_bias_Q3_deng_q25", "mean_bias_Q3_deng_q975",
                          "mean_bias_Q4_deng", "mean_bias_Q4_deng_q25", "mean_bias_Q4_deng_q975",
                          "mean_bias_Q5_deng", "mean_bias_Q5_deng_q25", "mean_bias_Q5_deng_q975",
                          "mean_bias_Q6_deng", "mean_bias_Q6_deng_q25", "mean_bias_Q6_deng_q975",
                          
                          "mean_pi_qin", "mean_pi_qin_lowerCI", "mean_pi_qin_upperCI",
                          "mean_bias_Q1_qin", "mean_bias_Q1_qin_q25", "mean_bias_Q1_qin_q975",
                          "mean_bias_Q2_qin", "mean_bias_Q2_qin_q25", "mean_bias_Q2_qin_q975",
                          "mean_bias_Q3_qin", "mean_bias_Q3_qin_q25", "mean_bias_Q3_qin_q975",
                          "mean_bias_Q4_qin", "mean_bias_Q4_qin_q25", "mean_bias_Q4_qin_q975",
                          "mean_bias_Q5_qin", "mean_bias_Q5_qin_q25", "mean_bias_Q5_qin_q975",
                          "mean_bias_Q6_qin", "mean_bias_Q6_qin_q25", "mean_bias_Q6_qin_q975")
  
  return(as.data.frame(vec_results))
}

    
