library(refund)
library(ncvreg)
PFLR_Yao_estimate <- function(y, x, z, locations, sn_range, lambda_range){
  
  lent <- dim(x)[2]
  pn <- dim(z)[2]
  
  sn_length <- length(sn_range)
  lambda_length <- length(lambda_range)
  
  BIC <- matrix(nrow=sn_length, ncol=lambda_length)
  
  fpca_x <- fpca.sc(x, npc=max(sn_range))
  
  scores <- fpca_x$scores
  efunctions <- fpca_x$efunctions
  
  results.all <- NULL
  
  for(i in 1:sn_length){
    sn.tmp <- sn_range[i]
    
    scores.tmp <- scores[, 1:sn.tmp]
    efunctions.tmp <- efunctions[, 1:sn.tmp]
    
    for(j in 1:lambda_length){
      lambda.tmp <- lambda_range[j]
      
      scalar_all<- cbind(scores.tmp, z)
      
      fit <- ncvfit(scalar_all, y, penalty='SCAD', lambda=lambda.tmp, penalty.factor = c( rep(0, sn.tmp), rep(1, pn) ))
      
      coe_beta_est <- fit$beta[1:sn.tmp]
      
      beta_est <- efunctions.tmp%*%coe_beta_est
      
      gamma_est <- fit$beta[-c(1:sn.tmp)]
      
      idx_est <- 1-(gamma_est==0)
      
      loss <- fit$loss
      
      BIC[i, j]=log(loss) + log(n)*(sn.tmp+sum(idx_est))
      
      results.tmp <- list(gamma_est= gamma_est,
                           beta_est=beta_est, idx=idx_est, sn=sn.tmp,lambda =lambda.tmp, BIC=BIC[i,j] )
      
      results.all <- c(results.all, list(results.tmp))
      
    }
  }
  tuning_position <- which(BIC==min(BIC), arr.ind=T)
  results_position <- as.numeric(which(t(BIC)==min(BIC)))
  
  if(length(results_position) >1){
    tuning_position <- tuning_position[1, ]
    results_position <- results_position[1]
    }
 
  
  final_sn <- sn_range[tuning_position[1]]
  final_lambda <- lambda_range[tuning_position [2]]
  
  final_results <- results.all[[results_position]]
  
  list(final_results=final_results, results.all=results.all)
  
}