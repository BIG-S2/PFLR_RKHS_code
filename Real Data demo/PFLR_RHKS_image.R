PFLR_RKHS_image <- function(y, x, z, K_gram, Kx_gram, Sigma_mat,n_basis, lambda_range, sp_range, init_gamma, max_iter,dim_1,dim_2,lambda_ridge){
  # K_gram: matrix of the RK
  # Sigma_mat: matrix of Sigma
  # n_basis: number of basis for the coefficient c
  # lambda: tuning parameter wrt the penalty of functional variable
  # sp_range: sparsity level
  # init_gamma: initial value of gamma 
  # max_iter: max iteration
  n <- length(y)
  
  lent <- dim(x)[2]
  
  sp_length <- length(sp_range)
  lambda_length <- length(lambda_range)
  
  results.sp <- NULL
  HBIC.sp <- NULL
  
  HBIC3.sp <- NULL
  
  for(i in 1:sp_length){

    t_0 <- Sys.time()
    
    sp_tmp <- sp_range[i]
    
    GCV.sp <- NULL
    gamma_est.sp <- NULL
    coe_beta_est.sp <- NULL
    idx.sp <- NULL
    
    for(j in 1:lambda_length){
      
      lambda_tmp <- lambda_range[j]
      
      P <- solve( Sigma_mat/(n*lambda_tmp)+diag( n  ), diag(n)  )
      
      d <- t(z)%*%P%*%( y-z%*%init_gamma )/n
      
      active_idx <- NULL
      iter <- 0
      
      update_gamma <- init_gamma
      
      while( iter <max_iter ){
        iter <- iter+1
        pre_active_idx <- active_idx
        arr <- abs( update_gamma+d )
        thre <- sort(arr, decreasing = T)[sp_tmp]
        
        active_idx <-  which( arr >= thre ) 
        inactive_idx <- which( arr < thre ) 
        
        idx <- arr >= thre
        
        active_z <- as.matrix(z[, active_idx])
        inactive_z <- z[, inactive_idx]
        
        Inv <- solve( t(active_z)%*%P%*%active_z  + lambda_ridge * diag(  dim(active_z)[2] ) )
        
        update_gamma[inactive_idx] <- 0
        update_gamma[active_idx] <- Inv%*%t(active_z)%*%P%*%y
        
        d[active_idx] <- 0
        d[inactive_idx] <- t(inactive_z)%*%P%*%( y-active_z%*%update_gamma[active_idx] )/n
        
        if( identical(active_idx, pre_active_idx)==T  ){break}
      }
            
      res_tmp <- y - z%*%update_gamma
      coe_beta_est <- P%*%res_tmp/(n*lambda_tmp)
      
      #beta_est <- as.numeric( t( coe_beta_est)%*%Kx_gram )
      
      gcv_tmp <- sum( (P%*%res_tmp)^2)/( sum(diag(P)) )^2  #  omit n
      
      GCV.sp <- c(GCV.sp, gcv_tmp)
      
      gamma_est.sp <- cbind( gamma_est.sp, update_gamma )
      coe_beta_est.sp <- cbind( coe_beta_est.sp, coe_beta_est )
      idx.sp <-cbind(idx.sp, idx)
    }
    
    lambda_position<- which( GCV.sp==min(GCV.sp) )
    gamma_est_tmp <- gamma_est.sp [, lambda_position]
    coe_beta_est_tmp <- coe_beta_est.sp[, lambda_position]
    idx_tmp <- idx.sp [, lambda_position]
    
    #residual <- y - x%*%beta_est_tmp/lent - z%*%gamma_est_tmp
    
    residual <- y - Sigma_mat%*%coe_beta_est_tmp - z%*%gamma_est_tmp
    
    #HBIC.tmp <-  log(sum(residual^2) ) + sp_tmp*log(log(n))*log(pn)/n
    HBIC.tmp <-  log(sum(residual^2) ) + sp_tmp*log(log(n))/n
        
    HBIC3.tmp <-  log(sum(residual^2) ) + 5*sp_tmp*log(log(n))/n
    
    results.tmp <- list(  gamma_est= gamma_est_tmp, coe_beta_est=coe_beta_est_tmp, idx=idx_tmp, sp=sp_tmp,
                          lambda =lambda_range[lambda_position], iter=iter, hbic=HBIC.tmp, hbic3=HBIC3.tmp )
    HBIC.sp <- c( HBIC.sp, HBIC.tmp )
    HBIC3.sp <- c( HBIC3.sp, HBIC3.tmp )
    
    results.sp <- c(results.sp, list(results.tmp))
    
    print(Sys.time()-t_0)
    
  }
  
  
  sp_position <- which( HBIC.sp==min(HBIC.sp))
  sp3_position <- which( HBIC3.sp==min(HBIC3.sp))
  
  final_results <- results.sp[[sp_position]]
  final_results_3 <- results.sp[[sp3_position]]
  
  final_results$beta_est <- matrix(as.numeric( t( final_results$coe_beta_est)%*%Kx_gram), dim_1, dim_2)
  
  list( final_results=final_results, final_results_3=final_results_3, results.sp.all=results.sp )
}