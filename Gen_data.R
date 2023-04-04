library(MASS)
Gen_data <- function(n, pn,  k, locations, nu, rho1,rho2, eigen_center, sd_noise){
  # n: sample size
  # pn: number of scalar variables
  # k: number of basis
  # locations: observation points of the functional variable
  # nu: smoothness of the functional variable
  # rho1: correlation of the FPC scores and the scalar variables
  # rho2: correlation of the scalar variables
  # eigen_center: to generate disordered representative basis
  # sd_noise: sd of the error
  lent <- length(locations)
  
  ######################################
  ## xi1-xi4 and  z ####################
  ######################################
  pn4=pn+4;   ## pn is the number of scalar variates
  sigma=matrix(0,pn4,pn4)  ## generate a zeros matrix
  sigma[1,1]=16; sigma[2,2]=16*2^(-nu); sigma[3,3]=16*3^(-nu); sigma[4,4]=16*4^(-nu);
  
  ## the four leading FPC scores correlate with z
  for(i in 1:4){
    for(j in 5:pn4){
      sigma[i,j]=rho1^(abs(i-j+4)+1)*sqrt(sigma[i,i]);  
      sigma[j,i]=sigma[i,j];
    }
  }
  
  for(i in 5:pn4){
    for(j in 5:pn4){
      sigma[i,j]=rho2^(abs(i-j));} 
  }
  
  xi_z=mvrnorm(n=n,rep(0,pn4),sigma)  
  xi_14=xi_z[,1:4]; z=xi_z[,5:pn4] 
  
  #####################################
  ## generate  xi and x ###############
  #####################################
  phi=function(k,t){
    if(k%%2==0)  return(sqrt(2)*sin((k-1)*pi*t))  else  return(sqrt(2)*cos(k*pi*t));
  }
  
  
  xi=apply(matrix(1:k,1,k), 2 , function(x)rnorm( n,0,4*( abs( (x-eigen_center))+1 )^(-nu/2) ))   ### n*k
  
  xi[, eigen_center:(eigen_center+3)]=xi_14
  
  j_trans=as.vector( matrix(rep(c(1:k),lent),lent,k,byrow=TRUE) )
  basis= matrix( mapply( function(x, t){ phi(x,t)}, j_trans, rep(locations,k) ), k, lent, byrow=TRUE)
  
  x <- xi%*%basis  # n*t
  
  ##### generate beta ###########
  beta_coe <- apply(as.matrix(c(1:k)), 1, function(x){ 4*(-1)^(x+1)/x^2} )
  
  beta <- as.numeric( beta_coe%*%basis  ) 
  
  
  ### generate response #####
  gamma <- c(3, 1.5, 1, 2.5, 2, rep(0, pn-5) )
  
  y <- as.numeric( x%*%beta )/lent+as.numeric( z%*%gamma ) + rnorm(n, 0, sd_noise)
  
  K_gram <- t(basis)%*%basis
  
  
  list(y=y, x=x, z=z, beta_true=beta, gamma=gamma, K_gram=K_gram, idx_true = c(rep(1,5), rep(0, pn-5)))
  
}