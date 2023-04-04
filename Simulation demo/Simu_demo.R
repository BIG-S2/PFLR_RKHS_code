source('Gen_data.R')
source('PFLR_RKHS_estimate.R')



n <- 200
pn <- 150
nu <- 2
rho1 <- 0.2
rho2 <- 0.3

eigen_center <- 1
sd_noise <-1
kp=1 ## for Gaussian kernel

k <- 50
lent <- 100
locations <- seq(0, 1, length=lent)

B=4

data <- Gen_data(n, pn,  k, locations, nu, rho1,rho2, eigen_center, sd_noise, B)
y <- data$y
x <- data$x
z <- data$z

type='Gaussian'
K_gram <- outer( locations, locations, function(x, y) exp(-(x-y)^2/kp) )


nbasis <- n

lambda_range.1 <-  seq(1e-5, 0.1, length=50) 
sp_range <- c(1:50)
init_gamma <- rep(0, pn)
max_iter <- 200

    
Sigma_mat <- x%*%K_gram%*%t(x)/(lent^2)
    
    
fit_proposed <-PFLR_RKHS_estimate(y, x, z, K_gram, Sigma_mat,locations,n_basis, lambda_range.1, sp_range, init_gamma, max_iter)






