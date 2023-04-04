source('PFLR_RKHS_estimate.R')   ### the proposed method
source('PFLR_Yao_estimate.R')    ### the method proposed by Kong et al. (2016)

####################
### Load data ######
####################
load("Demo_data.Rdata")

# y: the scalar response
# x: the functional predictior, n*lent with n denotes the sample size and lent denotes the number of observations in the functional domain
# z: the high-dimensional scalar covariates, n*pn with pn denotes the number of scalar covariates
# Model
###################
# y = z^\top \gamma + \int x(t) \beta(t) +\epsilon(t)
#

n <- length(y) # sample size
pn <- dim(z)[2] 
lent <- dim(x)[2]
locations <- seq(0, 1, length=lent)


###############################
## tuning parameters settings #
###############################
nbasis <- n                       # the basis for the RHKS estimates of the functional coefficient, default is n.

kp=1
K_gram <- outer( locations, locations, function(x, y) exp(-(x-y)^2/kp) )  ## Gaussian Kernel

Sigma_mat <- x%*%K_gram%*%t(x)/(lent^2)

###########

lambda_range.1 <-  seq(0.01, 0.02, length=5) ## penalty tuning for the functional coefficient
sp_range <- c(3,4,5,6,7)                     ## sparsity level for the scalar coefficients

init_gamma <- rep(0, pn)   ### Initial value for gamma
max_iter <- 200



sn_range <- c(2,3,4,5,6)   ## number of Fpcs for Yao's method
lambda_range.2 <- seq(0.2, 0.35, length=5)  ## penalty tuning for the scad penalty of Yao's method

##############
t0 <- Sys.time()
fit_proposed <-PFLR_RKHS_estimate(y, x, z, K_gram, Sigma_mat,locations,n_basis, lambda_range.1, sp_range, init_gamma, max_iter)
t_proposed <- difftime(Sys.time(), t0, units='secs')

#fit_proposed$final_results: results for the selected tuning
#fit_proposed$final_results$gamma_est: estimate for gamma
#fit_proposed$final_results$beta_est: esetimate for beta
#fit_proposed$final_results$idx: True denotes that the scalar is selected
#fit_proposed$results.sp.all: resutls for all the sparsity level (sp_range)


t1 <- Sys.time()
fit_Yao <- PFLR_Yao_estimate(y, x, z, locations, sn_range, lambda_range.2)
t_Yao <- difftime(Sys.time(), t1, units='secs')

#fit_Yao$final_results: results for the selected tuning
#fit_Yao$final_results$gamma_est: estimate for gamma
#fit_Yao$final_results$beta_est: esetimate for beta
#fit_Yao$final_results$idx: True denotes that the scalar is selected
#fit_Yao$results.sp.all: resutls for all the lambda_range.2
