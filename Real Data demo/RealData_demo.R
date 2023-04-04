

######## Images ########
final_image <- as.matrix( read.table(('./Data/Hippocampus_left.txt'), header=F) )
gene_data <- as.matrix( read.csv('./Data/SNP_select_combine.csv', header=T) )[, -1]


##############
source('PFLR_RHKS_image.R')

dim_1 <- 100
dim_2 <- 150

dim <- dim_1*dim_2

locations_1 <- seq(0,1, length=dim_1)
locations_2 <- seq(0,1, length=dim_2)

locations <- cbind( rep(locations_1, dim_2), rep(locations_2, each=dim_1) )


#### dim*dim matrix####

kp <- 1

K_gram <- outer( c(1:dim), c(1:dim), function(x, y) exp(-( (locations[x,1]-locations[y,1])^2 +(locations[x,2]-locations[y,2])^2  )/kp)  )

###### extract variables that need to be regressed out #########

D <- as.matrix( read.csv('./Data/final_clinical_m12_used.csv', header=T) )

n <- dim(D)[1] # sample size

D[, c(4,6,10:15)] <- scale(D[, c(4,6,10:15)], center = T, scale = T)  # scale the continuous variable

y <- D[, 15]  # MMSE score

D <- D[, -c(1, 7:9,15)]  # exclude APOE4 and diagnosis

inv_DD <- solve( crossprod(D, D), diag( dim(D)[2] ) )
Projection_D <- diag(n) - tcrossprod( tcrossprod(D, inv_DD), D)



###############################
## tuning parameters settings #
###############################

lent <-dim


nbasis <- n       # the basis for the RHKS estimates of the functional coefficient, default is n.

###########

lambda_range.1 <-  seq(0.01, 0.02, length=5) ## penalty tuning for the functional coefficient

sp_range2 <- c(50:100)

max_iter <- 50

lambda_ridge <- 1e-4



##################################################
## add clinical covariates and run the model centeralization #####
##################################################
x_nonProject <- final_image
x <- Projection_D%*%x_nonProject

x_center <- scale(x, center = T, scale = F)


########## Create the Sigma matrix ########
Sigma_mat <- x_center%*%K_gram%*%t(x_center)/(dim^2)

###### Create Kx_gram #####################
Kx_gram <- x_center %*%K_gram/dim  


###############################
## tuning parameters settings #
###############################

lent <-dim


nbasis <- n       # the basis for the RHKS estimates of the functional coefficient, default is n.
  
###################
y_nonProject <- y
z_nonProject <- scale(gene_data, center=T, scale=T) 


y <- Projection_D%*%y_nonProject
z <- Projection_D%*%z_nonProject

y_center <- y - mean(y)

z_center <- scale(z, center = T, scale = F)

pn <- dim(z)[2]

init_gamma <- rep(0, pn)   ### Initial value for gamma

###################
fit_proposed <- PFLR_RKHS_image(y_center, x_center, z_center, K_gram, Kx_gram, Sigma_mat,n_basis, 
                                 lambda_range.1, sp_range2, init_gamma, max_iter,dim_1,dim_2,lambda_ridge)


