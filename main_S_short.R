args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])

set.seed(seed)
source("functions_S.R")
usePackage("splines2")
usePackage("TruncatedNormal")
usePackage("mvtnorm")
usePackage("matrixStats")
source("exchangable_cov.R")
options(warn=0)

library(EnvStats)

start_time_low <- 50
start_time_upper <- 90 # Normal
interval_time_min <- 1
interval_time_exp_rate <- 20 # 1+Exp(rate)
num_visits_mean <- 10 # Poisson
N_cont_covars <- 2 # N(0,1) continous covariates
N_binary_covars <- 2
p_binary_covars <- c(1,0.5)
random_effect_var <- 1
residual_var <- 0.5

true_fixed_effect <- matrix(c(+0.4,+0.4,-0.5,
                              -0.5,-0.5,+0.1,
                              +0.1,+0.1,+0.5),nrow = 3,ncol = 3,byrow=TRUE)
nX <- 3
# true_fixed_effect <- matrix(c(0,0,0),1,3)
a0 <- 1; b0 <- 70; d0 <- 5
mode1<-50; range_L1 <- 30; range_R1 <- 100
mean1 <- 70; sd1 <- 5; mean2 <- 100; sd2 <- 5; p1 <- 0.4; p2 <- 0.6

N <- 250

dataset_num <- 1

CI_repeat <- array(0,dim=c(dataset_num,1201,8))
turning <- array(0,dim=c(dataset_num,6))
CI_covariate_repeat <- array(0,dim=c(dataset_num,nrow(true_fixed_effect),7))
true_turning <- rep(0, dataset_num)
Q50 <- array(0,dim=c(dataset_num,6))
true_Q50 <- rep(0, dataset_num)



#coef_repeat_S <- array(0,dim=c(dataset_num,5000,4+4))

RE_repeat <- array(0,dim=c(dataset_num,N,7))
# RE + fixed intercept
offset_repeat <- array(0,dim=c(dataset_num,N,7))

sigmay_repeat <- array(0,dim=c(dataset_num,7))
sigmaw_repeat <- array(0,dim=c(dataset_num,7))

for(di in 1:dataset_num){
  
cat(di);cat(":\n")
  # Generate a set of data

start_times <- runif(N, start_time_low, start_time_upper)
visits <- pmax(rpois(N, num_visits_mean),1)
df <- data.frame()
for(i in 1:N){
  interval <- rexp(visits[i]-1,interval_time_exp_rate) + 1
  df <- rbind(df, data.frame(
    ageori = start_times[i] + c(0, cumsum(interval)),
    id = rep(i, visits[i])
  ))
}
df <- df[dplyr::between(df$ageori,0,120),]
dfi <- 8
qknot <- (1:(dfi-3))/(dfi-2)
VIF <- 0.1
boundary.knot <- c(0,120)
t01 <- (df$ageori-boundary.knot[1])/(boundary.knot[2]-boundary.knot[1])
t01 <- t01[dplyr::between(t01,0,1)]
knot <- betaKDE(t01,s=VIF,q=qknot)$quantile
knot <- (boundary.knot[2]-boundary.knot[1])*knot+boundary.knot[1]

X <- cbind(X1=rbinom(N,1,p_binary_covars[1]),
           X2=rbinom(N,1,p_binary_covars[2]),
           X3=rnorm(N))
X_names <- colnames(X)
# X <- matrix(0, nrow(X), ncol(X))
# X <- matrix(rep(1,N),N,1)
# colnames(X) <- c('intercept')
df <- cbind(df, X[df$id,])

Y <- as.matrix(df[,X_names]) %*% true_fixed_effect
truthRE <- matrix(rnorm(N*3,sd=sqrt(random_effect_var)),nrow=N,ncol=ncol(Y))
Y <- Y + matrix(rnorm(length(Y),sd=sqrt(residual_var)),nrow=nrow(Y),ncol=ncol(Y))
Y <- Y + truthRE[df$id,]

coef00 <- c(0,0,c(0.1,5,8,0.1)/100,0,0)
B00 <- splines2::ibs(df$ageori,knots=knot,degree=2,intercept=TRUE,Boundary.knots=c(0,120))

Y[,1] <- Y[,1] + B00 %*% coef00
Y[,2] <- Y[,2] + f_sshape(df$ageori,mode1,range_L1,range_R1)
Y[,3] <- Y[,3] + f_wiggle(df$ageori,mean1,sd1,mean2,sd2,p1,p2)
colnames(Y) <- c('Y1','Y2','Y3')

mis <- missing_pattern(nrow(Y),3,1/5,1/10)
for(i in 1:nrow(Y))
  for(j in 1:ncol(Y))
    if(mis[i,j]) Y[i,j] <- NA

Y <- Y[,1]
mis <- mis[,1]
truthRE0 <- truthRE[df$id,1]
df <- cbind(df,Y,truthRE0)
# Loading Data ---------------------------------------------------
df <- df

K <- 1 # Number of biomarkers
# All non-age covariates
X <- as.matrix(df[,X_names],ncol=nX) 
Y <- as.matrix(df[,c('Y')],ncol=1) # Biomarkers array
t <- df$ageori # Age in original scale
dfi <- 8 # DoF of Spline
qknot <- (1:(dfi-3))/(dfi-2) # Quantiles to determine knots
VIF <- 0.1 # Variance inflation factor for BETAKDE

library(ggplot2)

# Construct biomarker-specific design matrix
covar.list <- as.list(rep(NA,K))
knot.list <- as.list(rep(NA,K))
boundary.knot <- c(30,90)#range(t)
#remove <- 2# Removing the first .. and last .. basis
for(i in 1:K){
  # Calculate knot points
  t01 <- (t-boundary.knot[1])/(boundary.knot[2]-boundary.knot[1])
  t01 <- t01[dplyr::between(t01,0,1)]
  knot.list[[i]] <- betaKDE(t01,s=VIF,q=qknot)$quantile
  knot.list[[i]] <- (boundary.knot[2]-boundary.knot[1])*knot.list[[i]]+boundary.knot[1]
  B <- ibs(pmin(pmax(t,min(boundary.knot)),max(boundary.knot)), 
           knots=knot.list[[i]], Boundary.knots = boundary.knot, 
           degree=2, intercept=TRUE) # IBSpline Basis
  #B <- iSpline(t, knots=knot.list[[i]], Boundary.knots = boundary.knot, 
  #             degree=2, intercept=TRUE) # IBSpline Basis
  #B <- B/min(apply(B,2,max))
  B <- B[,(3):(ncol(B)-2)]
  covar.list[[i]] <- cbind(X,B)
}


# Create consecutive pseudo-IDs for each individual for easy coding
unique.IDs <- sort(unique(df$id))
df$ID <- match(df$id,unique.IDs)
# Pre-calculate longitudinal sample size 
# for each individual-biomarker combination
long_ss <- matrix(0,nrow=length(unique.IDs),ncol=K)
for(i in 1:length(unique.IDs))
  for(j in 1:K){
    long_ss[i,j] <- sum((df$ID==i)&(!is.na(Y[,j])))
  }

long_all_ss <- rep(0,length(unique.IDs))
for(i in seq_along(unique.IDs)){
  long_all_ss[i] <- sum(df$ID==i)
}

R <- 1e4 # Set Number of Iterations
Burnin <- R/2 # Set Number of Burn-ins



# Set Priors -----------------------------------------------------
# Beta parameter: Coefficients for adjusting covariates
beta.prior <- list(mean=rep(0,ncol(X)),
                   variance=diag(rep(10000,ncol(X))),
                   #variance=10000,
                   precision=NULL)
beta.prior$precision <- solve(beta.prior$variance)
# Gamma parameter: Coefficients for splines
gamma.prior <- list(mean=rep(0,ncol(B)),
                    variance=NULL)
coef.prior <- list(mean=c(beta.prior$mean, gamma.prior$mean),
                   variance=NULL,
                   precision=NULL)
#coef.prior$precision <- solve(coef.prior$variance)


# Set initial guess ----------------------------------------------

# Fixed Effect of X & All-positive Spline Coefs
coefs <- array(0,c(ncol(covar.list[[1]]),ncol(Y),R)) 
nX <- ncol(X)
#inflex_prob <- array(0,c(ncol(covar.list[[1]])-nX,ncol(Y),R))
# Variance of Biomarkers
sigmays <- rep(0,R)
sigmaws <- rep(0,R)
pens <- array(0, c(2, R))
REs <- array(0, c(dim(long_ss),R))
offsets <- array(0, c(dim(long_ss),R))

coefs[1:nX,,1] <- 
  t(rtmvnorm(ncol(Y),mu=beta.prior$mean,
             sigma=beta.prior$variance))
coefs[(nX+1):ncol(covar.list[[1]]),,1] <- 
  t(rtmvnorm(ncol(Y),mu=gamma.prior$mean,
             sigma=penalty_Matrix(ncol(B),
                                  smooth.sigma = 1,flat.sigma = 1)$V,
             lb=rep(0,ncol(B))))
sigmays[] <- residual_var#1/rgamma(1,shape=3,rate=0.5)
pens[,1] <- c(0.1,0.5)
sigmaws[] <- random_effect_var#1/rgamma(1,shape=3,rate=0.5)
# REs[,,1] <- array(rnorm(prod(dim(long_ss))),
#                   dim=dim(long_ss))

# Prior density for penalties
lpd <- function(s) 
  log(2)+dnorm(s[1], log=TRUE) +
  log(2)+dnorm(s[2], log=TRUE)
ls <- -2 # Log of Jump Standard Deviation
acc <- 0 # Accepted Proposals in one batch
lss <- ls # Sequence of LS for reference
w <- planck_taper(ncol(B), eps=0.5) # Window Function
w <- rep(1,ncol(B))
w <- NULL

M_coef <- lincon(nX+dfi-4,nX)
M_pen <- lincon(dfi-4,0)

# Perform MCMC ----------------------------------------------------
# i <- 1
for(i in 1:(R-1)){
 if((i+1)%%(R/10)==0) cat(sprintf("%03d%% ",(i+1)/(R/100)))
 if((i+1)%%(R/1)==0) cat("\n")
  verbose=FALSE
  prec <- block_Matrix(beta.prior$precision,
                       penalty_Matrix(ncol(B),pens[1,i],pens[2,i],
                                      weight=w)$prec)
  sigmays[i+1] <- update_sigmay(covar.list,Y,as.matrix(REs[df$ID,,i],ncol=1),
  as.matrix(coefs[,,i],ncol=1),
  3,0.5)
  u <- update_coef(covar.list,nX,Y,as.matrix(REs[df$ID,,i],ncol=1),
                   sigmays[i+1],sigmaws[i],df$ID,
                   coef.prior$mean,
                   prec,M_coef,verbose,samples=1)
  coefs[,,i+1] <- aperm(u$res,c(2,3,1))
  REs[,,i+1] <- update_W(covar.list,Y,as.matrix(coefs[,,i+1],ncol=K),long_ss,
                         df$ID,sigmays[i+1],sigmaws[i])
  #inflex_prob[,,i+1] <- u$inflex_prob
  # current <- "Penalty"
  new_pens <- update_pens(gamma=as.matrix(coefs[(nX+1):ncol(covar.list[[1]]),,i+1],ncol=1),
                          mu=gamma.prior$mean,
                          lambda=pens[,i],
                          lpd=lpd,
                          ls=ls,
                          weight=w,
                          Ms=M_pen,
                          verbose=verbose
  )
  pens[,i+1] <- new_pens$new
  acc <- acc + new_pens$acc_status
  if(i%%50==0){
    delta <- min(0.1, 1/sqrt(i/50))
    rate <- acc/50
    if(rate >= 0.234){
      ls <- ls + delta
    }
    else{
      ls <- ls - delta
    }
    acc <- 0
    lss <- c(lss, ls)
  }
  # pens[,i+1] <- pens[,i]
  sigmaws[i+1] <- update_sigmaw(REs[,,i+1],3,0.5)
}

for(i in 1:R){
  for(k in 1:K){
  offsets[,k,i] <- REs[,k,i] + coefs[1,k,i]
  }
}

ages <- seq(0,120,by=0.1)
points <- array(0,c(length(ages),R-Burnin))

indice <- seq(Burnin+1,R,1)
spline.basis <- splines2::ibs(pmin(pmax(ages,min(boundary.knot)),max(boundary.knot)), 
                              knots=knot.list[[1]], Boundary.knots = boundary.knot, 
                              degree=2, intercept=TRUE)
spline.basis <- spline.basis[,3:(dfi-2)]
points <- spline.basis %*% coefs[-(1:nX),1,indice]
est <- apply(points,1,function(x) c(mean(x),
                                    coda::HPDinterval(coda::as.mcmc(x))))
var_est <- apply(points,1,var)
est <- data.frame(t(est))
colnames(est) <- c("avg","lower","upper")
est$truth <- spline.basis %*% coef00[3:6]
est$age <- ages
est$MSE <- (est$avg - est$truth)^2+var_est
est$bias2 <- (est$avg - est$truth)^2
est$var <- var_est

inflects <- apply(points, 2, function(x){
  ages[max(which(diff(x,differences=2)>0))+1]
})

turning[di,1:2] <- HDInterval::hdi(inflects)
turning[di,3] <- mean(inflects)
true_turning[di] <- ages[max(which(diff(est$truth,differences=2)>0))+1]
turning[di,5] <- (turning[di,3] - true_turning[di])^2
turning[di,6] <- var(inflects)
turning[di,4] <- sum(turning[di,5:6])

Q50s <- apply(points, 2, function(x){
  ages[min(which(x>=max(x)/2))]
})
Q50[di,1:2] <- HDInterval::hdi(Q50s)
Q50[di,3] <- mean(Q50s)
true_Q50[di] <- ages[min(which(est$truth>=max(est$truth)/2))]
Q50[di,5] <- (Q50[di,3] - true_Q50[di])^2
Q50[di,6] <- var(Q50s)
Q50[di,4] <- sum(Q50[di,5:6])

CI_repeat[di,,] <- as.matrix(est)

CI_covariate_repeat[di,,1:3] <- t(apply(coefs[1:nX,1,indice],1,
                                        function(x) c(mean(x),coda::HPDinterval(coda::as.mcmc(x)))))
CI_covariate_repeat[di,,4] <- c(0.4,-0.5,0.1)
CI_covariate_repeat[di,,7] <- t(apply(coefs[1:nX,1,indice],1,var))
CI_covariate_repeat[di,,6] <- (CI_covariate_repeat[di,,1]-CI_covariate_repeat[di,,4])^2
CI_covariate_repeat[di,,5] <- CI_covariate_repeat[di,,6] + CI_covariate_repeat[di,,7]



#coef_repeat_S[di,,] <- t(coefs[,1,indice])
RE_repeat[di,,1:3] <- t(apply(REs,c(1,2),
                              function(x) c(mean(x),coda::HPDinterval(coda::as.mcmc(x))))[,,1])
RE_repeat[di,,4] <- truthRE[,1]
RE_repeat[di,,7] <- t(apply(REs[,,indice],1,var))
RE_repeat[di,,6] <- (RE_repeat[di,,1]-RE_repeat[di,,4])^2
RE_repeat[di,,5] <- RE_repeat[di,,6] + RE_repeat[di,,7]

offset_repeat[di,,1:3] <- t(apply(offsets,c(1,2),
                              function(x) c(mean(x),coda::HPDinterval(coda::as.mcmc(x))))[,,1])
offset_repeat[di,,4] <- truthRE[,1] + 0.4
offset_repeat[di,,7] <- t(apply(offsets[,,indice],1,var))
offset_repeat[di,,6] <- (offset_repeat[di,,1]-offset_repeat[di,,4])^2
offset_repeat[di,,5] <- offset_repeat[di,,6] + offset_repeat[di,,7]

sigmay_repeat[di,1] <- mean(sigmays[indice])
sigmay_repeat[di,2:3] <- coda::HPDinterval(coda::as.mcmc(sigmays[indice]))
sigmay_repeat[di,4] <- residual_var
sigmay_repeat[di,6:7] <- c((sigmay_repeat[di,1]-sigmay_repeat[di,4])^2,
                           var(sigmays[indice]))
sigmay_repeat[di,5] <- sum(sigmay_repeat[di,6:7])

sigmaw_repeat[di,1] <- mean(sigmaws[indice])
sigmaw_repeat[di,2:3] <- coda::HPDinterval(coda::as.mcmc(sigmaws[indice]))
sigmaw_repeat[di,4] <- random_effect_var
sigmaw_repeat[di,6:7] <- c((sigmaw_repeat[di,1]-sigmaw_repeat[di,4])^2,
                           var(sigmaws[indice]))
sigmaw_repeat[di,5] <- sum(sigmaw_repeat[di,6:7])

}
save(CI_repeat,turning,true_turning,CI_covariate_repeat,RE_repeat,offset_repeat,
     Q50, true_Q50,
     sigmay_repeat,sigmaw_repeat,file=sprintf('sS_CIs_%03d.rda',seed))
covered <- apply(CI_repeat,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate <- apply(covered, 2, mean)

