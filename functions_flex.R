# Functions
# Create PT window
planck_taper <- function(dimension, eps=0.1) {
  N <- dimension-1
  w <- rep(0, dimension)
  w[1] <- 0
  if(eps*N>=1)
    for(i in seq(1,eps*N,by=1)){
      w[i+1] <- (1+exp(eps*N/i-eps*N/(eps*N-i)))^(-1)
    }
  for(i in seq(ceiling(eps*N),N/2)) w[i+1] <- 1
  for(i in seq(0, N/2)) w[N-i+1] <- w[i+1]
  return(w)
}

# Create random-walk + window prior matrix
penalty_Matrix <- function(dimension,smooth.sigma=Inf,flat.sigma=Inf,weight=NULL)
{
  if(is.null(weight)){
    weight <- (0:(dimension-1))*((dimension-1):0)
  }
  weight <- weight + 0.001
  weight <- weight/max(weight)
  prec <- matrix(0,dimension,dimension)
  for(j in 1:(dimension-1)){
    u <- rep(0,dimension); u[j] <- 1; u[j+1] <- -1;
    prec <- prec + u%*%t(u)/(smooth.sigma/dimension/dimension)
  }
  prec <- prec + diag(1/flat.sigma/weight)
  V <- solve(prec)
  return(list(V=V,prec=prec))
}

# Create block diagonal matrices
block_Matrix <- function(M1,M2){
  rbind(cbind(M1,matrix(0,nrow=nrow(M1),ncol=ncol(M2))),
        cbind(matrix(0,nrow=nrow(M2),ncol=ncol(M1)),M2))
}

betaKDE <- function(z, s, q){
  # If s is NULL, use empirical quantiles
  # If s is between 0 and 1, use BETA-KDE
  # Otherwise, terminate the program with error
  if(is.null(s)){
    quan <- quantile(z, q)
    names(quan)=paste(round(q*100,1),"%",sep='')
    return(list(density=NULL,CDF=NULL,quantile=quan))
  }
  if((s<=0)||(s>=1)) 
    stop("Variance inflation factor should be NULL or between 0 and 1")
  a <- 1/s
  dens <- function(x){
    res <- 0
    for(t in z){
      res <- res + dbeta(x,a*t+1,a*(1-t)+1)
    }
    res <- res/length(z)
    return(res)
  }
  CDF <- function(x){
    res <- 0
    for(t in z){
      res <- res + pbeta(x,a*t+1,a*(1-t)+1)
    }
    res <- res/length(z)
    return(res)
  }
  
  quan <- rep(0, length(q))
  for(i in 1:length(q)){
    quan[i] <- uniroot(function(x){CDF(x)-q[i]},
                       c(0,1))$root
  }
  names(quan)=paste(round(q*100,1),"%",sep='')
  return(list(density=dens,CDF=CDF,quantile=quan))
}

# Define true progession functions of T

# 1. Sigmoid

f_sigmoid <- function(t, a=1, b=0, d=1){
  return(a/(1+exp(-(t-b)/d)))
} 

# 2. General S-shape

f_sshape <- function(t, mode=0.5, range_L=0, range_R=1){
  return(ptri(t,range_L,range_R,mode))
}

# 3. Wiggled

f_wiggle <- function(t, mean1=0, sd1=1, mean2=1, sd2=1, p1=0.5, p2=0.5){
  return((pnorm(t,mean1,sd1)*p1+pnorm(t,mean2,sd2)*p2)/(p1+p2))
}

missing_pattern <- function(N,P,total_missing_prob=1/3,random_missing_prob=1/2){
  mat <- matrix(rep(T,N*P),nrow=N,ncol=P)
  missing_rows <- as.logical(rbinom(N,1,total_missing_prob))
  missing_elem <- matrix(as.logical(rbinom(N*P,1,total_missing_prob)),nrow=N,ncol=P)
  missing_elem[missing_rows,] <- TRUE
  return(missing_elem)
}

Linf_norm <- function(f1,f2,L=-Inf,R=+Inf){
  return(optim(0,function(x) abs(f1(x)-f2(x)),lower = L,upper = R,control=list(fnscale=-1)))$value
}
# General package loader
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
rtMVN <- function(mean,Sigma,posit=1:length(mean),
                  SShape=FALSE){
  if(is.null(posit)){
    return(as.vector(mvtnorm::rmvnorm(1,mean,Sigma)))
  }
  if(!(SShape)){
    while(TRUE){
      x <- as.vector(mvtnorm::rmvnorm(1,mean,Sigma))
      if(all(x[posit]>=0)) return(x)
    }
  }
  else
  {
    while(TRUE){
      x <- as.vector(mvtnorm::rmvnorm(1,mean,Sigma))
      x_res <- x[posit]
      if(!all(x[posit]>=0)) next
      x_comp <- x_res[-1] >= x_res[-length(x_res)]
      if(length(x_comp)==1) return(x)
      if(all(x_comp[-1]<=x_comp[-length(x_comp)])) return(x)
    }
  }
}
# Updating Parameters --------------------------------------------
# 0816Change: Included RE variables, individual-biomarker specific random effect
# 0821Change: Included separate spline covariates the same number as biomarkers

# Update BETA and GAMMA together----
# RE is the replicated random effects
# RE should always be the same dimension as Y
update_coef <- function(covars.list,nX,Y,RE,sy,sw,id,prior.mean,prior.precision,samples=1,burnin=5){
  res <- array(0,c(samples,ncol(covars.list[[1]]),ncol(Y)))
  unique_id <- unique(id)
  for(k in 1:ncol(Y)){
    non_mis <- !is.na(Y[,k])
    CC <- covars.list[[k]][non_mis,]
    comp_id <- id[non_mis]
    block_size <- integer(length(unique_id))
    for(j in 1:length(block_size))
      block_size[j] <- sum(comp_id==j)
    lik_prec <- exchangable_cov(sy,sw,block_size)$Prec
    precision <- prior.precision + t(CC)%*%lik_prec%*%CC
    variance <- solve(precision)
    
    mu <- as.vector(precision%*%prior.mean)
    mu <- mu + t(CC)%*%lik_prec%*%Y[,k][non_mis]
    
    mu <- as.vector(variance%*%mu)
    # res[,k] <- rtmvnorm(1,mu,variance,
    #                     lb=c(rep(-Inf,nX),rep(0,ncol(covars.list[[k]])-nX))
    # )
    # res[,k] <- rtMVN(mu,variance,(nX+1):length(mu),FALSE)
    res[,,k] <- hdtg::harmonicHMC(samples,burnin,mu,chol(variance),
                                 diag(length(mu))[-(1:nX),],rep(0,length(mu)-nX),
                                 rep(0.1,length(mu)),precFlg=FALSE)$samples[(burnin+1):(burnin+samples),]
  }
  return(res)
}
# Update SIGMA_Y----
update_sigmay <- function(covars.list,Y,RE,coefs,prior.shape,prior.scale){
  residual <- array(0, dim(Y))
  for(k in 1:ncol(Y)){
    residual[,k] <- Y[,k]-RE[,k]-drop(covars.list[[k]]%*%coefs[,k])
  }
  return(1/rgamma(1,
                  shape=prior.shape + sum(!is.na(Y))/2,
                  rate=prior.scale + sum(residual^2,na.rm=TRUE)/2
  ))
}
# Update PENALTIES----
# Old version: NORMAL jump
update_pens0 <- function(gamma,mu,pen,lpd,ls,weight){
  p <- dim(gamma)[1]
  V <- penalty_Matrix(p, pen[1], pen[2],weight)$V
  ll <- sum(dtmvnorm(t(gamma),mu,V,lb=rep(0,p),log=TRUE)) +
    lpd(pen)
  pen_pro <- pen + rnorm(2, 0, exp(ls))
  V2 <- penalty_Matrix(p, pen_pro[1], pen_pro[2],weight)$V
  if(pen_pro[1]>0 & pen_pro[2]>0){
    ll2 <- sum(dtmvnorm(t(gamma),mu,V2,lb=rep(0,p),log=TRUE)) +
      lpd(pen_pro)
  }
  else
    ll2 <- -Inf
  diff <- ll - ll2
  temp <- rexp(1)
  if(temp>diff){
    acc_status <- 1
    new <- pen_pro
  }
  else
  {
    acc_status <- 0
    new <- pen
  }
  return(list(acc_status=acc_status,
              new=new))
}
# New Version: Log-normal jump
update_pens <- function(
    gamma, # Matrix containing K gamma coefficients
    mu, # Prior Mean of Gamma
    lambda, # Current value of penalty parameters, 1 for smooth, 2 for flat 
    lpd, # Log posterior density function of lambda
    ls, # Log of standard deviation of proposal distribution
    weight # Window Function; NULL means default quadratic window
){
  p <- dim(gamma)[1]
  V <- penalty_Matrix(p, lambda[1], lambda[2], weight)$V # Current variance
  # Log-density at current value
  ll <- sum(dtmvnorm(t(gamma),mu,V,lb=rep(0,p),log=TRUE)) + lpd(lambda)
  # Propose new state
  # Since lambda is positive, we use log-normal jumps
  new <- exp(log(lambda) + rnorm(2, sd=exp(ls)))
  
  V2 <- penalty_Matrix(p, new[1], new[2], weight)$V # New variance
  # Log-density at new value
  ll2 <- sum(dtmvnorm(t(gamma),mu,V2,lb=rep(0,p),log=TRUE)) + lpd(new)
  
  diff <- ll - sum(log(new)) - ll2 + sum(log(lambda))
  # diff <- ll - ll2
  temp <- rexp(1)
  if(temp>diff){
    acc_status <- 1
    state <- new
  }
  else
  {
    acc_status <- 0
    state <- lambda
  }
  return(list(acc_status=acc_status,
              new=state))
}

# Update Random Effect VARIANCE----
# W is the non-replicated version of RE
# W should contain N*K elements, where N is the number of individuals
# and K is the number of biomarkers
update_sigmaw <- function(W,prior.shape,prior.scale){
  post.shape <- prior.shape + length(W)/2
  post.scale <- prior.scale + sum(W^2)/2
  return(1/rgamma(1,shape=post.shape,
                  rate=post.scale))
}
# Update Individual Random Intercept----
update_W <- function(covars.list,Y,coefs,long_ss,ID,
                     sigmay, sigmaw)
{
  residual <- array(0, dim(Y))
  for(k in 1:ncol(Y)){
    residual[,k] <- Y[,k] - drop(covars.list[[k]]%*%coefs[,k])
  }
  sds <- (1/sigmaw + long_ss/sigmay)^(-1/2)
  means <- matrix(0,nrow=nrow(long_ss),ncol=ncol(long_ss))
  for(k in 1:ncol(long_ss)){
    means[,k] <- tapply(residual[,k],ID,sum,na.rm=TRUE)
  }
  means <- means / (sigmay/sigmaw + long_ss)
  norms <- array(rnorm(prod(dim(long_ss))),
                 dim=dim(long_ss))
  return(norms*sds+means)
}