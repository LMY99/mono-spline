lwmean <- function(a,b,w=1/2){
  require(matrixStats)
  require(Rmpfr)
  a0=a+log(w);b0=b+log(1-w)
  mat=rbind(a0,b0)
  return(colLogSumExps(mat))
}

logpmvnorm <- function(lb,ub,mu,Sigma,Nmax=1e3){
  require(matrixStats)
  require(VGAM)
  require(Rmpfr)
  precision=100
  a=lb-mu; b=ub-mu
  #a=mpfr(a,precision);b=mpfr(b,precision)
  # Reordering
  # double_inf=which((a==-Inf)&(b==+Inf))
  # no_inf=which(!is.infinite(b-a))
  # no_inf=no_inf[order((b-a)[no_inf],decreasing=FALSE)]
  # a_inf=which((a==-Inf)&(b!=+Inf))
  # b_inf=which((a!=-Inf)&(b==+Inf))
  # score=rep(NA,length(a))
  # score[a_inf]=b[a_inf]; score[b_inf]=-a[b_inf]
  # one_inf=c(a_inf,b_inf)
  # one_inf=one_inf[order(score[one_inf],decreasing=FALSE)]
  # perm=c(no_inf,one_inf,double_inf)
  # a=a[perm]; b=b[perm]; Sigma=Sigma[perm,perm];
  # Monte-Carlo
  C=t(chol(Sigma))
  m=ncol(C)
  #d=rep(0,m);e=rep(0,m);f=rep(0,m);y=rep(0,m-1)
  d=array(0,dim=c(Nmax,m))
  e=array(0,dim=c(Nmax,m))
  f=array(0,dim=c(Nmax,m))
  y=array(0,dim=c(Nmax,m-1))
  temp=array(0,dim=c(Nmax,m-1))
  #d=mpfr(d,precision);e=mpfr(e,precision);f=mpfr(f,precision)
  #y=mpfr(y,precision);temp=mpfr(temp,precision)
  d[,1]=ifelse(a[1]==-Inf,-Inf,pnorm(a[1]/C[1,1],log=T));
  e[,1]=ifelse(b[1]==+Inf,0,pnorm(b[1]/C[1,1],log=T));
  f[,1]=e[1,1]+log1mexp(e[1,1]-d[1,1])
  
  # for(N in 1:Nmax){
  #   w=runif(m-1)
  #   for(i in 2:m){
  #     #print(c(e[i-1],d[i-1]))
  #     y[i-1]=qnorm(exp(d[i-1])+w[i-1]*(exp(e[i-1])-exp(d[i-1])),log=F)
  #     #if(N==1) print(c(w[i-1],y[i-1]))
  #     d[i]=ifelse(a[i]==-Inf,-Inf,
  #       pnorm((a[i]-sum(C[i,1:(i-1)]*y[1:(i-1)]))/C[i,i],log=T)          
  #     )
  #     e[i]=ifelse(b[i]==+Inf,0,
  #       pnorm((b[i]-sum(C[i,1:(i-1)]*y[1:(i-1)]))/C[i,i],log=T)        
  #     )
  #     #if(N==1) print(c(d[i],e[i]))
  #     f[i]=e[i]+log1mexp(e[i]-d[i])+f[i-1]
  #     #if(N==1) print(f[i])
  #   }
  #   #print(y)
  #   #if(N==1) print("***")
  #   #delta=f[m]+log1p(-exp(Intsum-f[m]))-log(N)
  #   #Intsum=logSumExp(c(Intsum,delta))
  #   ff[N]=f[m]
  # }
  #cat(c(sum(is.na(d[,1])),sum(is.na(e[,1])),sum(is.na(f[,1]))))
  w=array(runif((m-1)*Nmax),dim=c(Nmax,m-1))
  for(i in 2:m){
    temp[,i-1]=lwmean(d[,i-1],e[,i-1],w[,i-1])
    temp[temp[,i-1]>=0,i-1]=0
    y[,i-1]=qnorm(temp[,i-1],log=T)
    y[y[,i-1]==+Inf,i-1]=1e3
    y[y[,i-1]==-Inf,i-1]=-1e3
    if(a[i]==-Inf)
      d[,i]=-Inf
    else
      d[,i]=pnorm((a[i]-colSums(C[i,1:(i-1)]*t(y[,1:(i-1)])))/C[i,i],log=T)
    if(b[i]==+Inf)
      e[,i]=0
    else
      e[,i]=pnorm((b[i]-colSums(C[i,1:(i-1)]*t(y[,1:(i-1)])))/C[i,i],log=T)
    diff=log1mexp(e[,i]-d[,i])
    diff=ifelse(diff!=-Inf,diff,log1mexp(.Machine$double.xmin))
    f[,i]=e[,i]+diff+f[,i-1]
    #cat(c(sum(is.na(d[,i])),sum(is.na(e[,i])),sum(is.na(f[,i]))))
  }
  res=logSumExp(f[,m],na.rm=TRUE)-log(sum(!is.na(f[,m])))
  if(is.na(res)||res==-Inf){ 
    res=(-Inf)
    View(f); View(d); View(e); View(y); View(temp)
    stop("Error")
    attr(res,"rel_error")=0
    return(res)
  }
  else{
    v=cbind(f[,m],res)
    idx=v[,1]<v[,2]
    idx[is.na(idx)]=FALSE
    v[idx,]=v[idx,c(2,1)]
    ad=v[,1]+log1mexp(v[,1]-v[,2])
    ad=ad*2
    var0=logSumExp(ad,na.rm=TRUE)-log(sum(!is.na(ad)))
    var0=var0-log(sum(!is.na(ad)))
    sd0=var0/2
    attr(res,"rel_error")=exp(sd0-res)
  }
  #cat(sum(is.na(ad))," ",sum(is.na(f[,m])), "\n")
  return(res)
}

hdtg_S <- function(n,mu,sigma){
  p <- length(mu)
  FF <- matrix(0,nrow=p+1,ncol=p)
  for(i in 1:p){
    FF[i,i] <- 1
    if(i>1) FF[i,i-1] <- -1
  }
  FF[p+1,p] <- 1
  logprob <- rep(0,p)
  for(i in 1:p){
    Fmat <- FF
    if(i<p) Fmat[(i+1):p,] <- -Fmat[(i+1):p,]
    g <- rep(0,p+1)
    new.S <- Fmat%*%sigma%*%t(Fmat); new.S[p+1,p+1] <- new.S[p+1,p+1]*(1+1e-3)
    new.mu <- Fmat%*%mu
    logprob[i] <- logpmvnorm(rep(0,p+1),rep(+Inf,p+1),new.mu,new.S)
  }
  logprob <- logprob - max(logprob)
  inflex.points <- sample(1:p,n,replace=TRUE,prob=exp(logprob))
  result <- matrix(0,nrow=n,ncol=p)
  for(i in 1:p){
    pnums <- sum(inflex.points==i)
    if(pnums==0) next
    Fmat <- FF
    if(i<p) Fmat[(i+1):p,] <- -Fmat[(i+1):p,]
    Fmat <- FF
    if(i<p) Fmat[(i+1):p,] <- -Fmat[(i+1):p,]
    g <- rep(0,p+1)
    init <- rep(0,p)
    init[i] <- 0.1
    init[-i] <- 0.1*runif(p-1)
    if(i>1) init[1:(i-1)] <- sort(init[1:(i-1)],decreasing=FALSE)
    if(i<p) init[(i+1):p] <- sort(init[(i+1):p],decreasing=TRUE)
    result[inflex.points==i,] <- hdtg::harmonicHMC(n=pnums,burnin=100,
                                                   mean=mu,choleskyFactor=chol(sigma),
                                                   F=Fmat,g=g,init=init,precFlg=FALSE)$samples[101:(100+pnums),]
  }
  if(n==1) return(as.vector(result))
  else return(result)
}

toy_simu <- function(dim, 
                     true.mean=NULL,
                     prior.mean=rep(0,dim),
                     prior.var=diag(dim),
                     lik.var=diag(dim),
                     restrict="None",
                     iter=100,
                     mcmc.size=1000,
                     seed=Sys.time(),
                     verbose=TRUE
                     )
{
  require(hdtg)
  require(mvtnorm)
  set.seed(seed)
  post.var <- solve(solve(prior.var)+solve(lik.var))
  result <- list()
  result$samples <- array(0,dim=c(iter,mcmc.size,dim))
  result$true.mean <- matrix(0,nrow=iter,ncol=dim)
  for(i in 1:iter){
    if(verbose) cat(sprintf("%d ",i))
    if(is.null(true.mean)){
      if(restrict=='None'){
        mu <- rmvnorm(n=1,mean=prior.mean,sigma=prior.var)
      }
      if(restrict=='Flex'){
        mu <- harmonicHMC(1,100,prior.mean,chol(prior.var),F=diag(dim),g=rep(0,dim),
                          init=rep(1,dim),precFlg=FALSE)$samples[101,]
      }
      if(restrict=='S'){
        mu <- hdtg_S(1,prior.mean,prior.var)
      }
    }
    else{
      mu <- true.mean
    }
    result$true.mean[i,] <- mu
    y <- rmvnorm(n=1,mu,lik.var)[1,]
    post.mean <- post.var %*% (solve(prior.var,prior.mean)+solve(lik.var,y))
    if(restrict=='None'){
      result$samples[i,,] <- rmvnorm(mcmc.size,post.mean,post.var)
    }
    if(restrict=='Flex'){
      result$samples[i,,] <- harmonicHMC(mcmc.size,100,post.mean,chol(post.var),F=diag(dim),g=rep(0,dim),
                                  init=rep(1,dim),precFlg=FALSE)$samples[101:(100+mcmc.size),]
    }
    if(restrict=='S'){
      result$samples[i,,] <- hdtg_S(mcmc.size,post.mean,post.var)
    }
  }
  result$CI <- array(0,c(iter,2,dim))
  for(i in 1:iter){
    result$CI[i,1,] <- apply(result$samples[i,,],2,quantile,probs=0.025)
    result$CI[i,2,] <- apply(result$samples[i,,],2,quantile,probs=0.975)
  }
  result$coverage <- apply((result$CI[,1,]-result$true.mean)*(result$CI[,2,]-result$true.mean)<=0,2,mean)
  return(result)
}