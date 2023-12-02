# testing
source('toy.R')
M <- matrix(0,4,4)
M[,] <- 0.5
diag(M) <- 1
# Bayesian Coverage: True para generated from prior
print(toy_simu(4,NULL,prior.var=diag(4)*100,lik.var=M,verbose=FALSE,restrict='S')$coverage)
print(toy_simu(4,NULL,prior.var=diag(4)*100,lik.var=M,verbose=FALSE,restrict='Flex')$coverage)
print(toy_simu(4,NULL,prior.var=diag(4)*100,lik.var=M,verbose=FALSE,restrict='None')$coverage)
# Frequentist Coverage: True para fixed
print(toy_simu(4,c(6, 8, 10, 4),prior.var=diag(4)*100,lik.var=M,verbose=FALSE,restrict='S')$coverage)
print(toy_simu(4,c(6, 8, 10, 4),prior.var=diag(4)*100,lik.var=M,verbose=FALSE,restrict='Flex')$coverage)
print(toy_simu(4,c(6, 8, 10, 4),prior.var=diag(4)*100,lik.var=M,verbose=FALSE,restrict='None')$coverage)