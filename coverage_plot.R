load("flex_CIs.rda")
flex_CIs <- CI_repeat
covered_flex <- apply(flex_CIs,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate_flex <- apply(covered_flex, 2, mean)
sink("Covariate and RE coverage.txt")
cat("Target: 0.95\n")
cat("Flexible model(Covariate):\n")
cat(apply(apply(CI_covariate_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0),2,mean))
cat("\nFlexible model(RE):\n")
cat(mean(apply(RE_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nFlexible model(Offset):\n")
cat(mean(apply(offset_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))

load("S_CIs.rda")
s_CIs <- CI_repeat
covered_S <- apply(s_CIs,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate_S <- apply(covered_S, 2, mean)
cat("\nS-Shaped model(Covariate):\n")
cat(apply(apply(CI_covariate_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0),2,mean))
cat("\nS-Shaped model(RE):\n")
cat(mean(apply(RE_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nS-Shaped model(Offset):\n")
cat(mean(apply(offset_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))

sink()

ages <- seq(0,120,by=0.1)

length_flex <- apply(flex_CIs,c(1,2),function(x) abs(x[3]-x[2]))
length_S <- apply(s_CIs,c(1,2),function(x) abs(x[3]-x[2]))

pdf("coverage_comparison_trueSpline_hdtg.pdf",width=14)
plot(ages,cover_rate_flex,xlab='age',ylab='coverage',ylim=c(0,1),type='l',
     col='red',main='Pointwise Credible Band Coverage')
lines(ages,cover_rate_S,col='blue')
abline(h=0.95,lty=2)
abline(v=mean(true_turning),lty=2,col='lightblue')
legend('bottomright',legend=c("Flexible","S-shape","Target","True Inflection Point(Avg)"),col=c("red","blue","black","lightblue"),
       lty=c(1,1,2,2))
text(x=mean(true_turning),y=0,labels=sprintf("Inflection point coverage: %.3f",
                                             mean((turning[,1]-true_turning)*(turning[,2]-true_turning)<=0)))

plot(ages,colMeans(length_flex),type='l',xlab='Age',ylab='Width',main='Pointwise Credible Band Width',
     col='red')
lines(ages,colMeans(length_S),col='blue')
legend('bottomright',legend=c("Flexible","S-shape"),col=c("red","blue"),
       lty=c(1,1))
plot(ages,colMeans(length_flex)/colMeans(length_S),type='l',xlab='Age',ylab='CI Width Ratio',
     main='Relative Efficiency of S-shape vs Flex')
abline(h=1,lty=2)

plot(ages,colMeans(flex_CIs[,,6]),type='l',xlab='Age',ylab='MSE',main='Average posterior MSE',
     col='red')
lines(ages,colMeans(s_CIs[,,6]),col='blue')
legend('bottomright',legend=c("Flexible","S-shape"),col=c("red","blue"),
       lty=c(1,1))
dev.off()

cat("")