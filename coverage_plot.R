load("flex_CIs.rda")
flex_CIs <- CI_repeat
covered_flex <- apply(flex_CIs,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate_flex <- apply(covered_flex, 2, mean)
sink("coverage.txt")
cat("Target: 0.95\n")
cat("Flexible model(Covariate):")
cat(apply(apply(CI_covariate_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0),2,mean))
cat("\nFlexible model(RE):")
cat(mean(apply(RE_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nFlexible model(Offset):")
cat(mean(apply(offset_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nFlexible model(Residual variance):")
cat(mean((sigmay_repeat[,4]-sigmay_repeat[,2])*(sigmay_repeat[,4]-sigmay_repeat[,3])<=0))
cat("\nFlexible model(Random effect variance):")
cat(mean((sigmaw_repeat[,4]-sigmaw_repeat[,2])*(sigmaw_repeat[,4]-sigmaw_repeat[,3])<=0))

load("S_CIs.rda")
s_CIs <- CI_repeat
covered_S <- apply(s_CIs,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate_S <- apply(covered_S, 2, mean)
cat("\nS-Shaped model(Covariate):")
cat(apply(apply(CI_covariate_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0),2,mean))
cat("\nS-Shaped model(RE):")
cat(mean(apply(RE_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nS-Shaped model(Offset):")
cat(mean(apply(offset_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nS-Shaped model(Residual variance):")
cat(mean((sigmay_repeat[,4]-sigmay_repeat[,2])*(sigmay_repeat[,4]-sigmay_repeat[,3])<=0))
cat("\nS-Shaped model(Random effect variance):")
cat(mean((sigmaw_repeat[,4]-sigmaw_repeat[,2])*(sigmaw_repeat[,4]-sigmaw_repeat[,3])<=0))
cat("\nS-Shaped model(Inflection point):")
cat("\nTrue value: ");cat(mean(true_turning));
cat("\nCoverage: ");cat(mean((turning[,1]-true_turning)*(turning[,2]-true_turning)<=0))

load("para_CIs.rda")
para_CIs <- CI_repeat
covered_para <- apply(para_CIs,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate_para <- apply(covered_para, 2, mean)
cat("\nParametric model(Covariate):")
cat(apply(apply(CI_covariate_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0),2,mean))
cat("\nParametric model(RE):")
cat(mean(apply(RE_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nParametric model(Offset):")
cat(mean(apply(offset_repeat,c(1,2),function(x)(x[4]-x[2])*(x[4]-x[3])<=0)))
cat("\nParametric model(Residual variance):")
cat(mean((sigmay_repeat[,4]-sigmay_repeat[,2])*(sigmay_repeat[,4]-sigmay_repeat[,3])<=0))
cat("\nParametric model(Random effect variance):")
cat(mean((sigmaw_repeat[,4]-sigmaw_repeat[,2])*(sigmaw_repeat[,4]-sigmaw_repeat[,3])<=0))
cat("\nParametric model(Inflection point):")
cat("\nTrue value: ");cat(mean(true_turning));
cat("\nCoverage: ");cat(mean((turning[,2]-true_turning)*(turning[,3]-true_turning)<=0))

sink()

ages <- seq(0,120,by=0.1)

length_flex <- apply(flex_CIs,c(1,2),function(x) abs(x[3]-x[2]))
length_S <- apply(s_CIs,c(1,2),function(x) abs(x[3]-x[2]))
length_para <- apply(para_CIs,c(1,2),function(x) abs(x[3]-x[2]))

pdf("coverage_comparison.pdf",width=14)
plot(ages,cover_rate_flex,xlab='age',ylab='coverage',ylim=c(0,1),type='l',
     col='red',main='Pointwise Credible Band Coverage')
lines(ages,cover_rate_S,col='blue')
lines(ages,cover_rate_para,col='green')
abline(h=0.95,lty=2)
abline(v=mean(true_turning),lty=2,col='lightblue')
legend('bottomright',legend=c("Flexible","S-shape","Parametric","Target","True Inflection Point(Avg)"),
       col=c("red","blue","green","black","lightblue"),
       lty=c(1,1,2,2))

plot(ages,colMeans(length_flex),type='l',xlab='Age',ylab='Width',main='Pointwise Credible Band Width',
     col='red')
lines(ages,colMeans(length_S),col='blue')
lines(ages,colMeans(length_para),col='green')
legend('bottomright',legend=c("Flexible","S-shape","Parametric"),col=c("red","blue","green"),
       lty=c(1,1))
plot(ages,colMeans(length_flex)/colMeans(length_S),col='blue',type='l',xlab='Age',ylab='CI Width Ratio',
     main='Relative Efficiency vs Flex')
lines(ages,colMeans(length_flex)/colMeans(length_para),col='green',
     main='Relative Efficiency vs Flex')
abline(h=1,lty=2)

plot(ages,colMeans(flex_CIs[,,6]),type='l',xlab='Age',ylab='MSE',main='Average posterior MSE',
     col='red')
lines(ages,colMeans(s_CIs[,,6]),col='blue')
lines(ages,colMeans(para_CIs[,,6]),col='green')
legend('bottomright',legend=c("Flexible","S-shape","Parametric"),col=c("red","blue","green"),
       lty=c(1,1))

matplot(ages,t(flex_CIs[1:5,,1]),col='red',type='l',lty=1,lwd=0.5
        ,xlab='Age',ylab='Biomarker',main='Model Truth VS Flexible Posterior Mean Estimates',)
lines(ages,colMeans(flex_CIs[,,4]),lwd=2,lty=2)

matplot(ages,t(s_CIs[1:5,,1]),col='red',type='l',lty=1,lwd=0.5
        ,xlab='Age',ylab='Biomarker',main='Model Truth VS S-Shape Posterior Mean Estimates',)
lines(ages,colMeans(s_CIs[,,4]),lwd=2,lty=2)

matplot(ages,t(para_CIs[1:5,,1]),col='red',type='l',lty=1,lwd=0.5
        ,xlab='Age',ylab='Biomarker',main='Model Truth VS Parametric Posterior Mean Estimates',)
lines(ages,colMeans(para_CIs[,,4]),lwd=2,lty=2)

dev.off()

cat("")