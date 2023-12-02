load("flex_CIs.rda")
flex_CIs <- CI_repeat
covered_flex <- apply(flex_CIs,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate_flex <- apply(covered_flex, 2, mean)
load("S_CIs.rda")
s_CIs <- CI_repeat
covered_S <- apply(s_CIs,c(1,2),function(x) (x[4]-x[2])*(x[4]-x[3])<=0)
cover_rate_S <- apply(covered_S, 2, mean)
ages <- seq(0,120,by=0.1)
pdf("coverage_comparison_trueSpline_hdtg.pdf",width=14)
plot(ages,cover_rate_flex,xlab='age',ylab='coverage',ylim=c(0,1),type='l',
     col='red')
lines(ages,cover_rate_S,col='blue')
abline(h=0.95,lty=2)
legend(100,0.6,legend=c("Flexible","S-shape"),col=c("red","blue"),
       lty=c(1,1))
dev.off()