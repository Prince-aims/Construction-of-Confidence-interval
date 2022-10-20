# Likelihood
binomLCI<-function(y, n, alf) {
  pihat <- y/n
  chisq <- qchisq(1 - alf, 1)
  num <- -y*log(y)-(n - y)*log(n - y) + n*log(n) + 0.5*chisq
  fn <- function(theta) y*log(theta) + (n - y)*log(1 -theta) + num
  solL <- uniroot(fn, c(0,pihat))$root
  solR <- uniroot(fn, c(pihat,1))$root
  sol <- cbind(solL, solR)
  sol
}
binomLCI(1, 15, 0.05)
#loglikelihood curve with confidence intervals 
y=1
n=15
thetahat<-y/n
chisq <- qchisq(1 -0.05, 1)
LLcurve<-function(theta){return(
  y*log(theta) + (n - y)*log(1 -theta)
)}
p<-seq(0,1,length=100000)
plot(p,LLcurve(p),xlab=expression(theta),ylab = "Log-likelihood",ylim=LLcurve(thetahat)+c(-6,0),type="l")
#curve(LL,from=0.0001,to=0.7,n=10000,xlab=expression(theta),ylab = "LL",ylim=LL(thetahat)+c(-6,0))
abline(h= LLcurve(thetahat)-0.5*chisq,lty="dashed",col='red')
abline(v=binomLCI(1, 15, 0.05),col="yellow",pch=c(10))
title("Likelihood Curve with LBCI")
legend("bottomright",legend = c("Cutline","LCI"),
       pch = c(19, 15),
       col = c("red", "yellow"),lwd = 3)
LLcurve(thetahat)
LLcurve(thetahat)-0.5*chisq

