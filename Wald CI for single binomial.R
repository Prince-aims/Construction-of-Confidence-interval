# Wald
binomWCI <- function(y, n, alf) {
  pihat <- y/n
  SE <- sqrt(pihat*(1-pihat)/n)
  pihat + c(-1,1)*qnorm(1 - 0.5*alf)*SE
}
binomWCI(1,15, 0.05)

#loglikelihood curve with confidence intervals 
y=1
n=15
thetahat<-y/n
chisq <- qchisq(1-0.05, 1)
LLcurve<-function(theta){return(
  y*log(theta) + (n - y)*log(1 -theta)
)}
p<-seq(-0.1,1,length=100000)
plot(p,LLcurve(p),xlab=expression(theta),ylab = "Log-likelihood",ylim=LLcurve(thetahat)+c(-5,0),type="l")
abline(v=binomWCI(1, 15, 0.05),col="blue")
title("logLikelihood Curve with WCI")
legend("topright",legend = c("waldCI"),
       pch = c(19, 15),
       col = c("blue"),lwd = 3)


