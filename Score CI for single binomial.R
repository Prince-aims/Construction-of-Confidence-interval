# Score
binomSCI <- function(y, n, alf) {
  pihat <- y/n
  zz <- qnorm(1 - 0.5*alf)
  aa <- 1 + zz*zz/n
  bb <- 2*pihat + zz*zz/n
  cc <- pihat*pihat
  fn <- function(pi) aa*pi*pi - bb*pi + cc
  solL <- uniroot(fn, c(0, pihat))$root
  solR <- uniroot(fn, c(pihat, 1))$root
  sol <- cbind(solL, solR)
  sol
}
binomSCI(1, 15, 0.05)

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
abline(v=binomSCI(1, 15, 0.05),col="green")
title("logLikelihood Curve with SCI")
legend("bottomright",legend = c("ScoreCI"),
       pch = c(19, 15),
       col = c("green"),lwd = 3)

