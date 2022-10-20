
# Likelihood
binomLCI <- function(y, n, alf) {
  pihat <- y/n
  chisq <- qchisq(1 - alf, 1)
  num <- -y*log(pihat) -(n - y)*log(1 -pihat)  + 0.5*chisq
  fn <- function(theta) y*log(theta) + (n - y)*log(1 -theta) + num
  solL <- uniroot(fn, c(0,pihat))$root
  solR <- uniroot(fn, c(pihat,1))$root
  sol <- cbind(solL, solR)
  sol
}
binomLCI(1, 15, 0.05)
-y*log(y)-(n - y)*log(n - y) + n*log(n)
# Wald
binomWCI <- function(y, n, alf) {
  pihat <- y/n
  SE <- sqrt(pihat*(1-pihat)/n)
  pihat + c(-1,1)*qnorm(1 - 0.5*alf)*SE
}
binomWCI(1,15, 0.05)
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
LLcurve<-function(theta){return(
  y*log(theta) + (n - y)*log(1 -theta)
)}
p<-seq(-0.1,1,length=100000)
plot(p,LLcurve(p),xlab=expression(theta),ylab = "Log-likelihood",ylim=LLcurve(thetahat)+c(-6,0),type="l")
abline(h= LLcurve(thetahat)-0.5*chisq,lty="dashed",col='red')
abline(v=binomLCI(1, 15, 0.05),col='yellow')
abline(v=binomSCI(1, 15, 0.05),col='green')
abline(v=binomWCI(1, 15, 0.05),col='blue')
legend("bottomright",
       legend = c("cutline","waldCI","Likelihood based CI","ScoreCI"),
       pch = c(19, 15, 17, 18),
       col = c("red", "blue", "yellow","green"),lwd = 3)
title("log-likelihood curve with Wald CI, 
      Likelihood based CI and ScoreCI")









