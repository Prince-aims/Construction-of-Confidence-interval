# HWE Likelihood nased CI
HWELCI <- function(n1,n2,n3,N,alf) {
  thetahat <- (2*n1+n2)/(2*N)
  chisq <- qchisq(1 - alf, 1)
  num <- -(2*n1+n2)*log(thetahat)-(n2+2*n3)*log(1-thetahat)+0.5*chisq
  fn <- function(theta) (2*n1+n2)*log(theta)+(n2+2*n3)*log(1-theta)+ num
  solL <- uniroot(fn, c(0,thetahat))$root
  solR <- uniroot(fn, c(thetahat,1))$root
  sol <- cbind(solL, solR)
  sol
}
HWELCI(125,34,10,169,0.05)

# HWE Wald confidence interval
HWEWCI<-function(n1,n2,N,alf){
  thetahat <-(2*n1+n2)/(2*N)
  SE <- sqrt(thetahat*(1-thetahat)/(2*N))
  thetahat+ c(-1,1)*qnorm(1 - 0.5*alf)*SE
}
HWEWCI(125,34,169,0.05)

###HWE Score confidence interval
HWESCI<-function(n1,n2,N,conf)
{thetahat <-(2*n1+n2)/(2*N)
z<-qnorm(1 - 0.5*conf)
score<-(thetahat+z^2/2/(2*N)+c(-1,1)*z*sqrt(thetahat*(1-thetahat)/(2*N)+z^2/4/(2*N)^2))/(1+z^2/(2*N))
return(score)}
HWESCI(125,34,169,0.05)	

#drawing log-likelihood curve
n1=125;n2=34;n3=10;N=169
thetahat <- (2*n1+n2)/(2*N)
thetahat
SE <- sqrt(thetahat*(1-thetahat)/(2*N))
SE
LLcurve<-function(theta){return(
  n1*log(theta^2)+n2*log(2*theta*(1-theta))+n3*log((1-theta)^2)
)}
p<-seq(0,1,length=100)
plot(p,LLcurve(p),xlab=expression(theta),ylab = "Log-likelihood",ylim=LLcurve(thetahat)+c(-100,0),type="l")
abline(h= LLcurve(thetahat)-0.5*chisq,lty="dashed",col='black')
abline(v=HWEWCI(125,34,169,0.05),col='blue')
abline(v=HWELCI(125,34,10,169,0.05),col='green',pch=c(19))
abline(v=HWESCI(125,34,169,0.05),col='red')
legend("bottomleft", legend = c("cutline","waldCI","Likelihood based CI","Score CI"),
pch = c(19, 16, 17),col = c("black","blue", "green","red"),lwd = 3)
title("Log-likelihood curve with HWEWCI, HWELCI and HWESCI")
#cutline Value
chisq <- qchisq(1 -0.05, 1)
LLcurve(thetahat)-0.5*chisq 

