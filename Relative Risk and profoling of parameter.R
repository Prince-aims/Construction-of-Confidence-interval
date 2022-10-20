#Relative Risk (Risk Ratio) of Pirenzepine drug compared to Trithiozine drug
y1 <- 23; y2 <-18; n1 <-30; n2 <-31
(RRest <- (y1*n2)/(y2*n1));

# Modified-Wald method for Relative Risk (Risk Ratio) of Pirenzepine drug
#compared to Trithiozine drug
RR.MWCI <- function(y1, n1, y2, n2, alf) {
  theta1hat <- y1/n1
  theta2hat <- y2/n2
  z <- qnorm(1 - 0.5*alf)
  exp(log(theta1hat/theta2hat) +
        c(-1, 1)*z*sqrt((1 - theta1hat)/(y1) + (1-theta2hat)/(y2)))
}
RR.MWCI(23, 30, 18, 31, 0.05) 
#Score Confidence interval for Relative Risk (Risk Ratio) of Pirenzepine drug compared to
#Trithiozine drug
#install.packages("PropCIs")
library(PropCIs)
riskscoreci(y1, n1, y2, n2, conf.level = 0.95)

# Function to Plot Modified Wald CI curve.
PWaldcurveRR <- function(phi, y1, n1, y2, n2) {
  phiHat <- (y1*n2)/(y2*n1)
  bot <- (n1 - y1)/(n1*y1) + (n2 - y2)/(n2*y2)
  top <- (log(phi) - log(phiHat))^2
  -0.5*top/bot
}
#Function to plot Profile likelihood curve 
PLcurveRR <- function(phi, y1, n1, y2, n2) {
  y <- y1 + y2
  n <- n1 + n2
  theta1hat <- y1/n1
  theta2hat <- y2/n2
  LLhat <- y1*log(theta1hat) + y2*log(theta2hat) + (n1 - y1)*log(1 - theta1hat) +
    (n2 - y2)*log(1 - theta2hat)
  del <- (y1 + n2 + phi*(y2 + n1))^2 - 4*phi*n*y
  top <- y1 + n2 + phi*(y2 + n1)-sqrt(del)
  bot <- 2*n*phi
  pi2tilde <- top/bot
  pi1tilde <- phi*pi2tilde
  LLtilde <- y1*log(pi1tilde) + y2*log(pi2tilde) + (n1 - y1)*log(1 - pi1tilde) +
    (n2 - y2)*log(1 - pi2tilde)
  LLtilde - LLhat
}
#Intializing scalling values 
phia <- seq(0.04, 11, length = 1000)
scLL <- PLcurveRR(phia, 23, 30, 18, 31)
phiW <- seq(0.2,11, length = 1000)
sWLL <- PWaldcurveRR(phiW, 23, 30, 18, 31)
xrange <- range(phia,phiW, 6)
yrange <- range(scLL,sWLL, -4)
xx <- rep(1.3203, 2)
yy <- seq(-64.74816, 0, length = 2)
#Plotting PLcurve
plot(phia, scLL, type = "l", lwd = 3,
     xlab = expression(phi * " (Relative Risk)"),
     ylab = "Scaled Constrained Log Likelihood",xlim = xrange, ylim = yrange)
par(new = T)
plot(phiW, sWLL, type = "l", lwd = 3, lty = 3, col = "blue",
     xlab = "", ylab = "", xlim = xrange, ylim = yrange)
par(new = T)
plot(xx, yy, type = "l", lwd = 1, lty = 3, xlim = xrange, ylim = yrange,
     xlab = "", ylab = "")
#Setting the points of PLCI and Wald CI
points(1.3203, 0, pch = 19, cex = 2,col=c("red"))
points(0.9226, -0.5*qchisq(0.95,1), pch = 19, cex = 1.5)
points(1.8895, -0.5*qchisq(0.95,1), pch = 19, cex = 1.5)
points(0.929656, -0.5*qchisq(0.95,1), pch = 17, cex = 1.5, col = c("blue"))
points(1.960705, -0.5*qchisq(0.95,1), pch = 17, cex = 1.5, col = c("blue"))
legend("bottom",
       legend = c("95% PLCI (0.93,1.96)",
                  "95% MWCI (0.92,1.89)"), pch = c(17, 19),
       col = c("blue", "black"))
legend("topright", c("PL Curve", "Wald Approx Curve"), lty = c(1, 3),
       lwd = 3, col = c("black", "blue"))
title("Scaled Constrained LL: Relative Risk Parameter")