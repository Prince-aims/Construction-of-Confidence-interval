
  
############################## BHH BOD DATASET###########

#initialization 
time <- c(1, 2, 3, 5, 7, 10)
yBOD <- c(109, 149, 149, 191, 213, 224)
(BOD2 <- data.frame(time, yBOD))
xrange <- range(0, 10.5)
yrange <- range(0, 250)
par(mfrow = c(1, 1))
plot(BOD2$time, BOD2$yBOD, pch = 19, cex = 1.5, col = "blue",
     xlim = xrange, ylim = yrange, xlab = "time (days)",
     ylab = "Biochemical Oxygen Demand (mg/L)")

#defining function
SE2 <- function(x, th1, th2) th1*(1 - exp(-th2*x))
#attaching data 
attach(BOD2)
#Calculating theta2 using linear regression model
lm(log(1 - yBOD/250) ~ time - 1)

#Fitting model and finding estimated theta1 and theta2
fit2 <- nls(yBOD ~ SE2(time, theta1, theta2), data = BOD2,
            start = list(theta1 = 250, theta2 = 0.2571), trace = F)
summary(fit2)

theta1e <- coef(fit2)[1]
theta2e <- coef(fit2)[2]

#ploting  
#par(mfrow = c(1, 1))
#plot(BOD2$time, BOD2$yBOD, pch = 19, cex = 1.5, col = "blue",
     #xlim = xrange, ylim = yrange, xlab = "time (days)",
     #ylab = "Biochemical Oxygen Demand (mg/L)")
xx <- seq(0, 10.5, length = 1000)
yy <- SE2(xx, theta1e, theta2e)
lines(xx, yy, lwd = 3, col = "black")
title("Box, Hunter & Hunter BOD Data and SE2 Fit")

#Calculating confidence interval for theta1 and theta2.
confint(fit2)
#Calculating confidence interval using Wald test(linear approximation) for theta1 and theta2
theta1e + c(-1, 1)*qt(0.975, 4)*12.3545
theta2e + c(-1, 1)*qt(0.975, 4)*0.1046
#ploting the confidence region using profiling.S
par(mfrow = c(2, 1))
plot(profile(fit2))

