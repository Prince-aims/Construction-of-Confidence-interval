# Odds Scale
# Likelihood
binomOddsLCI <- function(y, n, alf) {
  thetahat <- y/n
  kaphat <- thetahat/(1 - thetahat)
  chisq <- qchisq(1-alf, 1)
  num <- -y*log(y) - (n - y)*log(n - y) + n*log(n) + 0.5*chisq
  fn <- function(kap) y*log(kap) - n*log(1 + kap) + num
  solL <- uniroot(fn, c(0, kaphat))$root
  solR <- uniroot(fn, c(kaphat, 20*kaphat))$root
  intermed <- cbind(solL, solR)
  1/(1 + 1/intermed)
}
binomOddsLCI(1, 15, 0.05) 

#WCI
binomOddsWCI <- function(y, n, alf) {
  thetahat <- y/n
  kaphat <- thetahat/(1 - thetahat)
  SE <- (1 + kaphat)*sqrt(kaphat/n)
  intermed <- kaphat + c(-1, 1)*qnorm(1 - 0.5*alf)*SE
  result <- 1/(1 + 1/intermed)
  result
}
binomOddsWCI(1, 15, 0.05) 

# Log-Odds Scale.

# Likelihood
binomLogOddsLCI <- function(y, n, alf) {
  thetahat <- y/n
  lamhat <- log(thetahat/(1 - thetahat))
  chisq <- qchisq(1 - alf, 1)
  num <- -y*log(y) - (n - y)*log(n - y) + n*log(n) + 0.5*chisq
  fn <- function(lam) y*lam - n*log(1 + exp(lam)) + num
  solL <- uniroot(fn, c(0, lamhat))$root
  solR <- uniroot(fn, c(lamhat, 20*lamhat))$root
  intermed <- cbind(solL, solR)
  sol <- 1/(1 + exp(-intermed))
  sol
}
binomLogOddsLCI(1, 15, 0.05) 

#WCI
binomLogOddsWCI <- function(y, n, alf) {
  thetahat <- y/n
  lamhat <- log(thetahat/(1 - thetahat))
  SE <- (1 + exp(lamhat))/(sqrt(n*exp(lamhat)))
  intermed <- lamhat + c(-1, 1)*qnorm(1 - 0.5*alf)*SE
  sol <- 1/(1 + exp(-intermed))
  sol
}
binomLogOddsWCI(1, 15, 0.05) # 0.009305387 0.351990342

