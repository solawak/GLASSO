library(MASS)
library(glasso)
p <- 8
n <- 10
trueK <- matrix(0, nrow = p, ncol = p)
trueK[,1] <- -1 / sqrt(p)
trueK[1,] <- -1 / sqrt(p)
diag(trueK) <- 1
trueSigma <- solve(trueK)

set.seed(124)
z <- mvrnorm(n=n, mu = rep(0,p), Sigma = trueSigma)
S <- t(z) %*% z/n

lasso <- function(beta, rho, V,u) {
  for (j in 1:ncol(V)) {
    x <- u[j] - V[-j,j]%*%beta[-j]
    beta[j] <- sign(x)*max(abs(x)-rho,0)/V[j,j]
  }
  return(beta)
}

GLASSO_ <- function(S, rho=1.2, t=0.001, max_iter=100000) {
  p <- ncol(S)
  W = S + rho*diag(1,p,p)
  Theta = solve(S,diag(1,p,p))
  stop = FALSE
  for (iter in 1:max_iter) {
    for (i in 1:p) {
      beta <- -Theta[-i,i]/Theta[i,i] #-Theta12/Theta22
      beta <- lasso(beta, rho, W[-i,-i], S[-i,i])
      Theta[i,i] <- 1/(W[i,i] - W[-i,i]%*%beta) # 1/(W22-W12*Beta)
      Theta[-i,i] <- -Theta[i,i]*beta #Theta12 = -Theta22*Beta
      Theta[-i,i] <- -Theta[i,i]*beta #Theta21 = -Theta22*Beta
      W_lag <- W
      W[-i,i] <- W[-i,-i]%*%beta #W12 = W11*Beta
      W[i,-i] <- W[-i,-i]%*%beta #W21 = W11*Beta
      if (mean((abs(W-W_lag))) < t*(sum(abs(S)) - sum(abs(diag(S))))/(p^2 - p)) {
        stop = TRUE
        print("Break")
        print(iter)
        break
      }
    } 
    if (stop) {
      break
    }
  }
  return(list(Theta, W))
}

X<- GLASSO_(S)
mdl <- glasso(S,rho=1.2,approx=FALSE, penalize.diagonal = FALSE)
# Estimated inverse covariance matrix
X[1]
mdl$wi

