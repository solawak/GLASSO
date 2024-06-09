# Katarzyna Solawa, Małgorzata Kaczkowska, Martyna Czajkowska
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

# Implementacja LASSO 
lasso <- function(beta, rho, V,u) {
  for (j in 1:ncol(V)) {
    x <- u[j] - V[-j,j]%*%beta[-j]
    beta[j] <- sign(x)*max(abs(x)-rho,0)/V[j,j]
  }
  return(beta)
}

# Implementacja GLASSO, funkcja zwraca macierz odwrotną do macierzy kowariancji i macierz kowariancji
GLASSO_ <- function(S, rho=1.2, t=0.001, max_iter=100000,penalize_diag =TRUE) {
  p <- ncol(S)
  W = S + rho*diag(1,p,p)*penalize_diag
  Theta = matrix(NA, p,p)
  #inicjujemy beta = 0
  betas <- matrix(0,p-1,p)
  for (iter in 1:max_iter) {
    # zapisujemy poprzednią wartość W jako w_old
    W_old <- W
    for (i in 1:p) {
      #kolejne kroki lasso liczymy na W_old
      betas[,i] <- lasso(betas[,i], rho, W_old[-i,-i], S[-i,i])
      W[-i,i] <- W_old[-i,-i]%*%betas[,i] #W12 = W11*Beta
      W[i,-i] <- W_old[-i,-i]%*%betas[,i] #W21 = W11*Beta
    }
    if (mean((abs(W-W_old))) < t*(sum(abs(S)) - sum(abs(diag(S))))/(p^2 - p)) { # kryterium stopu
      break
    }
  }
  
  for (i in 1:p) {  
    Theta[i,i] <- 1/(W[i,i] - W[-i,i]%*%betas[,i]) # 1/(W22-W12*Beta)
    Theta[-i,i] <- -Theta[i,i]*betas[,i] #Theta12 = -Theta22*Beta
    Theta[i,-i] <- -Theta[i,i]*betas[,i] #Theta21 = -Theta22*Beta
  }
  return(list(Theta, W))
}

# Porównanie dla penalize_diag = FALSE
X <- GLASSO_(S, rho=1.2, penalize_diag = FALSE) 
mdl <- glasso(S,rho=1.2,approx=FALSE, penalize.diagonal = FALSE) # wbudowana funkcja
# Estymowana macierz odwrotna do macierzy kowariancji
X[1]
mdl$wi
# Estymowana macierz kowariancji
X[2]
mdl$w
# Wniosek: otrzymaliśmy wyniki zbliżone do funkcji glasso

# Porównanie dla penalize_diag = TRUE
X <- GLASSO_(S, rho=1.2, penalize_diag = TRUE) 
mdl <- glasso(S,rho=1.2,approx=FALSE, penalize.diagonal = TRUE) # wbudowana funkcja
# Estymowana macierz odwrotna do macierzy kowariancji 
X[1]
mdl$wi
# Estymowana macierz kowariancji
X[2]
mdl$w
# Wniosek: otrzymaliśmy wyniki zbliżone do funkcji glasso
