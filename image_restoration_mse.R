library("png")
library(grid)
library(utils)
setwd("C:/Users/kasia/OneDrive/Dokumenty/GitHub/GLASSO/")

get_neighbors <- function(i, j, N, M, d) {
  # Define the potential neighbors
  if (d==1) {
    potential_neighbors <- list(c(i-1, j), c(i+1, j), c(i, j-1), c(i, j+1))
  }else {
    potential_neighbors <- list(c(i-1, j), c(i+1, j), c(i, j-1), c(i, j+1), c(i-1,j-1), c(i-1,j+1), c(i+1,j-1), c(i+1,j+1))
  }
  
  
  # Filter out neighbors that are outside the matrix
  valid_neighbors <- lapply(potential_neighbors, function(coord) {
    if (coord[1] >= 1 && coord[1] <= N && coord[2] >= 1 && coord[2] <= M) {
      return(coord)
    } else {
      return(NULL)
    }
  })
  # Remove NULL elements
  valid_neighbors <- valid_neighbors[!sapply(valid_neighbors, is.null)]
  return(valid_neighbors)
}


calculate_conditional_distribution2 <- function(i,j,x,y,col_pal,beta,alpha, lambda,sigma, method, d) {
  neighbors <- get_neighbors(i,j, nrow(x), ncol(x), d)
  conditional_distribution <- rep(0,length(col_pal))
  # col_pal <- 0:255 # paleta kolorow
  #iertujemy po palecie, czyli po kolejnych propozycjach wartosci pixela
  for (k in 1:length(col_pal)) {
    if (method == "potts") {
      prior_term <- sum(unlist(lapply(neighbors, function(coord) {
        return(as.integer(1 - (col_pal[k] == x[coord[1], coord[2]])))# sum(1-I(proposed_x == xij))
      })))
    }else {
      prior_term <- sum(unlist(lapply(neighbors, function(coord) {
        return(max((lambda*(col_pal[k]-x[coord[1], coord[2]]))^2,alpha))# sum(max((lambda*(proposed_x-xij))^2,alpha))
      })))
    }
    
    likelihood_term <-  (y[i, j] - col_pal[k] )^2 / (2 * sigma^2)# (y_ij - x_proposed)^2 / (2*sigma^2)
    conditional_distribution[k] <- exp(-beta * prior_term - likelihood_term)
  }
  conditional_distribution = conditional_distribution/sum(conditional_distribution)
  return(conditional_distribution)
}


calculate_conditional_distribution <- function(i,j,x,y,col_pal,beta,alpha, lambda,sigma, method, d) {
  neighbors <- get_neighbors(i,j, nrow(x), ncol(x), d)
  xs <- unlist(lapply(neighbors, function(coord) {
    return(x[coord[1], coord[2]])}))
  
  if (method == "potts") {
    prior_term <- sapply(col_pal, function(k) {
      sum(sapply(xs, function(x) {
        as.integer(1 - (k == x))
      }))
    })
    
  }else if(method == "2"){
    prior_term <- sapply(col_pal, function(k) {
      sum(sapply(xs, function(x) {
        max((lambda*(col_pal-x))^2,alpha)
      }))
    })
  }
  
  likelihood_term <-  (y[i, j] - col_pal)^2 / (2 * sigma^2)# (y_ij - x_proposed)^2 / (2*sigma^2)
  conditional_distribution <- exp(-beta * prior_term - likelihood_term)
  conditional_distribution <- conditional_distribution/sum(conditional_distribution)
  return(conditional_distribution)
}


denoise <- function(y, iter_n=1, beta=1,alpha=1, lambda=1,sigma=1, method = "potts", d=1) {
  x <- y
  list_x <- list()
  # paleta <- unique(c(y))
  paleta <- seq(0,1,length = 256)
  pb <- txtProgressBar(min = 0, max = iter_n*nrow(x)*ncol(x), style = 3)
  k  <- 0
  for (iter in 1:iter_n) {
    # cat(paste0("Iteracja ",iter), "\n")
    for (i in 1:nrow(x)) {
      if (i%%50==0) {
        # cat(paste0("Iteracja wiersza ",i), "\n")
      }
      for (j in 1:ncol(x)) {
        conditional_distribution <- calculate_conditional_distribution(i,j,x,y,paleta,beta,alpha,lambda,sigma,method, d)
        # conditional_distribution <- calculate_conditional_distribution2(i,j,x,y,paleta,beta,alpha,lambda,sigma,method, d) starsza wersja
        x[i,j] <- sample(paleta, size = 1, prob = conditional_distribution)
        k <-k+1
        setTxtProgressBar(pb, k)
      }
    }
    close(pb)
    list_x[[iter]] <- x
  }
  return(list_x)
}


img <- readPNG("lena_noisy.png")
denoised <- denoise(img[1:50,1:50],alpha=2,lambda=2,sigma=5,iter_n=2, d=2, method="2")

#srednia z ostanich N iteracji
# im_d <- Reduce("+",tail(denoised,N) )/N
im_d <- Reduce("+",tail(denoised,2) )/2

par(mfrow=c(1,2))
image(t(img[1:50,1:50])[,ncol(img[1:50,1:50]):1], axes = FALSE, col = grey(seq(0,1,length = 256)))
image(t(im_d)[,ncol(im_d):1], axes = FALSE, col = grey(seq(0,1,length = 256)))


