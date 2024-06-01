library("png")
library(grid)
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


calculate_conditional_distribution <- function(i,j,x,y,col_pal,beta,alpha, lambda,sigma, method, d) {
  neighbors <- get_neighbors(i,j, nrow(x), ncol(x), d)
  conditional_distribution <- rep(0,length(col_pal))
  # col_pal <- 0:255 # paleta kolorow
  #iertujemy po palecie, czyli po kolejnych propozycjach wartości pixela
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


denoise <- function(y, iter_n=1, beta=1,alpha=1, lambda=1,sigma=1, method = "potts", d=1) {
  x <- y
  paleta <- unique(c(y))
  for (iter in 1:iter_n) {
    T = iter_n/(iter_n+1-iter)
    cat(paste0("Iteracja ",iter), "\n")
    for (i in 1:nrow(x)) {
      if (i%%50==0) {
        cat(paste0("Iteracja wiersza ",i), "\n")
      }
      for (j in 1:ncol(x)) {
        conditional_distribution <- calculate_conditional_distribution(i,j,x,y,paleta,beta,alpha,lambda,sigma,method, d)
        new_pixel_value <- sample(paleta, size = 1, prob = conditional_distribution)
        delta_E <- conditional_distribution[which(paleta == x[i,j])] -conditional_distribution[which(paleta == new_pixel_value)]
        # symulowane wyzarzanie
        if (max(delta*exp(T),1)> runif(1)) {
          x[i,j] <- new_pixel_value
        }
      }
    }
  }
  return(x)
}


img <- readPNG("lena_noisy.png")

denoised <- denoise(img,beta=1,sigma=1,iter_n=1, d=2)
par(mfrow=c(1,2))
image(t(img)[,ncol(img):1], axes = FALSE, col = grey(seq(0,1,length = 256)))
image(t(denoised)[,ncol(denoised):1], axes = FALSE, col = grey(seq(0,1,length = 256)))

#trzeba ustawic dobra lambde i alfe, lamda=1,alpha=1 zwraca szum XD
denoised2 <- denoise(img,alpha=2,lambda=5,sigma=1,iter_n=10, d=1, method="2")
par(mfrow=c(1,2))
image(t(img)[,ncol(img):1], axes = FALSE, col = grey(seq(0,1,length = 256)))
image(t(denoised2)[,ncol(denoised2):1], axes = FALSE, col = grey(seq(0,1,length = 256)))

png(filename = "denoised2.png", width = nrow(denoised2), height = ncol(denoised2))
grid.raster(denoised2, interpolate=FALSE)
dev.off()


