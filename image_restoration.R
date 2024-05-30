install.packages("png")
library("png")

img <- readPNG("C:/Users/kasia/OneDrive/Dokumenty/GitHub/GLASSO/lena_noisy.png")
image(t(img)[,ncol(img):1],col = c("black", "white"),axes=FALSE,useRaster=TRUE)


get_neighbors <- function(i, j, N, M) {
  # Define the potential neighbors
  potential_neighbors <- list(c(i-1, j), c(i+1, j), c(i, j-1), c(i, j+1))
  
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


calculate_conditional_distribution <- function(i,j,x,y,beta,sigma) {
  neighbors <- get_neighbors(i,j, nrow(x), ncol(x))
  conditional_distribution <- rep(0,2)
  col_pal <- c(0,1) # paleta kolorow
  for (k in 1:length(col_pal)) {
    prior_term <- sum(unlist(lapply(neighbors, function(coord) {
      return(as.integer(1 - (col_pal[k] == x[coord[1], coord[2]])))# sum(1-I(proposed_x == xij))
    })))
    
    likelihood_term <-  (y[i, j] - col_pal[k] )^2 / (2 * sigma^2)# (y_ij - x_proposed)^2 / (2*sigma^2)
    conditional_distribution[k] <- exp(-beta * prior_term - likelihood_term)
  }
  conditional_distribution = conditional_distribution/sum(conditional_distribution)
  return(conditional_distribution)
}


runif(1)
denoise <- function(y, beta, sigma, iter_n, method = "potts") {
  x <- y
  for (iter in 1:iter_n) {
    T = 1 - (iter-1)/iter_n
    cat(paste0("Iteracja ",iter), "\n")
    for (i in 1:nrow(x)) {
      if (i%%10==0) {
        cat(paste0("Iteracja wiersza ",i), "\n")
      }
      for (j in 1:ncol(x)) {
        conditional_distribution <- calculate_conditional_distribution(i,j,x,y,beta,sigma)
        new_pixel_value <- sample(0:1, size = 1, prob = conditional_distribution)
        energy_old <- -log(conditional_distribution[which(0:1 == x[i,j])])
        energy_new <- -log(conditional_distribution[which(0:1 == new_pixel_value)])
        delta_E <- energy_new -energy_old
        # symulowane wyzarzanie
        if (max(exp(-delta_E/T),1)> runif(1)) {
          x[i,j] <- new_pixel_value
        }
      }
    }
  }
  return(x)
}

denoised <- denoise(img,1,1,2)
par(mfrow=c(1,2))
image(t(img)[,ncol(img):1],col = c("black", "white"),axes=FALSE,useRaster=TRUE)
image(t(denoised)[,ncol(img):1],col = c("black", "white"),axes=FALSE,useRaster=TRUE)

