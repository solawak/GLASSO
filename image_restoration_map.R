library("png")
library(grid)
library(utils)

get_neighbors <- function(i, j, N, M, d) {
  # Definiujemy potencjalnych sasiadow
  if (d==1) {
    potential_neighbors <- list(c(i-1, j), c(i+1, j), c(i, j-1), c(i, j+1))
  }else {
    potential_neighbors <- list(c(i-1, j), c(i+1, j), c(i, j-1), c(i, j+1), c(i-1,j-1), c(i-1,j+1), c(i+1,j-1), c(i+1,j+1))
  }
  
  # Odfiltrowujemy sasiadow z poza macierzy
  valid_neighbors <- lapply(potential_neighbors, function(coord) {
    if (coord[1] >= 1 && coord[1] <= N && coord[2] >= 1 && coord[2] <= M) {
      return(coord)
    } else {
      return(NULL)
    }
  })
  valid_neighbors <- valid_neighbors[!sapply(valid_neighbors, is.null)]
  return(valid_neighbors)
}


calculate_conditional_distribution <- function(i,j,x,y,col_pal,beta,alpha, lambda,sigma, method, d) {
  neighbors <- get_neighbors(i,j, nrow(x), ncol(x), d)
  conditional_distribution <- rep(0,length(col_pal))
  #iterujemy po palecie, czyli po kolejnych propozycjach wartosci pixela
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

#wybierając metode:
#      "2" - należy usatwic parametry alpha,lambda i sigma
#      "potts" - należy usatwic parametry beta i sigma
denoise <- function(y, iter_n=1, beta=1,alpha=1, lambda=1,sigma=1, method = "potts", d=1) {
  x <- y
  paleta <- seq(0,1,length = 256)
  pb <- txtProgressBar(min = 0, max = iter_n*nrow(x)*ncol(x), style = 3)
  k  <- 0
  for (iter in 1:iter_n) {
    T = iter_n/(iter_n+1-iter)
    for (i in 1:nrow(x)) {
      for (j in 1:ncol(x)) {
        # wyznacz rozklad warunkowy
        conditional_distribution <- calculate_conditional_distribution(i,j,x,y,paleta,beta,alpha,lambda,sigma,method, d)
        new_pixel_value <- sample(paleta, size = 1, prob = conditional_distribution)
        delta_E <- conditional_distribution[which(paleta == x[i,j])] -conditional_distribution[which(paleta == new_pixel_value)]
        # symulowane wyzarzanie
        if (max(delta_E*exp(T),1)> runif(1)) {
          x[i,j] <- new_pixel_value
        }
        k <-k+1
        setTxtProgressBar(pb, k)
      }
    }
    close(pb)
  }
  return(x)
}

setwd("your_path")

img <- readPNG("image.png")
#wczytane obrazy (jezeli trzeba) przeksztalcamy do jednej macierzy z warosciami w [0,1]
img <- 0.2126*img[,,1] + 0.7152*img[,,2] + 0.0722*img[,,3]
denoised <- denoise(img[1:10,1:10],lambda=5,alpha=1.2,sigma=0.8,iter_n=1, d=2, method="potts")

#srednia z ostanich N iteracji
# im_d <- Reduce("+",tail(denoised,N) )/N
im_d <- Reduce("+",tail(denoised,5) )/5
par(mfrow=c(1,2))
image(t(img)[,ncol(img):1], axes = FALSE, col = grey(seq(0,1,length = 256)))
image(t(im_d)[,ncol(im_d):1], axes = FALSE, col = grey(seq(0,1,length = 256)))


img <- readPNG("lena_noisy.png")
im_d <- Reduce("+",tail(denoised,5) )/5
denoised <- denoise(img[1:10,1:10],alpha=0.1,lambda=0.1,sigma=1,iter_n=1, d=2, method="2")
image(t(img)[,ncol(img):1], axes = FALSE, col = grey(seq(0,1,length = 256)))
image(t(im_d)[,ncol(im_d):1], axes = FALSE, col = grey(seq(0,1,length = 256)))


