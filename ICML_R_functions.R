require(filling)
require(CVXR)
require(ggplot2)

fillDiag <- function(G, approx_method='auto') {
  ## This function fills the diagonal of G as in Algorithm 1
  
  ## For approx method you can pass either auto which will default
  ## to nuclear norm minimization for small matrix or SVT for larger
  ## or you can pass SVT or nuclear directly to force
  
  diag(G) <- NA
  if (approx_method == 'auto') {
    if (dim(G)[1] > 150) {
      new_diag <- diag(fill.SVT(G)$X)
    }
    if (dim(G)[1] <= 150) {
      new_diag <- diag(fill.nuclear(G)$X)
    }
  }
  
  if (approx_method == 'SVT') {
    new_diag <- diag(fill.SVT(G)$X)
  }
  
  if (approx_method == 'nuclear') {
    new_diag <- diag(fill.nuclear(G)$X)
  }
  
  diag(G) <- new_diag
  return(G)
}

computeAR <- function(G, rank=-1, compute_diag=TRUE) {
  ## It returns the rank k embeddings as well as various summary
  ## statistics
  
  ## Pass rank=-1 for full rank approx
  
  ## The compute_diag flag should be set to false if the diagonal
  ## is already pre-computed, otherwise it will replace by the nuclear
  ## norm minimizing diagonal
  
  if (compute_diag == TRUE) {
    G <- fillDiag(G)
  }
  
  eigs <- eigen(G)
  Vecs <- eigs$vectors
  
  D <- eigs$values
  ordering <- order(-abs(D))
  
  # To get low rank, truncate smallest eigenvalues to 0
  if (rank != -1) {
    D[ordering[-c(1:rank)]] <- 0
  }
  
  D_pos <- which(D > 0)
  D_neg <- which(D < 0)
  
  A <- 0
  Acomp <- 0
  An <- 0
  
  # Deal with R's terrible broadcasting
  
  ## Note that eigenvectors are column-wise in R
  if (length(D_pos) == 1) {
    A <- Vecs[,D_pos] * sqrt(D[D_pos])
    Acomp <- A %*% t(A)
    An <- A / sqrt(A^2)
  }
  
  if (length(D_neg) == 1) {
    R <- Vecs[,D_neg] * sqrt(abs(D[D_neg]))
    Rcomp <- R %*% t(R)
    Rn <- R / sqrt(R^2)
  }
  
  if (length(D_pos) > 1) {
    A <- as.matrix(Vecs[,D_pos]) %*% as.matrix(sqrt(diag(D[D_pos])))
    Acomp <- A %*% t(A)
    An <- A / sqrt(rowSums(A^2))
  }
  
  R <- 0
  Rcomp <- 0
  Rn <- 0
  if (length(D_neg) > 1) {
    R <-as.matrix(Vecs[,D_neg]) %*% as.matrix(sqrt(diag(abs(D[D_neg]))))
    Rcomp <- R %*% t(R)
    Rn <- R / sqrt(rowSums(R^2))
  }
  
  # Compute variance explained 
  G_approx <- Acomp - Rcomp
  
  delta <- G - G_approx
  print(delta)
  diag(delta) <- NA
  G_tmp <- G
  diag(G_tmp) <- NA
  var_explained <- 1 - mean(delta^2, na.rm=TRUE) /mean(G_tmp^2, na.rm=TRUE)
  R_fraction <- sum(D[D_neg]^2) / sum(D^2)
  
  results <- list('A' = A,
                  'R' = R,
                  'cos_sim_A' = An %*% t(An),
                  'cos_sim_R' = Rn %*% t(Rn),
                  'var_explained' = var_explained,
                  'R_fraction' = R_fraction)
  return(results)
}

# Test
raw <- matrix(rnorm(n=20*20), nrow=20, ncol=20)
G <- (raw + t(raw))/2

# Example: Compute full results for single matrix
r <- computeAR(G)

# Example: Compute variance explained for all ranks
G_filled <- fillDiag(G)
possible_ranks <- c(3:dim(G)[1])
vars <- c()
R_frac <- c()

for (k in possible_ranks) {
  tmp <- computeAR(G_filled, rank=k, compute_diag=FALSE)
  vars <- c(vars, tmp[['var_explained']])
  R_frac <- c(R_frac, tmp[['R_fraction']])
}

plot_df <- melt(data.frame(possible_ranks, vars, R_frac), id='possible_ranks')
ggplot(data=plot_df, aes(x=possible_ranks, y=value, colour=variable)) +
  geom_line() +
  xlab("Rank Chosen") +
  ylab("Value")