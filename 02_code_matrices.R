
################
# Diamonds 2   # 
################
# script 2


#########################
#                       #
#     Load libraries    #
#                       #
#########################




# Input network (undirected)

# Define input matrices: pos, neg, or/and both

# Define matrices IVs



dyad_based <- function(g, gs, gn){
  
  # positive net
  adjg <- as_adj(g)
  adjm <- as.matrix(adjg)
  # negat net
  adjgn <- as_adj(gn)
  adjmn <- as.matrix(adjgn)
  # signed
  adjg <- as_adj(gs)
  
  adjms <- as.matrix(adjg)
  
  # R mat has 1 in a cell if any tie, pos or neg
  R_mat <- matrix(0, nrow = vcount(g), ncol = vcount(g))
  Ncpp_mat <- matrix(0, nrow = vcount(g), ncol = vcount(g))
  
  Ncnn_mat <- matrix(0, nrow = vcount(g), ncol = vcount(g))
  Ncn_mat <- matrix(0, nrow = vcount(g), ncol = vcount(g))
  
  for (i in 1:vcount(g)) { # counting it twice ij and ji
    for (j in 1:vcount(g)) {
      if(adjm[i,j] == 1 | adjmn[i,j] == 1){R_mat[i,j] = 1}
      else{R_mat[i,j]=0}
      Ncpp_mat[i,j] <- length(intersect(neighbors(g, i), neighbors(g, j))) #~ number of pos cn
      Ncnn_mat[i,j] <- length(intersect(neighbors(gn, i), neighbors(gn, j)))
      Ncn_mat[i,j] <- length(intersect(neighbors(gs, i), neighbors(gs, j)))
      ## this is important:
      diag(Ncpp_mat) = 0
      diag(Ncnn_mat) = 0
      diag(Ncn_mat) = 0
    }
  }
  # mixed +/_  common neighbour of the two nodes
  pn_mat<- matrix(0, nrow = vcount(g), ncol = vcount(g))
  for (i in 1:vcount(g)) { # counting it twice ij and ji
    for (j in 1:vcount(g)) {
      al <- Ncpp_mat[i,j] + Ncnn_mat[i,j] # number of ++ and -- cn
      if (#Ncpp_mat[i,j]!=0 & # here is the mistake
          Ncn_mat[i,j] > al
      ){
        pn_mat[i,j] <- Ncn_mat[i,j] - al # number of =/- cn
        
      }else{pn_mat[i,j] <-0}
    }
  }
  
  d_mat <-matrix(0, nrow = vcount(g), ncol = vcount(g))
  d4p_mat <- matrix(0, nrow = vcount(g), ncol = vcount(g))
  d2p2n_mat <- matrix(0, nrow = vcount(g), ncol = vcount(g))
  for (i in 1:vcount(g)) { # counting it twice ij and ji
    for (j in 1:vcount(g)) {
      for (i in 1:vcount(g)) { # counting it twice ij and ji
        for (j in 1:vcount(g)) {
          if (Ncpp_mat[i,j] != 0 & pn_mat[i,j] != 0) {
            d_mat[i,j] = 1
          }
          if (Ncpp_mat[i,j] >= 2) {
            d4p_mat[i,j] = 1
          }
          if (Ncpp_mat[i,j] != 0 & Ncnn_mat[i,j] != 0) {
            d2p2n_mat[i,j] = 1
          }
        }
      }
      }
      }

  out = list(R_mat, adjm, adjmn, Ncpp_mat, Ncnn_mat, pn_mat,
             d_mat, d4p_mat, d2p2n_mat) # 9
  # DV: any tie, pos tie, negative tie
  # IV: ++ CN, -- CN, +- CN  -how many common neighbors
  # diamond 3p1n, diamond 4p, diamond 2p2n - binary variables
  names(out) <- c("any_tie", # any tie, pos or neg ==1
                  "pos_tie", "neg_tie", # just matrices 
                  "Nc_pp", "Nc_nn", "Nc_pn", # N common neigh
                  "d_3pn", "d_4p", "d_2p2n") # binary, with AC tie or not
  return(out)
}

# Names of outcomes of the function above
outcomes <- c("any_tie", # any tie, pos or neg ==1
              "pos_tie", "neg_tie", # just matrices 
              "Nc_pp", "Nc_nn", "Nc_pn", # N common neigh
              "d_3pn", "d_4p", "d_2p2n") # binary, with AC tie or not


# any_tie or pos network should be DV, neg network DV

dyad_based2 <- function(g, gs, gn) {
  # Convert adjacency matrices to dense matrices
  adjm <- as.matrix(as_adj(g))
  adjmn <- as.matrix(as_adj(gn))
  adjms <- as.matrix(as_adj(gs))
  
  # Compute common neighbors using matrix multiplications
  Ncpp_mat <- adjm %*% t(adjm)
  Ncnn_mat <- adjmn %*% t(adjmn)
  Ncn_mat <- adjms %*% t(adjms)
  
  # Remove diagonals
  diag(Ncpp_mat) <- 0
  diag(Ncnn_mat) <- 0
  diag(Ncn_mat) <- 0
  
  # R_mat: Any tie (positive or negative)
  R_mat <- (adjm | adjmn) * 1
  
  # Mixed +- common neighbors
  pn_mat <- pmax(Ncn_mat - (Ncpp_mat + Ncnn_mat), 0)
  
  # Diamond matrices
  d_mat <- (Ncpp_mat != 0 & pn_mat != 0) * 1
  d4p_mat <- (Ncpp_mat >= 2) * 1
  d2p2n_mat <- (Ncpp_mat != 0 & Ncnn_mat != 0) * 1
  
  # Output
  out <- list(
    any_tie = R_mat, pos_tie = adjm, neg_tie = adjmn, 
    Nc_pp = Ncpp_mat, Nc_nn = Ncnn_mat, Nc_pn = pn_mat,
    d_3pn = d_mat, d4pn = d4p_mat, d_2p2n = d2p2n_mat
  )
  
  return(out)
}


########################################################
# Matrices of OT and TOTs (DV)
#######################################################

open_triads_mats <- function(g, gn){
  # Number of nodes
  n <- vcount(g)
  
  # Initialize OTmat and TOTmat as zero matrices
  OTmat <- matrix(0, n, n)
  TOTmat <- matrix(0, n, n)
  
  # Get adjacency matrices of the graphs
  adj_g <- as_adj(g, sparse = FALSE, type = "both")   # Adjacency matrix of positive ties network
  adj_gn <- as_adj(gn, sparse = FALSE, type = "both") # Adjacency matrix of negative ties network
  
  # Check pairs of nodes
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Check if there is no direct tie between nodes i and j in g
      if (adj_g[i, j] == 0) {
        # Find common neighbors in g
        common_neighbors <- which(adj_g[i, ] == 1 & adj_g[j, ] == 1)
        
        # OTmat: 1 if nodes i and j have a common neighbor in g
        if (length(common_neighbors) > 0) {
          OTmat[i, j] <- 1
          OTmat[j, i] <- 1
        }
        
        # TOTmat: 1 if nodes i and j have a common neighbor in g AND no tie in gn
        if (length(common_neighbors) > 0 && adj_gn[i, j] == 0) {
          TOTmat[i, j] <- 1
          TOTmat[j, i] <- 1
        }
      }
    }
  }
  
  # Result
  res1 <- list(OTmat, TOTmat)
  names(res1) <- c("open_triads", "true_open_triads")
  return(res1)
}








####################
# END of the script
####################




