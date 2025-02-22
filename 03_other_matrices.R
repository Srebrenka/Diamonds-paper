################
# Diamonds  3  # 
################
# script 3

# Constructing other matrices 
# Inspecting their distributions and correlations



#########################
#                       #
#     Load libraries    #
#                       #
#########################
library(blockmodeling)
library(sna)
# install xUCINET from archive
library(xUCINET)
#library(devtools)
#install_github("jason-morgan/ina")
library(ina)
library(statnet)
library(intergraph)

##########################
# Saving all matrices
########################

G = c(nets_pos, rec_pos, vil_pos)
GN = c(nets_neg, rec_neg, vil_neg)
GS = c(nets, rec_net, villages)




# What is the highest neg degree of CN in ++, otherwise 
# median neg degree of all nodes with at least two positive ties


constructCNneg_matrix <- function(g, gn) {
  n <- vcount(g)  # Number of nodes in the network
  
  # Calculate degree of nodes in gn
  degree_gn <- igraph::degree(gn)
  
  # Find nodes in g with at least two ties
  nodes_with_two_ties <- which(igraph::degree(g) >= 2)
  
  # Calculate the median degree in gn for these nodes
  median_degree_gn <- median(degree_gn[nodes_with_two_ties])
  
  # Initialize the result matrix
  result_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Populate the matrix
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Find common neighbors in g
      common_neighbors <- intersect(neighbors(g, i), neighbors(g, j))
      
      if (length(common_neighbors) > 0) {
        # Get the degrees of common neighbors in gn
        common_degrees <- degree_gn[common_neighbors]
        # Assign the highest degree
        result_matrix[i, j] <- max(common_degrees, na.rm = TRUE)
      } else {
        # Assign the median degree in gn
        result_matrix[i, j] <- median_degree_gn
      }
      
      # Symmetry: Copy the value to the lower triangle of the matrix
      result_matrix[j, i] <- result_matrix[i, j]
    }
  }
  
  return(result_matrix)
}

Other_eff <- function(g, gn){

  degree_similarity_mat <- abs(outer(igraph::degree(g), igraph::degree(g), "-"))
  degree_g <- igraph::degree(g) # Degree in positive ties network
  degree_gn <- igraph::degree(gn) # Degree in negative ties network
  
  # Initialize Degree matrices
  degree_mat_g <- outer(degree_g, degree_g, "+")
  degree_mat_gn <- outer(degree_gn, degree_gn, "+")

  bet_mat_g <- outer(igraph::betweenness(g), igraph::betweenness(g))
  ###### FOR degree of common neighbor
  
  n <- vcount(g) 
  # Step 1: Calculate degree for all nodes
  degree_vals <- igraph::degree(g)
  
  # Step 2: Create a function to calculate the highest degree of common neighbors
  get_common_neighbors_degree <- function(g, u, v) {
    # Find common neighbors of nodes u and v
    common_neighbors <- intersect(neighbors(g, u), neighbors(g, v))
    
    if (length(common_neighbors) > 0) {
      # Return the highest degree from common neighbors
      return(max(degree_vals[common_neighbors]))
    } else {
      # If no common neighbors, calculate the median degree for nodes with at least two ties
      nodes_with_at_least_two_ties <- V(g)[degree_vals >= 2]
      return(median(degree_vals[nodes_with_at_least_two_ties]))
    }
  }
  
  # Step 3: Construct the degree matrix

  degreeCN_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Loop through all pairs of nodes
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      degreeCN_matrix[i, j] <- get_common_neighbors_degree(g, i, j)
      degreeCN_matrix[j, i] <- degreeCN_matrix[i, j]  # Symmetric matrix
    }
  }
  
  
  # Step 1: Calculate betweenness for all nodes
  betweenness_vals <- igraph::betweenness(g)
  
  # Step 2: Create a function to calculate the highest betweenness of common neighbors
  get_common_neighbors_betweenness <- function(g, u, v) {
    # Find common neighbors of nodes u and v
    common_neighbors <- intersect(neighbors(g, u), neighbors(g, v))
    
    if (length(common_neighbors) > 0) {
      # Return the highest betweenness from common neighbors
      return(max(betweenness_vals[common_neighbors]))
    } else {
      # If no common neighbors, calculate the median betweenness for nodes with at least two ties
      nodes_with_at_least_two_ties <- V(g)[degree_vals >= 2]
      return(median(betweenness_vals[nodes_with_at_least_two_ties]))
    }
  }
  
  # Step 3: Construct the betweenness matrix
  betCN_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Loop through all pairs of nodes
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      betCN_matrix[i, j] <- get_common_neighbors_betweenness(g, i, j)
      betCN_matrix[j, i] <- betCN_matrix[i, j]  # Symmetric matrix
    }
  }

  # CN negative degree
  negDegCN_matrix <- constructCNneg_matrix(g, gn)
  

  out2 <- list(degree_similarity_mat, degree_mat_g,
               degree_mat_gn,
               bet_mat_g,
               degreeCN_matrix,  
               betCN_matrix, negDegCN_matrix)
  names(out2) <- c( "degree_similarity", "degree_pos", "degree_neg",
                   "bet_mat_pos", "degreeCN", "betCN", "negDegCN")
  return(out2)
}

other_effects <- c( "degree_similarity", "degree_pos", "degree_neg",
                    "bet_mat_pos", "degreeCN", "betCN", "negDegCN")
outcomes

OE <- list() # there will be 7 in each list

for(z in 1:length(G)){
  OE[[z]] <- Other_eff(G[[z]], GN[[z]])
  print(z)
}

DVs <- list() # 2
for(z in 1:length(G)){
  DVs[[z]] <- open_triads_mats(G[[z]], GN[[z]])
  print(z)
}
length(DVs[[1]])
# run 
Main.mats<- list() # 9
for(z in 1:length(G)){
  g <- G[[z]]
  gs <-GS[[z]]
  gn <- GN[[z]]
  Main.mats[[z]]<- dyad_based2(g=g, gs= gs, gn = gn)
  print(z)
}


# Other diamonds
# sum matrices 8 and 9 in main nets

otherD <- list() # 1
for(z in 1:length(G)){
  otherD[[z]] <- Main.mats[[z]][[8]] + Main.mats[[z]][[9]]
}
length(otherD)
class(otherD[[2]])
##### 
# other common neighbours? ++ Npp - 1
# binarize them

NC_PPminus1<- list() # 1
for(z in 1:length(G)){
  mat <- Main.mats[[z]][[4]] - 1
  mat[mat < 0] <- 0
  NC_PPminus1[[z]] <- mat
}

Any_t_OLD <- Any_t_b

Any_t_b <- list() # it is only binary measure
for(z in 1:length(G)){
  mat <- MATS[[z]][[7]] + MATS[[z]][[8]] + MATS[[z]][[9]]
  mat[mat > 1] <- 1
  Any_t_b[[z]] <- mat
}


# 20 per each? - 21 later
MATS <- list()
for(z in 1:24){
  
  MATS[[z]] <- c(Main.mats[[z]], OE[[z]], 
                 list(otherD[[z]]), list(NC_PPminus1[[z]]),
                 DVs[[z]], list(Any_t_b[[z]])) 
  
}

Names_matrices2

mx <- MATS

on <- mx[[1]]

sum(!(on[[20]] == Any_t_OLD[[1]]))

length(mx)

Names_matrices2 <- c(Names_matrices, "any_tet")


saveRDS(Any_t_b, file = "Any_tetrad_b.rds")

max(Any_t_b[[6]])

length(MATS)
# saving
saveRDS(MATS, file = "MATS_d.rds") # saved 25.12

MATS <- readRDS("MATS_d.rds")


# sanity checks
test <- MATS[[15]]
for(i in 1:length(test)){
  print(i)
  mat <- test[[1]]
  # Check if the matrix is numeric
  is_numeric <- is.numeric(mat)
  print(is_numeric)
  # Check for infinite values
  has_infinite <- any(is.infinite(mat))
  print(has_infinite)
  # Check for NA values
  has_na <- any(is.na(mat))
  print(has_na)
  is_mat <- is.matrix(mat)
  # Check for negative values
  has_negative <- any(mat < 0, na.rm = TRUE)
  print(has_negative)
  
  
}

# Descriptive - occurence 

# Function to calculate the percentage of non-zero cells in the upper triangle
calculate_upper_triangle_percentage <- function(mat_list) {
  sapply(mat_list, function(mat) {
    # Extract the upper triangular part (excluding diagonal)
    upper_triangle <- mat[upper.tri(mat)]
    
    # Count the number of non-zero elements
    non_zero_count <- sum(upper_triangle == 0)
    
    # Calculate the total number of elements in the upper triangle
    total_count <- length(upper_triangle)
    
    # Calculate the percentage of non-zero cells
    percentage <- (non_zero_count / total_count) * 100
    return(percentage)
  })
}

# Apply the function to the list of matrices


Names_matrices <- c(outcomes, other_effects, "other diamonds", 
                    "NCpp_minus1", "open_t", "true_open_T", "any_tetrad")

percM <- list()
for(z in 1:length(MATS)){
  ML <- MATS[[z]]
  percentages <- calculate_upper_triangle_percentage(ML)
  percM[[z]]<- as.vector(round(percentages,1))
}

pml <- as.data.frame(do.call(rbind, percM))
colnames(pml) <- Names_matrices
pml$Nets <- NamesNets

write.xlsx(pml, "Distribution_of_d_matrices_values.xlsx") 



#####################################
# QAP correlations between matrices
#######################################



CMAT1 <- list()
PMAT1 <- list()

for(z in 1:length(G)){
  print(z)
  MATSx <- mx[[z]]
  
  for (k in 1:length(MATSx)){
    if (k %in% c(4, 5, 6, 21)){
      tmat <- MATSx[[k]]
      tmat[tmat > 1] <- 1
      MATSx[[k]] <-tmat
    }
 
  }

  VarN1 <- Names_matrices2
  cmat1 <- matrix(NA, nrow = length(VarN1), ncol = length(VarN1))
  pmat1 <- matrix(NA, nrow = length(VarN1), ncol = length(VarN1))
  colnames(cmat1) <- VarN1 
  rownames(cmat1) <- VarN1
  
  # calculate via matrix
  
  for (i in 1:(length(VarN1) - 1)) {       # Loop over rows, excluding the last row
    for (j in (i + 1):length(VarN1)) {     # Loop over columns starting from the element after the diagonal
      # Perform your operations on M[i, j]
      print(c(i,j))
      if(sum(MATSx[[i]])<2 | sum(MATSx[[j]])<2){
        cmat1[i, j] <-NA
        pmat1[i, j] <-NA
      }else{
        k <- xQAPCorrelation(MATSx[[i]], MATSx[[j]], NPERM = 1000, Directed = F, 
                             Loops = FALSE)
        cmat1[i, j] <- round(k$CorrCoef,2)
        pmat1[i, j] <- round(k$ClassicSign_2tailed,3)
        
      }}}
  CMAT1[[z]]<- cmat1
  PMAT1[[z]]<- pmat1
  print(paste("end_", z, sep=""))
}
saveRDS(CMAT1, file = "model_any_variables_qap_cor.rds")
saveRDS(PMAT1, file = "model_any_variables_qap_p_values.rds")
CMAT1[[1]]
### Summarize in on (two) matrices

# Create empty matrices to store the results
mean_cmat <- matrix(NA, nrow = 21, ncol = 21)
colnames(mean_cmat) <- Names_matrices2 # MEDIAN!!!
rownames(mean_cmat)<- Names_matrices2
max_abs_cmat <- matrix(NA, nrow = 21, ncol = 21)
colnames(max_abs_cmat) <- Names_matrices2
rownames(max_abs_cmat)<- Names_matrices2
min_abs_cmat <- matrix(NA, nrow = 21, ncol = 21)
colnames(min_abs_cmat) <- Names_matrices2
rownames(min_abs_cmat)<- Names_matrices2

# Loop over the upper triangle (i < j) to compute mean, max and min
for (i in 1:20) {  # Looping until 13 since 14th row has no upper triangle cells
  for (j in (i + 1):21) {  # Loop over columns, always above the diagonal
    # Extract the values from each matrix in the list for this cell [i,j]
    values <- sapply(CMAT1, function(mat) mat[i, j])
    
    # Remove NAs
    values <- values[!is.na(values)]
    
    if (length(values) > 0) {
      # Calculate the mean for this cell
      mean_cmat[i, j] <- round(median(values),2)  ## MEDIAN
      
      # Calculate the max absolute value for this cell
      max_abs_cmat[i, j] <- round(max(values),2)
      
      # Calculate the min absolute value for this cell
      min_abs_cmat[i, j] <- round(min(values),2)
    }
  }
}

n_pmat <- matrix(NA, nrow = 21, ncol = 21)
colnames(n_pmat) <- Names_matrices2
rownames(n_pmat)<- Names_matrices2

# Loop over the upper triangle (i < j) to compute mean, max and min
for (i in 1:20) {  # Looping until 13 since 14th row has no upper triangle cells
  for (j in (i + 1):21) {  # Loop over columns, always above the diagonal
    # Extract the values from each matrix in the list for this cell [i,j]
    values <- sapply(PMAT1, function(mat) mat[i, j])
    
    # Remove NAs
    values <- round(values[!is.na(values)],2)
    
    if (length(values) > 0) {
      # Calculate the mean for this cell
      n_pmat[i, j] <- count <- sum( values <= 0.05)
    }
  }
}

mean_cmatdf <- as.data.frame(mean_cmat)

write.xlsx(mean_cmatdf, "Qap_corr_of_matrices_Median.xlsx", rowNames =T)

n_pmatdf <- as.data.frame(n_pmat)
write.xlsx(n_pmatdf, "N_sig_pvalues_qap_corr_of_matrices.xlsx", rowNames =T)


# 1 means all correlation between -1 and 1 happened

range_of_cor <- round(abs(max_abs_cmat - min_abs_cmat)/2,2)
range_of_cordf<- as.data.frame(range_of_cor)
write.xlsx(range_of_cordf, "Range_qap_corr_of_matrices.xlsx", rowNames =T)


