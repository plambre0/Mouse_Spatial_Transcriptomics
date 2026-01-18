#Required libraries
library(scatterplot3d)
library(spatstat)
library(bioimagetools)
library(igraph)

higher_k_clusters_neuron_pca <-
  read.csv("C:/Users/Paolo/Downloads/higher_k_clusters_neuron_pca.csv")

#Params
#Grid resolution
pixel_size_grid  <- 0.05
#Radius around each point we consider
max_distance_r   <- 0.5
#Threshold for weighted graph where we combine gene expression clusters.
clustering_threshold <- 0.01

#Study window bounds
x_limits <- c(5.094 - .001, 7.224 + .001)
y_limits <- c(3.688 - .001, 7.681 + .001)
z_limits <- c(10.33 - .001, 12.43 + .001)

study_volume <- diff(x_limits) * diff(y_limits) * diff(z_limits)

#Pre-Processing
higher_k_clusters_neuron_pca$cluster_id_k.145 <- as.factor(higher_k_clusters_neuron_pca$cluster_id_k.145)
cluster_labels <- levels(higher_k_clusters_neuron_pca$cluster_id_k.145)
n_clusters     <- length(cluster_labels)

#Ripley's-K calculation for all expression cluster pairs
cross_k_results <- list()

for (i in 1:(n_clusters - 1)) {
  for (j in (i + 1):n_clusters) {
    
    #Temp vars for each cluster pair
    data_i <- subset(higher_k_clusters_neuron_pca, cluster_id_k.145 == cluster_labels[i])
    data_j <- subset(higher_k_clusters_neuron_pca, cluster_id_k.145 == cluster_labels[j])
    
    # Skip if either cluster is too small to have a spatial distribution
    if (nrow(data_i) < 2 || nrow(data_j) < 2) next
    
    #K calculation
    k_obj <- bioimagetools::K.cross.3D(
      data_i$y_ccf, data_i$z_ccf, data_i$x_ccf,
      data_j$y_ccf, data_j$z_ccf, data_j$x_ccf,
      psz   = pixel_size_grid,
      width = max_distance_r
    )
    
    #Correct K based on global study volume vs local spread
    local_bbox_vol <- (diff(range(data_i$y_ccf)) * diff(range(data_i$z_ccf)) * diff(range(data_i$x_ccf)))
    if(local_bbox_vol > 0) {
      k_obj$y <- k_obj$y * (study_volume / local_bbox_vol)
    }
    
    #Store res
    pair_id <- paste0(cluster_labels[i], "_", cluster_labels[j])
    cross_k_results[[pair_id]] <- k_obj
  }
}

#L statistic calc. H(r) = L(r) - r
integrate_l_deviation <- function(k_result) {
  r_vec <- k_result$x
  k_vec <- k_result$y
  
  # Only process finite and positive dist.
  valid_idx <- is.finite(r_vec) & is.finite(k_vec) & r_vec > 0
  if (sum(valid_idx) < 2) return(0) 
  
  r_valid <- r_vec[valid_idx]
  k_valid <- k_vec[valid_idx]
  
  # L(r) transformation
  l_vec <- (3 * k_valid / (4 * pi))^(1/3)
  
  # H(r) is the deviation from CSR (Complete Spatial Randomness)
  h_vec <- l_vec - r_valid
  
  # Numerical Integration via Trapezoidal Rule
  # Calculates the area between the H(r) curve and the x-axis
  n <- length(r_valid)
  dx <- diff(r_valid)
  mean_heights <- (h_vec[-1] + h_vec[-n]) / 2
  
  total_area <- sum(mean_heights * dx)
  
  return(total_area)
}

#L matrix
adj_matrix_l <- matrix(0, nrow = n_clusters, ncol = n_clusters,
                       dimnames = list(cluster_labels, cluster_labels))

for (pair_name in names(cross_k_results)) {
  ids <- strsplit(pair_name, "_")[[1]]
  peak_val <- integrate_l_deviation(cross_k_results[[pair_name]])
  
  #Thresholding
  final_weight <- ifelse(peak_val > clustering_threshold, peak_val, 0)
  
  adj_matrix_l[ids[1], ids[2]] <- final_weight
  adj_matrix_l[ids[2], ids[1]] <- final_weight
}

#Weighted graph for spatial interaction 
g_clusters <- igraph::graph_from_adjacency_matrix(
  adj_matrix_l, 
  mode     = "undirected", 
  weighted = TRUE, 
  diag     = FALSE
)

#Community detection
louvain_comm <- igraph::cluster_louvain(g_clusters, weights = igraph::E(g_clusters)$weight)
meta_membership <- igraph::membership(louvain_comm)

#Assign meta-cluster IDs to the original df.
higher_k_clusters_neuron_pca$meta_cluster <- as.factor(
  meta_membership[as.character(higher_k_clusters_neuron_pca$cluster_id_k.145)]
)

#Visualization
scatterplot3d::scatterplot3d(
  higher_k_clusters_neuron_pca$y_ccf,
  higher_k_clusters_neuron_pca$z_ccf,
  higher_k_clusters_neuron_pca$x_ccf,
  color = as.numeric(higher_k_clusters_neuron_pca$meta_cluster),
  pch   = 16,
  xlab  = "Y", ylab = "Z", zlab = "X",
  main  = "Neuron Meta-Clusters (Based on Spatial L-Function)"
)
