################################################################################
# Loading and Formatting
################################################################################

#replace with your directory
setwd("C:/Users/Paolo/Downloads/")

source("bivariate_ripley_k_3d.R")
library(readxl)
library(dplyr)
library(pheatmap)

cell_metadata <- read_excel("cell_metadata_rf_region_LRN_separate.xlsx")

higher_k_clusters_neuron_pca <- read.csv("higher_k_clusters_neuron_pca.csv")

metadata_joined <- higher_k_clusters_neuron_pca %>%
  select("cell_label", "x_ccf", "y_ccf", "z_ccf", "cluster_id_k.145") %>%
  left_join(cell_metadata %>% select("cell_label", 
                                     "parcellation_structure"), 
            by = 'cell_label')

summary(higher_k_clusters_neuron_pca$x_ccf)
summary(higher_k_clusters_neuron_pca$y_ccf)
summary(higher_k_clusters_neuron_pca$z_ccf)

cube_bounds <- c(
  min(higher_k_clusters_neuron_pca$x_ccf),
  max(higher_k_clusters_neuron_pca$x_ccf),
  min(higher_k_clusters_neuron_pca$y_ccf),
  max(higher_k_clusters_neuron_pca$y_ccf),
  min(higher_k_clusters_neuron_pca$z_ccf),
  max(higher_k_clusters_neuron_pca$z_ccf)
)

#Filter out clusters when n<30
cluster_ids <- names(table(higher_k_clusters_neuron_pca$cluster_id_k.145))
cluster_ids_large <- table(higher_k_clusters_neuron_pca$cluster_id_k.145)>30
cluster_ids_valid <- cluster_ids[cluster_ids_large]

neuron_pca_valid <- 
  higher_k_clusters_neuron_pca[higher_k_clusters_neuron_pca$cluster_id_k.145 
                               %in% 
                                 cluster_ids_valid,]

neuron_pca_valid_split <- 
  split(neuron_pca_valid, neuron_pca_valid$cluster_id_k.145)
neuron_pca_valid_split <- lapply(neuron_pca_valid_split, 
                                 function(x){cbind(x$x_ccf, x$y_ccf, x$z_ccf)})
metadata_filtered <- 
  metadata_joined[higher_k_clusters_neuron_pca$cluster_id_k.145 %in% 
                    cluster_ids_valid, ]

################################################################################
# Compute and Store Ripley's K Statistics
################################################################################

#Set sequence of r's for Ripley's K
r_values <- seq(0.05, 1.0, by = 0.05)

#Compute Bivariate K for all clusters
all_compare <- lapply(neuron_pca_valid_split, function(ref_points){
  lapply(neuron_pca_valid_split, function(eval_points_i){
    if(identical(ref_points, eval_points_i)){return(NULL)}
    bivariate_k_3d(ref_points, eval_points_i, r_values, cube_bounds)
  })
})

#Isolute structure cell coordinates
structs <- cell_metadata[, c("cell_label", 
                             "x_ccf", "y_ccf", "z_ccf", 
                             "parcellation_structure")]
structs <- split(structs, structs$parcellation_structure)
structs <- lapply(structs, function(x){cbind(x$x_ccf, x$y_ccf, x$z_ccf)})

#Compute Bivariate K for structures to clusters
results <- list()
for (struct_name in names(structs)){
  results[[struct_name]] <- list()
  for (split_name in names(neuron_pca_valid_split)){
    results[[struct_name]][[split_name]] <- 
      bivariate_k_3d(structs[[struct_name]],  
                     neuron_pca_valid_split[[split_name]], 
                     r_values, 
                     cube_bounds)
  }
}

#Store H for r=1
struct_clust_H <- lapply(results, function(x){
  lapply(x, function(x){
    h <- x$H
    h[length(h)]
  })
})
clust_clust_H <- lapply(all_compare, function(x){
  lapply(x, function(x){
    h <- x$H
    h[length(h)]
  })
})

################################################################################
# Construct Matrices, Visualize, Write to CSV
################################################################################
struct_clust_H_flat <- unlist(struct_clust_H)
struct_clust_H_matrix <- matrix(struct_clust_H_flat, 
                                nrow = length(struct_clust_H), 
                                byrow = TRUE)

rownames(struct_clust_H_matrix) <- names(struct_clust_H)
colnames(struct_clust_H_matrix) <- names(struct_clust_H[[1]])

pheatmap(struct_clust_H_matrix, 
         main = "H(r) for Clusters and Structure at r = 1")

clust_clust_H_flat <- sapply(clust_clust_H, function(x) {
  if (is.null(x)) NA else x
})

clust_clust_H_matrix <- matrix(
  clust_clust_H_flat,
  nrow = length(all_compare),
  byrow = TRUE
)

rownames(clust_clust_H_matrix) <- names(clust_clust_H_matrix)
colnames(clust_clust_H_matrix) <- names(clust_clust_H_matrix[[1]])

pheatmap(clust_clust_H_matrix, 
         main = "H(r) for Clusters at r = 1")

write.csv(struct_clust_H_matrix, "struct_clust_H_matrix.csv")
write.csv(clust_clust_H_matrix, "clust_clust_H_matrix.csv")