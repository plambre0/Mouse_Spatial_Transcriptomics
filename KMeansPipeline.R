library(dplyr)
library(tidyverse)
library(rgl)
library(hash)

set.seed(1922)
options(warn=-1)

data_path <- readline(prompt = "Enter the path to the expression data: ")
data <- read.csv(data_path, row.names = 1)

planes_dict <- hash() 
planes_dict[["1"]] <- 'X'
planes_dict[["2"]] <- 'Z'
planes_dict[["3"]] <- 'Y'

slice_dict <- hash() 
slice_dict[["X"]] <- 'Coronal'
slice_dict[["Z"]] <- 'Sagittal'
slice_dict[["Y"]] <- 'Axial'

cat("Select a plane:\n",
    "  1: Coronal\n",
    "  2: Sagittal\n",
    "  3: Axial\n")
slice <- readline('Enter corrosponding number: ')
plane <- planes_dict[[slice]]

k <- readline('Enter k for k-means clustering: ')

by_image <- split(data, data[[plane]])

normalized_images <- lapply(by_image, function(.x) {
  metadata <- .x[, c('voxRowNum','X','Y','Z','Structure.ID')]
  genes <- .x[, !names(.x) %in% c('voxRowNum','X','Y','Z','Structure.ID')]
  
  genes[genes == -1] <- NA
  normalized <- scale(genes)
  cbind(metadata, normalized)
})

processed_df <- do.call(rbind, normalized_images)
rownames(processed_df) <- NULL

genes_raw <- data %>%
  select(-voxRowNum, -X, -Y, -Z, -Structure.ID) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
genes_normalized <- processed_df %>%
  select(-voxRowNum, -X, -Y, -Z, -Structure.ID) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

kmeans_raw <- kmeans(genes_raw,as.integer(k))
kmeans_normalized <- kmeans(genes_normalized,as.integer(k))

raw_clustered <- cbind(cluster = kmeans_raw$cluster,data[,c("voxRowNum", "Structure.ID", "X", "Y", "Z")],genes_raw)
norm_clustered <- cbind(cluster = kmeans_normalized$cluster,data[,c("voxRowNum", "Structure.ID", "X", "Y", "Z")],genes_normalized)
write.csv(raw_clustered, paste(slice_dict[[plane]],'_raw_clusters.csv',sep = ''))
write.csv(norm_clustered, paste(slice_dict[[plane]],'_norm_clusters.csv',sep = ''))

plot3D::scatter3D(data$Z,data$X,data$Y,col=kmeans_raw$cluster, colkey = FALSE, xlab = 'Z', ylab = 'X', zlab = 'Y')
plot3D::scatter3D(processed_df$Z,processed_df$X,processed_df$Y,col=kmeans_normalized$cluster, colkey = FALSE, xlab = 'Z', ylab = 'X', zlab = 'Y')

cat("Cluster results for ", slice_dict[[plane]], " Slice by percent similarity:","\n")
cat("Structure ID's by clusters for raw data",'\n')
results_raw <- table(Structures = data$Structure.ID,Clusters = kmeans_raw$cluster)
print(results_raw)
cat("Structure ID's by clusters for normalized data",'\n')
results_norm <- table(Structures = data$Structure.ID,Clusters = kmeans_normalized$cluster)
print(results_norm)

write.csv(results_raw,"Cluster_Struct_Table_Raw_Med_Imp.csv")
write.csv(results_norm,"Cluster_Struct_Table_Norm_Med_Imp.csv")
