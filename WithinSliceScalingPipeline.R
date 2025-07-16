library(dplyr)
library(tidyr)
library(hash)
library(rgl)
library(vegan)

library(gstat)

#Bioconductor are required for variancePartition and preprocessCore.
#Use BiocManager::install to install both packages.
library(variancePartition)
library(preprocessCore)

options(warn = -1)

data_path <- readline(prompt = "Enter the path to the expression data: ")

planes_dict <- hash() 
planes_dict[["1"]] <- 'X'
planes_dict[["2"]] <- 'Z'
planes_dict[["3"]] <- 'Y'

cat("Select a plane:\n",
    "  1: Coronal\n",
    "  2: Sagittal\n",
    "  3: Axial\n")
slice <- readline('Enter corrosponding number: ')
plane <- planes_dict[[slice]]

data <- read.csv(data_path, row.names = 1)
data$Structure.ID <- as.factor(data$Structure.ID)

by_image <- split(data, data[[plane]])

normalized_images <- lapply(by_image, function(.x) {
  metadata <- .x %>% select(c('voxRowNum','X','Y','Z','Structure.ID'))
  genes <- .x %>% select(-c('voxRowNum','X','Y','Z','Structure.ID'))
  genes[genes==-1] <- NA
  normalized <- scale(genes)
  cbind(metadata, normalized)
})

processed_df <- do.call(rbind, normalized_images)
rownames(processed_df) <- NULL
processed_df <- processed_df %>% arrange(match(voxRowNum, data$voxRowNum))

genes <- t(as.matrix(data %>% select(-c('voxRowNum','X','Y','Z','Structure.ID'))))
genes_n <- t(as.matrix(processed_df %>% select(-c('voxRowNum','X','Y','Z','Structure.ID'))))
coords <- data %>% select(c('voxRowNum','X','Y','Z','Structure.ID'))
coords_n <- processed_df %>% select(c('voxRowNum','X','Y','Z','Structure.ID'))
rownames(coords) <- coords$voxRowNum
rownames(coords_n) <- coords_n$voxRowNum
colnames(genes) <- rownames(coords)
colnames(genes_n) <- rownames(coords_n)

genes <- apply(genes, 2, function(genes) {
  genes[is.na(genes)] <- median(genes, na.rm = TRUE)
  genes
})
genes_n <- apply(genes_n, 2, function(genes_n) {
  genes_n[is.na(genes_n)] <- median(genes_n, na.rm = TRUE)
  genes_n
})

write.csv(genes_n, paste0(basename(data_path),'_NormalizedMediamImp.csv'))

v_part <- variancePartition::fitExtractVarPartModel(
  exprObj = genes,
  formula = ~ X + Y + Z + Structure.ID,
  data = coords
)

print('Variance Partitioning results for non-normalized data:')
print(summary(v_part))

v_part_n <- variancePartition::fitExtractVarPartModel(
  exprObj = genes_n,
  formula = ~ X + Y + Z + Structure.ID,
  data = coords_n
)

print('Variance Partitioning results for normalized data:')
print(summary(v_part_n))

par(mfrow = c(2, 1))
print(plotVarPart(v_part) + ggtitle("Variance Partitioning for non-normalized data"))
print(plotVarPart(v_part_n) + ggtitle("Variance Partitioning for normalized data"))
