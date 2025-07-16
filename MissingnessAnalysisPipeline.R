library("dplyr")
library("tidyr")
library("rgl")
library("mgcv")

options(warn=-1)

data_path <- readline(prompt = "Enter the path to the expression data: ")
data <- read.csv(data_path, row.names = 1)
genes <- data %>% dplyr::select(-c(voxRowNum, X, Y, Z, Structure.ID))
coords <- data %>% dplyr::select(c(voxRowNum,X,Y,Z,Structure.ID))
genes[genes==-1] <- NA
genes[!is.na(genes)] <- 0
genes[is.na(genes)] <- 1
null_presence <- cbind(coords,genes)

null_count <- rowSums(genes)
colors <- heat.colors(length(null_count))[rank(null_count)]

rgl::plot3d(null_presence$Z, null_presence$X, null_presence$Y, col = colors, size=30)

alpha <- readline(prompt = "Enter alpha in decimal: ")

spatial_total_missingness <- summary(gam(rowSums(genes) ~ s(coords$X)+s(coords$Y)+s(coords$Z)))
print(spatial_total_missingness)

if(as.matrix(spatial_total_missingness$s.table)[1,4]<alpha){
  print('There is evidence X plane has significant relationship with missingness.')
} else{print('No evidence X plane has significant relationship with missingness.')}

if(as.matrix(spatial_total_missingness$s.table)[2,4]<alpha){
  print('There is evidence Y plane has significant relationship with missingness.')
} else{print('No evidence Y plane has significant relationship with missingness.')}

if(as.matrix(spatial_total_missingness$s.table)[3,4]<alpha){
  print('There is evidence Z plane has significant relationship with missingness.')
} else{print('No evidence Z plane has significant relationship with missingness.')}
