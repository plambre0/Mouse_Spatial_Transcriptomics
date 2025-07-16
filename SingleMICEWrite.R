library(mice)

set.seed(1922)

options(warn = -1)
data_path <- readline(prompt = "Enter the path to the expression data: ")
data <- read.csv(data_path, row.names = 1)

data_mice <- mice::mice(data,m=1)
data_complete <- mice::complete(data_mice)
filename <- readline(prompt = "Enter desired filename without extension: ")
write.csv(cat(filename,".csv"))