library(fields)
library(Matrix)

source("~/Documents/uvFAIMS/HE/read_data.R")

source("/home/jim/R/Library/MultiImage.R")

#' Load data
data.df = read_data()
X = scale(do.call(rbind, data.df$vecs)[,((52224/2)+1):52224], scale=FALSE)

X2 = Matrix(X, sparse=TRUE)

dim.faims = c(512, 51)
grid = make.surface.grid(list(x=seq(0, 1, length=dim.faims[1]),
                              y=seq(0, 1, length=dim.faims[2])))

keep = apply(X, 2, sd) > 0.006

X2[,!keep] = NA
#X2 = drop0(X2)

multiImage(grid, X2[1,], X2[2,], X2[3,], X2[4,], X2[5,])


