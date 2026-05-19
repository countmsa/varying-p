library(RcppAlgos)
library(numDeriv)

source("planC9_Bubble.R")

getab <- function(p=0.5) {
  a = p/(1-p) 
  c(a, 1)
}

Bubble5 = read.csv("Bubble5data.csv", header=TRUE)[,-6]

  

load(file="bubble5fit.Rdata")


F5 = fisherInfC.manyRX(tem5$par, nreps=5,  maxC=15, M=4, tol=1e-3, numDerv=FALSE) 
J5 = obsInfC.manyRX(tem5$par, cdata = as.matrix(Bubble5), M=4)
J5 = J5/48
#F5
#eigen(F5)$values

save(p5, par5, tem5, J5, F5, file="~/Downloads/Bubble5Fit/bubble5fit.Rdata")
#saveRDS(F5, file="F5.rds")

#F5 = readRDS("~/Downloads/Bubble5Fit/F5.rds")
