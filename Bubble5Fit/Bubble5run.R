library(RcppAlgos)
library(numDeriv)

source("count_helper_functions_bubble.R")
source("planC_Bubble.R")


getab <- function(p=0.5) {
  a = p/(1-p) 
  c(a, 1)
}

Bubble5 = read.csv("~/Downloads/Bubble5Fit/Bubble5data.csv", header=TRUE)[,-6]
head( Bubble5 )


p5 = EM(ipar=c(mean(as.matrix(Bubble5)), 0.9, 0.1), cdata=as.matrix(Bubble5), maxN=30, maxiter=10^3, atol=1e-3) 

par5 = c(p5$mle[c(1,3)],  getab(p5$mle[2]))


fn5 = create.Cloglik.manyRX(cdat=as.matrix(Bubble5), M=6)
#fn5(par5)


system.time({ tem5 = optim( par=par5, fn=fn5, control=list(fnscale=-1, maxit=1000) )  })

## p value
pchisq(2*(tem5$value - max(p5$loglik)), 1, lower.tail = FALSE) 


#F5 = fisherInfC.manyRX(tem5$par, nreps=5,  maxC=15, M=4, tol=1e-3, numDerv=FALSE) 

round(sqrt(var.tau(fisherInfpar=F5*48, par=tem5$par) ),3)

#save(p5, par5, tem5, J5, F5, file="bubble5fit.Rdata")








alpha_V = function(par0 = NULL){
  l = par0[1]
  t = par0[2]
  a = par0[3]
  b = par0[4]
  mu = a / (a+b)
  l* mu^2 * (( a + 1 )*(a + b)/ (a*(a + b + 1) )) / (l * mu + t)
}

tau_derivs=function(par=NULL ){
  l =par[1]; t = par[2];a=par[3];b=par[4]
  mu = a/(a+b)
  l_der = a ^ 2 / (a + b) ^ 2 * (1 + b / a / (a + b + 1)) / (l * a / (a + b) + t) - l * a ^ 3 / 
    (a + b) ^ 3 * (1 + b / a / (a + b + 1)) / (l * a / (a + b) + t) ^ 2
  t_der = -l * a ^ 2 / (a + b) ^ 2 * (1 + b / a / (a + b + 1)) / (l * a / (a + b) + t) ^ 2
  a_der = 2 * l * a / (a + b) ^ 2 * (1 + b / a / (a + b + 1)) / (l * a / (a + b) + t) - 
    2 * l * a ^ 2 / (a + b) ^ 3 * (1 + b / a / (a + b + 1)) / (l * a / (a + b) + t) + 
    l * a ^ 2 / (a + b) ^ 2 * (-b / a ^ 2 / (a + b + 1) - b / a / (a + b + 1) ^ 2) / 
    (l * a / (a + b) + t) - l * a ^ 2 / (a + b) ^ 2 * (1 + b / a / (a + b + 1)) / 
    (l * a / (a + b) + t) ^ 2 * (l / (a + b) - l * a / (a + b) ^ 2)
  b_der = -2 * l * a ^ 2 / (a + b) ^ 3 * (1 + b / a / (a + b + 1)) / 
    (l * a / (a + b) + t) + l * a ^ 2 / (a + b) ^ 2 * (1 / a / (a + b + 1) - b / a / (a + b + 1) ^ 2) / 
    (l * a / (a + b) + t) + l ^ 2 * a ^ 3 / (a + b) ^ 4 * (1 + b / a / (a + b + 1)) / (l * a / (a + b) + t) ^ 2
  return(c(l_der, t_der, a_der, b_der))
}

var.tau <- function(fisherInfpar=NULL, par=NULL) {
  dpar = tau_derivs(par)
  
  as.numeric( t(dpar) %*% solve(fisherInfpar) %*% (dpar) )
}







