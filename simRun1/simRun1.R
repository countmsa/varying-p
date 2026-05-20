source("simRun1_functions.R")
source("fisher.R")

### SIM PILOT

### Plan Order Y, NF, F, N, C
mle    = matrix(NA, 5,4)
loglik = matrix(NA, 5,2)
time =  numeric(5)
asym.sd = matrix(NA,5,4)
tau.sd = numeric(5)

par0=c(3, 0.1, 9, 3)
nreps = 4
nparts = 50

#cdat.full = rPPdata(nparts, nreps, par0)

load(file="simRun1Data.Rdata")
cdat = fdat + xdat


### STARTING VALUES
#lambda, p ,theta
em.temp = EM(ipar=c(mean(cdat), 0.9,0.1), cdata=cdat, maxiter=1000, atol=1e-4, maxN = 30)
EM.start = em.temp$mle

st.theta = EM.start[3]
st.lam = EM.start[1]
ms = EM.start[2]
as = (ms*(1-ms) - 0.01)/(0.01*(1 - (ms-1)/ms))
bs = max(as*(1-ms)/ms, 1.05)
mndat = mean(cdat)
mfdat = mean(fdat)
      
#### RUN the EM  -  Plan Order Y, NF, F, N, C  
#### it will take under 30 seconds to run all 5 EMs
time[1] = system.time({ tempY = EMY.ab(ab=c(mndat,mfdat,as,bs), zdat=zdat, ndat=ndat, fdat=fdat, nreps=nreps,maxiter=1000, atol=1e-2)}, gcFirst=TRUE)[3]
mle[1,]     = tempY$mle
loglik[1,1] = tempY$loglik
      
time[2] = system.time({ tempNF = EMX.ab(par=c(mndat,mfdat,as,bs), xdat=xdat , fdat=fdat, ndat=ndat, maxiter=1000, atol=1e-2) }, gcFirst=TRUE)[3]
mle[2,]     = tempNF$mle
loglik[2,1] = tempNF$loglik
 
time[4] = system.time({tempN = EMN(par=c(mndat,st.theta,as,bs), cdat=cdat, ndat = ndat, maxiter=1000,  atol=1e-2)}, gcFirst=TRUE)[3]
mle[4,]     = tempN$mle
loglik[4,1] = tempN$loglik
      
time[3] = system.time({ tempF = EMR(par=c(st.lam, st.theta, as,bs), xdat=xdat, fdat=fdat, M=4, maxiter=1000, atol=1e-2)}, gcFirst=TRUE)[3]
mle[3,]= tempF$mle
loglik[3,1] = tempF$loglik
      
time[5] = system.time({ tempC = EMC(par=c(st.lam, st.theta, as,bs), cdat=cdat, atol=1e-2, maxiter = 1000)}, gcFirst=TRUE)[3]
mle[5,]= tempC$mle
loglik[5,1] = tempC$loglik

print(mle)
      
######the fixed cases
n.lik = create.loglikn2(cdat,ndat)
loglik[4,2] = -optim(c(3,0.9,0.1), n.lik)$value
loglik[3,2]= EM.r(ipar=c(st.lam, 0.9,0.1), xdat=xdat, fdat=fdat, maxiter=100, atol=1e-4)$loglik
      
#GS
gs.lik = create.loglikGS(xdat, fdat,ndat)
mle.par = c(mean(ndat), mean(xdat)/mean(ndat), mean(fdat))
loglik[2,2] = gs.lik(mle.par)
loglik[5,2] = em.temp$loglik
      

##Asymptotic standard deviations for each parameter and tau
##It will take approx. 30 seconds - 1 minute for all Fisher informations to complete
FY = fisherInfY.ab(mle[1,], nreps=4)
asym.sd[1,] = sqrt(diag(solve(nparts*FY)))
tau.sd[1] = sqrt(var.tau(FY*nparts, par = mle[1,]))

FNF = fisherInfX2(mle[2,] ,maxX=12, nReps=4, M=5 )
asym.sd[2,] = sqrt(diag(solve(nparts*FNF)))
tau.sd[2] = sqrt(var.tau(FNF*nparts, par = mle[2,]))

FF = fisherInfR(mle[3,], nreps=4)
asym.sd[3,] = sqrt(diag(solve(nparts*FF)))
tau.sd[3] = sqrt(var.tau(FF*nparts, par = mle[3,]))

FN = fisherInfN(mle[4,], nreps=4)
asym.sd[4,] = sqrt(diag(solve(nparts*FN)))
tau.sd[4] = sqrt(var.tau(FN*nparts, par = mle[4,]))

FC = fisherInfC(mle[5,], nreps=nreps)
asym.sd[5,] = sqrt(diag(solve(nparts*FC)))
tau.sd[5] = sqrt(var.tau(FC*nparts, par = mle[5,]))

rownames(asym.sd) = c("Plan Y", "Plan NF", "Plan F", "Plan N", "Plan C")
colnames(asym.sd) = c("lambda", "theta", "alpha","beta")
rownames(mle) = c("Plan Y", "Plan NF", "Plan F", "Plan N", "Plan C")
colnames(mle) = c("lambda", "theta", "alpha","beta")
tau.sd = data.frame(tau.sd)
rownames(tau.sd) = c("Plan Y", "Plan NF", "Plan F", "Plan N", "Plan C")

cat("These are the parameter estimates:\n")
round(mle, 4)

cat("The asymptotic standard deviation for each parameter for each plan:\n")
round(asym.sd, 4)

cat("The asymptotic standard deviation for tau for each plan:\n")
round(tau.sd, 4)

# save(mle, loglik, time, file="simRun1Fit.Rdata")







