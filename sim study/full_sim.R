#!/usr/bin/env Rscript
##slurm sim 3
source("count3_functions.R")



### SIM STUDY
lamvec = c(1,3)
thetavec = c(0.01, 0.1)
alphavec = c(9,20)
betavec = c(2,3)
npartsvec = c(50,100,200)
nrepsvec = c(4,5)
parvec = expand.grid(lam = lamvec, theta = thetavec, p=pvec, s = npartsvec, r = nrepsvec)
nsim = 1000

mles = array(0,dim=c(4,4,nsim))
logliks = array(0,dim=c(4,2,nsim))

system.time({
  for(j in 1:nrow(parvec)){
    par0=c(parvec$lam[j], parvec$p[j] ,parvec$theta[j])
    nreps = parvec$r[j]
    nparts = parvec$s[j]
    SP.temp =  matrix(NA, nrow = nsim, ncol=3)
    for(i in 1:nsim){                                                        
      cdat.full = fixPdata(nparts,nreps, par0)

      cdat = cdat.full$cdat
      xdat = cdat.full$xdat
      ndat = cdat.full$ndat
      fdat = cdat.full$fdat
  
      
      ### STARTING VALUES
      #lambda, p ,theta
      em.temp = EM(ipar=c(mean(cdat), 0.9,0.1), cdata=cdat, maxiter=1000, atol=1e-4, maxN = 30)
      EM.start = em.temp$mle
      SP.temp[i,] = EM.start
      st.theta = EM.start[3]
      st.lam = EM.start[1]
      ms = EM.start[2]
      as = (ms*(1-ms) - 0.01)/(0.01*(1 - (ms-1)/ms))
      bs = max(as*(1-ms)/ms, 1.05)
      mndat = mean(cdat)
      mfdat = mean(fdat)
      
      ##### RUN the EM

      
      #NF
      temp = EMX.ab(par=c(mndat,mfdat,as,bs), xdat=xdat , fdat=fdat, ndat=ndat, maxiter=1000, atol=1e-2)
      mles[1,,i]= temp$mle
      logliks[1,1,i] = temp$loglik
      
      #N
      temp = EMN(par=c(mndat,st.theta,as,bs), cdat=cdat, ndat = ndat,maxiter=1000, atol=1e-2)
      mles[2,,i]= temp$mle
      logliks[2,1,i] = temp$loglik
      
      ##R 
      temp = EMR(par=c(st.lam, st.theta, as,bs), xdat=xdat, fdat=fdat, M=4, maxiter=1000, atol=1e-2)
      mles[3,,i]= temp$mle
      logliks[3,1,i] = temp$loglik
      
      #SP/C
      temp = EMC(par=c(st.lam, st.theta, as,bs), cdat=cdat, atol=1e-2, maxiter = 1000)
      mles[4,,i]= temp$mle
      logliks[4,1,i] = temp$loglik
      
      ######the fixed case
      n.lik = create.loglikn2(cdat,ndat)
      logliks[2,2,i] = -optim(c(3,0.9,0.1), n.lik)$value
      
      logliks[3,2,i]= EM.r(ipar=c(st.lam, 0.9,0.1), xdat=xdat, fdat=fdat, maxiter=1000, atol=1e-2)$loglik
      
      #GS
      gs.lik = create.loglikGS(xdat, fdat,ndat)
      mle.par = c(mean(ndat), mean(xdat)/mean(ndat), mean(fdat))
      logliks[1,2,i] = gs.lik(mle.par)
      
      logliks[4,2,i] = max(em.temp$loglik)
      
      print(i)
    }
    #saveRDS(mles, file.path("~/Downloads/sim_results_fixedp/",paste("mlep_", j, "_", nparts,"-",nreps, ".rds", sep='') ))
    #saveRDS(logliks, file.path("~/Downloads/sim_results_fixedp/",paste("loglikp_", j, "_", nparts,"-",nreps, ".rds", sep='') ))
  }
  
})





