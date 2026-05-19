
source("~/Downloads/simRun1/simRun1_functions.R")

getZN.bubble <- function(lcomb=NULL, M=NULL, nReps=NULL) {
  #znSet <- readRDS("znSet5.rds")
  #znSet <- readRDS("znSetx12r3M10.rds")
  #lCombFile = rankComb( rep(15,nReps))
  
  znSet <- readRDS("znBubble.rds")
    
    znSet = znSet[1:lcomb]
    znSet = lapply(znSet, function(z=NULL, M=NULL ) { 
      nSet = z$n <= min(z$n) + M
      z$znum = z$znum[nSet]
      z$n    = z$n[nSet]
      z
    }, M=M )
  znSet
}




uniqueCdat <- function(cdat=NULL, wt=NULL) {
  nReps = ncol(cdat)
  
  # sort and find the unique pairs of (ci)
  df = data.frame( t(apply(cdat,1,sort, decreasing=TRUE)) )
  if (is.null(wt)) uCount = as.matrix(aggregate( data.frame(num=rep(1,nrow(cdat))), df, sum))
  else uCount = as.matrix(aggregate( data.frame(num=wt), df, sum))
  
  # cdat = the cdat  
  # uNum = number of times each values occurs
  # n    =  the value of n for 
  list( cdat   = uCount[,1:nReps], cNum = uCount[,nReps+1])
}


create.Cloglik.bubble <- function(cdat=NULL, M=NULL, wt=NULL) {
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  nReps = ncol(cdat); 
  
  ## sort to find unique c rows
  uCount  = uniqueCdat(cdat=cdat, wt=wt) 
  
  ## The max possible C 
  cMax   = max(cdat); nReps = ncol(cdat); 
  rCombs = apply(uCount$cdat, 1, rankComb) 
  #rfx   = getRFX.bubble( rCombs, nReps )

  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs);
  xdat  = getOrdComb(cMax=cMax, nReps=nReps)[1:lcomb,]
  znSet = getZN.bubble(lcomb, M, nReps);
  fdat  = xdat
  
  zseq  = 0:nReps; nSeq = 0:(cMax+ M);
  ZNden = matrix(0, nrow=length(znSet), ncol= (cMax+1) + M )
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    F.den = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    for (i in 1:length(znSet) ) {
      zn = znSet[[i]]
      if (!is.null(zn)) {
        ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
      }  
    }
  
    logC.den <- log(sapply(rCombs, function(i=NULL) {
      fxnum = readRDS( paste("rfxr5Bubble/rfx5c", i, ".rds", sep='') )
      sum( colSums(ZNden[fxnum[,1],,drop=FALSE]*F.den[fxnum[,2]]*fxnum[,3])*dpois( nSeq, par1[1])  )
    }))
    
    return( sum(uCount$cNum*logC.den) )  
  }
}













obsInfC <- function(mle=NULL, cdata=NULL, M=4, wt=NULL) {

  
  nReps = ncol(cdata)
  getNFAB = create.getNFAB(cdat, M)
  
  NFAB  = getNFAB(mle) ## (N, F, a, b)
  jNFAB = jacobian(getNFAB, x=mle)
  sab   = sum(mle[3:4])
  
  H = matrix(0,4,4)
  H[1,1] =       NFAB[1]/mle[1]^2 
  H[1,]  = H[1,] - jNFAB[1,]/mle[1]
  
  H[2,2] = NFAB[2]/mle[2]^2 
  H[2,] = H[2,] - jNFAB[2,]/mle[2]
  
  H[3,3] = NFAB[1]*(trigamma( mle[3] ) - trigamma( sab )) 
  H[3,]  = H[3,] - ( jNFAB[3,]  - jNFAB[1,]*(digamma(mle[3]) - digamma( sab ) ) )
  
  H[4,4] = NFAB[1]*(trigamma( mle[4] ) - trigamma( sab )) 
  H[4,] = H[4,] - ( jNFAB[4,]  - jNFAB[1,]*(digamma(mle[4]) - digamma( sab ) ) )
  
  H[3,4] = H[3,4] - NFAB[1]*trigamma( sab)  
  H[4,3] = H[4,3] - NFAB[1]*trigamma( sab)  
  
  H
}









create.CDen.manyRX <- function(cdat=NULL, M=NULL) {
  ## large number of RX 
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  nReps = ncol(cdat); 
  
  ## The max possible C 
  cMax   = max(cdat); nReps = ncol(cdat); 
  rCombs = apply(cdat, 1, rankComb) 
  #rfx   = getRFX( rCombs, nReps )

  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs);
  xdat  = getOrdComb(cMax=cMax, nReps=nReps)[1:lcomb,]
  znSet = getZN(lcomb, M, nReps);
  fdat  = xdat
  
  zseq  = 0:nReps; nSeq = 0:(cMax+ M);
  ZNden = matrix(0, nrow=length(znSet), ncol= (cMax+1) + M )
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    F.den = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    for (i in 1:length(znSet) ) {
      zn = znSet[[i]]
      if (!is.null(zn)) {
        ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
      }  
    }
  
    ## memory overload if we try to load all at once
    ## instead read and use each rx file one at a time 
    CDen <- sapply(rCombs, function(i=NULL) {
      fxnum = readRDS( paste("rfxr5Bubble/rfx5c", i, ".rds", sep='') )
      sum( colSums(ZNden[fxnum[,1],,drop=FALSE]*F.den[fxnum[,2]]*fxnum[,3])*dpois( nSeq, par1[1])  )
    })
    
    return( CDen )  
  }
}


create.Cloglik.manyRX <- function(cdat=NULL, M=NULL, wt=NULL) {
  ## large number of RX
  
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  nReps = ncol(cdat); 
  
  ## sort to find unique c rows
  uCount  = uniqueCdat(cdat=cdat, wt=wt) 
  
  ## The max possible C 
  cMax   = max(cdat); nReps = ncol(cdat); 
  rCombs = apply(uCount$cdat, 1, rankComb) 
  #rfx   = getRFX( rCombs, nReps )
  
  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs);
  xdat  = getOrdComb(cMax=cMax, nReps=nReps)[1:lcomb,]
  znSet = getZN.bubble(lcomb, M, nReps);
  fdat  = xdat
  
  zseq  = 0:nReps; nSeq = 0:(cMax+ M);
  ZNden = matrix(0, nrow=length(znSet), ncol= (cMax+1) + M )
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    F.den = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    for (i in 1:length(znSet) ) {
      zn = znSet[[i]]
      if (!is.null(zn)) {
        ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
      }  
    }
    
    logC.den <- log(sapply(rCombs, function(i=NULL) {
      fxnum = readRDS( paste("rfxr5Bubble/rfx5c", i, ".rds", sep='') )
      sum( colSums(ZNden[fxnum[,1],,drop=FALSE]*F.den[fxnum[,2]]*fxnum[,3])*dpois( nSeq, par1[1])  )
    }))
    
    return( sum(uCount$cNum*logC.den) )  
  }
}

  
  
fisherInfC.manyRX <- function(par0=NULL, nreps=NULL,  maxC=12, tol=1e-3, M=4) {
  ## there is a large number of RX
  ##### output is expected fisher information 
  cdat = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE)   
  num = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE, keepResults=TRUE,
                     FUN=function(z) { as.numeric(factorial(length(z))/prod(factorial(table( z )))) }, FUN.VALUE=as.numeric(1) )   
  
  den.count = create.CDen.manyRX(cdat, M=M)

  den.df = den.count(par0)*num
  ord.wt = order(den.df, decreasing = TRUE)
  
  df2 = cdat[ord.wt,]
  wt2 = den.df[ord.wt]
  
  ord.wt2 = (cumsum(wt2)/sum(wt2) ) < (1-tol) ### find the leading terms 
  
  cdata = as.matrix(df2[ord.wt2,])
  den       = wt2[ord.wt2]
print(length(den))
  
  fn1 = create.Cloglik.manyRX(cdata, M, den)
  F0 = -1*hessian(fn1, par0)
  F0  
}  









obsInfC.manyRX <- function(mle=NULL, cdata=NULL, M=4, wt=NULL) {
  ##### output is obsInf/n
  
  nReps = ncol(cdata)
  getNFAB = create.getNFAB.manyRX(cdata, M, wt)
  
  NFAB  = getNFAB(mle) ## (N, F, a, b)
  jNFAB = jacobian(getNFAB, x=mle)
  sab   = sum(mle[3:4])
  
  H = matrix(0,4,4)
  H[1,1] =       NFAB[1]/mle[1]^2 
  H[1,]  = H[1,] - jNFAB[1,]/mle[1]
  
  H[2,2] = NFAB[2]/mle[2]^2 
  H[2,] = H[2,] - jNFAB[2,]/mle[2]
  
  H[3,3] = NFAB[1]*(trigamma( mle[3] ) - trigamma( sab )) 
  H[3,]  = H[3,] - ( jNFAB[3,]  - jNFAB[1,]*(digamma(mle[3]) - digamma( sab ) ) )
  
  H[4,4] = NFAB[1]*(trigamma( mle[4] ) - trigamma( sab )) 
  H[4,] = H[4,] - ( jNFAB[4,]  - jNFAB[1,]*(digamma(mle[4]) - digamma( sab ) ) )
  
  H[3,4] = H[3,4] - NFAB[1]*trigamma( sab)  
  H[4,3] = H[4,3] - NFAB[1]*trigamma( sab)  
  
  H
}






create.CDen.manyRX <- function(cdat=NULL, M=NULL) {
  ## large number of RX 
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  nReps = ncol(cdat); 
  
  ## The max possible C 
  cMax   = max(cdat); nReps = ncol(cdat); 
  rCombs = apply(cdat, 1, rankComb) 
  #rfx   = getRFX( rCombs, nReps )

  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs);
  xdat  = getOrdComb(cMax=cMax, nReps=nReps)[1:lcomb,]
  znSet = getZN.bubble(lcomb, M, nReps);
  fdat  = xdat
  
  zseq  = 0:nReps; nSeq = 0:(cMax+ M);
  ZNden = matrix(0, nrow=length(znSet), ncol= (cMax+1) + M )
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    F.den = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    for (i in 1:length(znSet) ) {
      zn = znSet[[i]]
      if (!is.null(zn)) {
        ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
      }  
    }
  
    ## memory overload if we try to load all at once
    ## instead read and use each rx file one at a time 
    CDen <- sapply(rCombs, function(i=NULL) {
      fxnum = readRDS( paste("rfxr5Bubble/rfx5c", i, ".rds", sep='') )
      sum( colSums(ZNden[fxnum[,1],,drop=FALSE]*F.den[fxnum[,2]]*fxnum[,3])*dpois( nSeq, par1[1])  )
    })
    
    return( CDen )  
  }
}


  
fisherInfC.manyRX <- function(par0=NULL, nreps=NULL,  maxC=12, tol=1e-3, M=4, numDerv=FALSE ) {
  ## there is a large number of RX
  ##### output is expected fisher information 
  cdat = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE)   
  num = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE, keepResults=TRUE,
                     FUN=function(z) { as.numeric(factorial(length(z))/prod(factorial(table( z )))) }, FUN.VALUE=as.numeric(1) )   
  
  den.count = create.CDen.manyRX(cdat, M=M)

  den.df = den.count(par0)*num
  ord.wt = order(den.df, decreasing = TRUE)
  
  df2 = cdat[ord.wt,]
  wt2 = den.df[ord.wt]
  
  ord.wt2 = (cumsum(wt2)/sum(wt2) ) < (1-tol) ### find the leading terms 
  
  cdata = as.matrix(df2[ord.wt2,])
  den       = wt2[ord.wt2]
print(length(den))
  
  if (numDerv) {
    fn1 = create.Cloglik.bubble(cdata, M, den)
    F0 = -1*hessian(fn1, par0)
  } else {
    F0 = obsInfC.manyRX(mle=par0, cdata=cdata, M=M, wt=den)
  }  
  F0 
}  


create.getNFAB.manyRX <- function(cdat=NULL, M=4,  wt=NULL) {
  nReps = ncol(cdat); 

  uData  = uniqueCdat(cdat=cdat, wt=wt) 
  rCombs = apply(uData$cdat, 1, rankComb) 
 # rfx    = getRFX( rCombs, nReps ) ### list of matrices (fi, xi, num)
  
  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs)
  xdat  = getOrdComb(cMax=max(uData$cdat), nReps=nReps)[1:lcomb,]
  znSet = getZN.bubble(lcomb, M, nReps)
  fdat  = xdat
  
  zseq     = 0:nReps;
  maxN1M   = max(uData$cdat)+1 + M
  nSeq     = 0:(maxN1M -1)
  sum.fdat = rowSums(fdat)
  
  function(par1=NULL) {
    Fden    = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    znTable = getZNTable2(par1[c(3,4)], znSet, nReps, maxN1M) 
    
    cfs12 = sapply(rCombs, function(i=NULL) {
      fxnum = readRDS( paste("rfxr5Bubble/rfx5c", i, ".rds", sep='') )
      
      if (nrow(fxnum) == 1) {
        xs12  = znTable[fxnum[,1], ,,drop=FALSE ]
        xfn.den =  t( matrix(xs12[,,1], ncol=length(nSeq), nrow=1, byrow=TRUE) *Fden[fxnum[,2]]*fxnum[,3])*dpois(nSeq, par1[1]) 
      } else {
        xs12  = znTable[fxnum[,1],,]
        xfn.den =  t( xs12[,,1]*Fden[fxnum[,2]]*fxnum[,3])*dpois(nSeq, par1[1]) 
      }
      cden = sum(xfn.den)
      xfn.cden = xfn.den/cden
      
      xf.cden = colSums(xfn.cden)
      sum.f   = sum( sum.fdat[fxnum[,2]] * xf.cden ) 
      
      n.cden = rowSums(xfn.cden)
      exp.n  = sum( nSeq*n.cden )
      
      s1  = sum( xs12[,,2]*t(xfn.cden) )
      s2  = sum( xs12[,,3]*t(xfn.cden) )
      
      c(exp.n, sum.f, s1, s2)
    } ) 
    sum.fs12 = colSums( t(cfs12)*uData$cNum)
    
    S12 = sum.fs12[c(3,4)] - sum.fs12[1]*digamma( sum(par1[c(3,4)]) + nReps )
    c(sum.fs12[c(1,2)], S12 )
  }
}

