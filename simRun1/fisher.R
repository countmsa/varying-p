######################## Fisher Informations ###############
library(pracma)
### Plan Y

dBBy <- function(y=NULL, n=NULL, alpha=NULL, beta=NULL) {
  exp( lchoose(n, y) + lbeta(y + alpha, n - y + beta) - lbeta(alpha, beta) )
}

fisherInfY.ab <- function(par=NULL, nreps=NULL) {
  ab = par[3:4]
  zseq = 0:nreps
  zden = dBBy(zseq, nreps, ab[1], ab[2])
  H = matrix(0,4,4)
  
  H[1,1] = 1/(par[1])
  H[2,2] = nreps/par[2]
  H[3,3] = par[1]*trigamma( ab[1])  - par[1]*sum( zden*trigamma(        zseq + ab[1] )) 
  H[4,4] = par[1]*trigamma( ab[2])  - par[1]*sum( zden*trigamma( nreps - zseq + ab[2])) 
  
  # all elements 
  H[3:4,3:4] = H[3:4,3:4] + par[1]*trigamma(nreps+ sum(ab)) - par[1]*trigamma(sum(ab)) 
  H
}


### Plan NF

XNDen2 <- function(lab=NULL, maxX=NULL, nReps=NULL, M=NULL) {
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  
  xdat = comboGeneral(v=0:maxX, m=nReps, repetition =TRUE)   
  num  = comboGeneral(v=0:maxX, m=nReps, repetition =TRUE, keepResults=TRUE,
                      FUN=function(z) { as.numeric(factorial(length(z))/prod(factorial(table( z )))) }, FUN.VALUE=as.numeric(1) )   
  
  rx     = apply(xdat,1, rankComb )
  num   = num[order(rx)]
  
  znSet  = getZN( max(rx), M, nReps )
  
  zseq   = 0:nReps;
  ZNden  = matrix(0, nrow=length(znSet), ncol= maxX  +M +1 )
  nseq   = 0:(maxX  +M)
  
  # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
  logdBB = logdBBm(zseq, nReps, lab[2], lab[3] )
  for (i in 1:length(znSet) ) {
    zn = znSet[[i]]
    ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
  }
  
  XNden = ZNden*outer(num, dpois(nseq, lab[1]) )
  return(XNden)
}


create.getAB2 <- function(maxX=NULL, nReps=NULL, M=NULL, wt=NULL, lambda=NULL) {
  
  if (is.null(wt)) stop("wt is null")
  
  maxRx  = rankComb( rep(maxX, nReps) )
  znSet  = getZN( maxRx, M, nReps )
  
  maxN1M   = maxX+M+1
  
  function(ab=NULL) {
    s12.xn = getZNTable2(ab, znSet, nReps, maxN1M) 
    #s12.xn   = array(0, dim=c(length(znSet), maxN1M, 3) ) 
    
    s1.xn = sum(s12.xn[,,2]*wt)
    s2.xn = sum(s12.xn[,,3]*wt)
    
    S12 = c(s1.xn, s2.xn) - lambda*digamma( sum(ab) + nReps )
    S12
  }
}


fisherInfX.lab2 <- function(lab=NULL, maxX=12, nReps=4, M=5) {
  ##### mle is for ab
  ##### output is obsInf/n
  
  ab = lab[-1]
  
  wt = XNDen2(lab, maxX, nReps, M)
  
  getAB = create.getAB2(maxX, nReps, M, wt, lab[1] )
  AB    = getAB(ab) ## (N, F, a, b)
  
  jAB   = jacobian(getAB, x=ab)
  sab   = sum(ab)
  
  H = matrix(0, 2, 2)
  H[1,1] = lab[1]*(trigamma( ab[1] ) - trigamma( sab )) 
  H[1,]  = H[1, ] - ( jAB[1,] )
  
  H[2,2] = lab[1]*(trigamma( ab[2] ) - trigamma( sab )) 
  H[2,]  = H[2,] - ( jAB[2,]  )
  
  H[1,2] = H[2,1] - lab[1]*trigamma( sab )  
  H[2,1] = H[2,1] - lab[1]*trigamma( sab )  
  
  H2 = matrix(0,3,3)
  H2[1,1]      = 1/lab[1]
  H2[2:3,2:3] = H
  
  H2
}




fisherInfX2 <- function(par=NULL, maxX=12, nReps=4, M=5) {
  
  Flab = fisherInfX.lab2(lab=par[c(1,3,4)], nReps=nReps, maxX=maxX,  M=M )
  H = matrix(0,4,4)
  H[-2, -2] = Flab
  H[2,2] = nReps/par[2]
  H
}



### Plan N

obsInfN.tab <- function(tab=NULL, cdat=NULL, ndat=NULL, wt=NULL) {
  if (is.null(wt)) wt = rep(1, length(ndat))
  getFAB = create.getFAB(cdat, ndat, wt)
  
  FAB  = getFAB(tab) ## (F, a, b)
  jFAB = jacobian(getFAB, tab)
  sab  = sum(tab[2:3])
  
  sum.ndat.wt = sum(wt*ndat)
  
  H = matrix(0,3,3)
  
  H[1,1] = FAB[1]/tab[1]^2 
  H[1,]  = H[1,] - jFAB[1,]/tab[1]
  
  H[2,2] = sum.ndat.wt*(trigamma( tab[2] ) - trigamma( sab )) 
  H[2,]  = H[2, ] - ( jFAB[2,] )
  
  H[3,3] = sum.ndat.wt*(trigamma( tab[3] ) - trigamma( sab )) 
  H[3,]  = H[3,] - ( jFAB[3,]  )
  
  H[2,3] = H[3,2] - sum.ndat.wt*trigamma( sab )  
  H[3,2] = H[3,2] - sum.ndat.wt*trigamma( sab )  
  
  H
}


obsInfN <- function(mle=NULL, cdat=NULL, ndat=NULL) {
  H3 = obsInfN.tab(tab=mle[c(2,3,4)], cdat, ndat)
  H        = matrix(0,4,4)
  H[1,1]   = sum(ndat)/mle[1]^2 
  H[-1,-1] = H3
  H
}


create.NDen <- function(cdat=NULL, ndat=NULL) {
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  nReps = ncol(cdat); 
  samSize  = length(ndat);
  sum.ndat = sum(ndat); 
  
  M      = max(ndat - apply(cdat,1,max))
  rCombs = apply(cdat, 1, rankComb) 
  rfx    = getRFX( rCombs, nReps )
  
  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs)
  xdat  = getOrdComb(cMax=max(cdat), nReps=nReps)[1:lcomb,]
  znSet = getZN(lcomb, M, nReps )
  fdat  = xdat
  
  sum.ndat = sum(ndat); 
  zseq = 0:nReps;
  ZNden    = matrix(0, nrow=length(znSet), ncol= (max(cdat)+1) + M )
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    logN.den = dpois(ndat, par1[1], log=TRUE)
    F.den    = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    for (i in 1:length(znSet) ) {
      zn = znSet[[i]]
      ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
    }
    
    logC.den = log(vapply(1:length(rfx), function(m=NULL) {
      fxnum = rfx[[m]]
      sum(ZNden[fxnum[,1], ndat[m]+1]*F.den[fxnum[,2]]*fxnum[,3])
    }, numeric(1) ))   
    
    return( exp(logN.den + logC.den) )  
  }
}


fisherInfN <- function(par=NULL, nreps=NULL, maxC=12, tol=1e-4, M=4 ) {
  
  cdat.tem = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE)   
  num.tem  = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE, keepResults=TRUE,
                          FUN=function(z) { as.numeric(factorial(length(z))/prod(factorial(table( z )))) }, FUN.VALUE=as.numeric(1) )   
  nseq = 0:(maxC + M)
  
  cdat = do.call(rbind, lapply(nseq, function(z) { cdat.tem }) )
  ndat = rep(nseq, each=nrow(cdat.tem))
  num  = rep(num.tem, length(nseq) )
  
  ## Subset by M
  Mx  = ndat - apply(cdat,1,max)
  
  cdat = cdat[Mx <= M, ]
  ndat = ndat[Mx <= M ]
  num  = num[Mx <= M ]
  
  den.count = create.NDen(cdat, ndat)
  
  den.df = den.count(par)*num
  ord.wt = order(den.df, decreasing = TRUE)
  
  df2 = cdat[ord.wt,]
  nf2 = ndat[ord.wt]
  wt2 = den.df[ord.wt]
  
  ord.wt2 = cumsum(wt2) < (1-tol) ### find the leading terms 
  
  cdat = as.matrix(df2[ord.wt2,])
  ndat = nf2[ord.wt2]
  den   = wt2[ord.wt2]
  
  H3       = obsInfN.tab(tab=par[2:4], cdat=cdat, ndat=ndat, wt=den)
  H        = matrix(0,4,4)
  H[1,1]   = 1/par[1] 
  H[-1,-1] = H3
  H
}


### Plan F

create.XDen <- function(xdat=NULL, M=5) {
  # output is f(lab)
  # xdat is a n x r matrix 
  # ni  goes from  max(xi1,...,xir) to max(xi1,...,xir) + M
  
  # Some useful variables
  zseq = 0:ncol(xdat) ; 
  nReps = ncol(xdat)
  
  #xdat -> ydat -> zn.dat 
  rXCombs = apply(xdat, 1, rankComb) 
  zn.dat  = getZNSet(combSet=rXCombs, M=M, nReps = ncol(xdat))
  
  function(lab=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    
    # ( (maxN +1)  x nrow(x)   ) matrix for f( x | n) 
    logdBB = logdBBm(zseq, nReps, lab[2], lab[3] )
    Xden = vapply(zn.dat, function(zn=NULL) {
      xden.n = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1)) 
      sum( xden.n*dpois(zn$n, lab[1]) )
    }, numeric(1) )   
    
    return( Xden )  
  }
}


fisherInfR <- function(par=NULL, nreps=NULL,  maxX=12, tol=1e-4, M=4, numDerv=FALSE ) {
  ## par = (lambda, theta, alpha, beta)
  ## only needs to constuct xdat
  
  xdat = comboGeneral(v=0:maxX, m=nreps, repetition =TRUE)   
  num = comboGeneral(v=0:maxX, m=nreps, repetition =TRUE, keepResults=TRUE,
                     FUN=function(z) { as.numeric(factorial(length(z))/prod(factorial(table( z )))) }, FUN.VALUE=as.numeric(1) )   
  
  den.count = create.XDen(xdat, M=M)
  
  den.df = den.count(par[c(1,3,4)])*num
  ord.wt = order(den.df, decreasing = TRUE)
  
  df2 = xdat[ord.wt,]
  wt2 = den.df[ord.wt]
  
  ord.wt2 = cumsum(wt2) < (1-tol) ### find the leading terms 
  
  xdat = as.matrix(df2[ord.wt2,])
  den       = wt2[ord.wt2]
  
  if (numDerv) {
    fdat = matrix(par[2], 1, nreps)  
    fn1  = create.Rloglik(xdat, fdat, M, den)
    F0   = -1*hessian(fn1, par0)
  } else {
    F3 = obsInfR.lab(lab=par[c(1,3,4)], xdat, M, wt=den)
    F0 = matrix(0,4,4)
    F0[2,2]   = nreps/par[2]
    F0[-2,-2] = F3
  }    
  F0  
} 


obsInfR.lab <- function(lab=NULL, xdat=NULL, M=4, wt=NULL) {
  ## The fisher for lab = (lambda, alpha, beta)
  
  getNAB = create.getNAB(xdat, M, wt)
  
  NAB  = getNAB(lab) ## (N, a, b)
  jNAB = jacobian(getNAB, lab)
  sab  = sum(lab[2:3])
  
  H = matrix(0,3,3)
  H[1,1] =  NAB[1]/lab[1]^2 
  H[1,] = H[1,]  - jNAB[1,]/lab[1]
  
  H[2,2]  = NAB[1]*(trigamma( lab[2] ) - trigamma( sab )) 
  H[2,] = H[2,] - ( jNAB[2,]  - jNAB[1,]*(digamma(lab[2]) - digamma( sab ) ) )
  
  H[3,3]  = NAB[1]*(trigamma( lab[3] ) - trigamma( sab )) 
  H[3,] = H[3,] - ( jNAB[3,]  - jNAB[1,]*(digamma(lab[3]) - digamma( sab ) ) )
  
  H[2,3] = H[2,3] - NAB[1]*trigamma( sab)  
  H[3,2] = H[3,2] - NAB[1]*trigamma( sab)  
  
  H
}

obsInfR <- function(mle=NULL, xdat=NULL, fdat=NULL, M=4, wt=NULL) {
  H3 = obsInfR.lab(lab=mle[c(1,3,4)], xdat, M, wt)
  
  H = matrix(0,4,4)
  H[2,2]   = sum(fdat)/mle[2]^2
  H[-2,-2] = H3
  
  H
}



### Plan C
obsInfC <- function(mle=NULL, cdata=NULL, M=4, wt=NULL) {
  ##### output is obsInf/n
  
  nReps = ncol(cdata)
  getNFAB = create.getNFAB(cdata, M,wt=wt)
  
  NFAB  = getNFAB(mle) ## (N, F, a, b)
  jNFAB = jacobian(getNFAB, x=mle)
  sab   = sum(mle[3:4])
  
  H = matrix(0,4,4)
  H[1,1] =       NFAB[1]/(mle[1]^2 )
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

create.CDen <- function(cdat=NULL, M=10) {
  # cdat is a n x r matrix  
  ## The max possible C 
  cMax   = max(cdat); nReps = ncol(cdat); 
  rCombs = apply(cdat, 1, rankComb) 
  rfx   = getRFX( rCombs, nReps )
  
  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs);
  xdat  = getOrdComb(cMax=cMax, nReps=nReps)[1:lcomb,]
  znSet = getZN(lcomb, M, nReps);
  fdat  = xdat
  
  zseq  = 0:nReps; nSeq = 0:(cMax+ M);
  nReps1 = nReps+1; nReps2  = nReps+2; nReps3  = nReps+3; ## n row
  ZNden = matrix(0, nrow=length(znSet), ncol= (cMax+1) + M )
  nReps23 = c(nReps2,nReps3)
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    F.den = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    # znTable = getZNTable2(par1[c(3,4)], znSet, nReps, cMax + M) 
    
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    for (i in 1:length(znSet) ) {
      zn = znSet[[i]]
      temp = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
      #denZN = exp( colSums(zn[-nReps23,,drop=FALSE]*logdBB))*zn[nReps2,]
      #denZN = aggregate(list(den=denZN), by= data.frame(nSet = zn[nReps3,]), sum)
      denZN = aggregate(list(den=temp), by= data.frame(nSet = zn$n), sum)
      ZNden[i, denZN$nSet+1] = denZN$den
    }
    
    Cden = sapply(rfx, function(fxnum=NULL) {
      sum( colSums(ZNden[fxnum[,1],,drop=FALSE]*F.den[fxnum[,2]]*fxnum[,3])*dpois( nSeq, par1[1])  )
    } )   
    return( Cden )  
  }
}


fisherInfC <- function(par0=NULL, nreps=NULL,  maxC=12, tol=1e-4, M=4, numDerv=FALSE ) {
  ##### output is expected fisher information 
  cdat = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE)   
  num = comboGeneral(v=0:maxC, m=nreps, repetition =TRUE, keepResults=TRUE,
                     FUN=function(z) { as.numeric(factorial(length(z))/prod(factorial(table( z )))) }, FUN.VALUE=as.numeric(1) )   
  
  den.count = create.CDen(cdat, M=M)
  
  den.df = den.count(par0)*num
  ord.wt = order(den.df, decreasing = TRUE)
  
  df2 = cdat[ord.wt,]
  wt2 = den.df[ord.wt]
  
  ord.wt2 = cumsum(wt2) < (1-tol) ### find the leading terms 
  
  cdata = as.matrix(df2[ord.wt2,])
  den       = wt2[ord.wt2]
  
  
  if (numDerv) {
    fn1 = create.Cloglik(cdata, M, den)
    F0 = -1*hessian(fn1, par0)
  } else {
    F0 = obsInfC(mle=par0, cdata=cdata, M=M, wt=den)
  }  
  F0  
}  

### for function of the parameters

tau_derivs=function(par){
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

mu_derivs=function(ab){
  a = ab[1]
  b = ab[2]
  a_der = 1 / (a + b) - a / (a + b) ^ 2
  b_der =  -a / (a + b) ^ 2
  return(c(a_der, b_der))
}

var.bb_derivs = function(ab){
  a = ab[1]
  b = ab[2]
  a_der = (a * b / (a + b) ^ 2 / (a + b + 1)) ^ (-1 / 0.2e1) * 
    (b / (a + b) ^ 2 / (a + b + 1) - 2 * a * b / (a + b) ^ 3 / (a + b + 1) -
       a * b / (a + b) ^ 2 / (a + b + 1) ^ 2) / 2
  b_der = (a * b / (a + b) ^ 2 / (a + b + 1)) ^ (-0.1e1 / 0.2e1) * 
    (a / (a + b) ^ 2 / (a + b + 1) - 2 * a * b / (a + b) ^ 3 / (a + b + 1) - 
       a * b / (a + b) ^ 2 / (a + b + 1) ^ 2) / 2
  return(c(a_der, b_der))
}

var.tau <- function(fisherInfpar=NULL, par=NULL) {
  dpar = tau_derivs(par)
  
  as.numeric( t(dpar) %*% solve(fisherInfpar) %*% (dpar) )
}
var.mu <- function(fisherInfpar=NULL, par=NULL) {
  dpar = mu_derivs(par)
  
  as.numeric( t(dpar) %*% solve(fisherInfpar) %*% (dpar) )
}
var.sig <- function(fisherInfpar=NULL, par=NULL) {
  dpar = var.bb_derivs(par)
  
  as.numeric( t(dpar) %*% solve(fisherInfpar) %*% (dpar) )
}


