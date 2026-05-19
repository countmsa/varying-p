##
## These are the functions that used by more than one plan 
##


getZN <- function(lcomb=NULL, M=NULL, nReps=NULL) {
  #znSet <- readRDS("znSet5.rds")
  #znSet <- readRDS("znSetx12r3M10.rds")
  lCombFile = rankComb( rep(25,nReps))
  
  if ( lcomb <= lCombFile & (nReps == 5) & M <= 6) {
    znSet <- readRDS("znBubble.rds")

    znSet = znSet[1:lcomb]                     
    znSet = lapply(znSet, function(z=NULL, M=NULL ) { 
      if (!is.null(z)) {
        nSet = z$n <= min(z$n) + M
        z$znum = z$znum[nSet]
        z$n    = z$n[nSet]
      }  
      z
    }, M=M )
  } else {
    ## OR maybe return NULL instead? easier for simulations 
    ## just use try? stop()
    if ( lcomb > lCombFile ) warning("There is a combination is larger than (",paste0( rep(12, nReps-1), collpase=',' ), 12, ")" )    
    if ( all(nReps != c(3,4)) ) warning("nReps is not one of {3,4}")
    if ( M > 10 ) warning("M has to be <= 10")
    znSet = NULL
  }  
  
  znSet
}


getRFX <- function(rcombs=NULL, nReps=NULL) {
  ## nReps can be equal to {3,4,5}
  ## The largeret combinations is rep(12,nReps)
  lComb = rankComb( rep(25,nReps))
  
  if ( (max(rcombs) <= lComb) & (nReps == 5) ) {
    rfx <- sapply(rcombs, function(i=NULL) {
      readRDS( paste("rfxr5Bubble/", "rfx", nReps, "c", i, ".rds", sep='') )
    })
  } else {
    ## OR maybe return NULL instead? easier for simulations 
    ## just use try? stop()
    if ( (max(rcombs) > lComb)) warning("There is a combination is larger than (",paste0( rep(12, nReps-1), collpase=',' ), 12, ")" )    
    if ( all(nReps != c(3,4,5)) ) warning("nReps is not one of {3,4,5}")
    rfx = NULL
  }
  rfx
}



getZNTable2 <- function(ab=NULL, znSet=NULL, nReps=NULL, maxN1M=NULL) {
  ZNden2   = array(0, dim=c(length(znSet), maxN1M, 3) ) 
  
  logdBB = logdBBm(0:nReps, nReps, ab[1], ab[2] )
  digam1 = digamma( ab[1] + 0:nReps )
  digam2 = digamma( ab[2] + nReps - 0:nReps )
  
  for (i in 1:length(znSet) ) {
    zn = znSet[[i]]
    # if (i== 1) stop('here')
    if (!is.null(zn)) {
      ZNden2[i, zn$n+1, ] = t(vapply( zn$znum, function(tem=NULL) { 
        den = exp(colSums(tem$z*logdBB))*tem$num
        denzn =  sum(den)
        con.den = den/denzn
        s1 = sum( colSums(tem$z*digam1)*con.den )
        s2 = sum( colSums(tem$z*digam2)*con.den )
        c(denzn, s1, s2)
      }, numeric(3)))
    }  
  }
  ZNden2
}




atiken <- function(loglik = NULL, k = NULL, eps = 1e-2) {
  continue = TRUE
  if (k > 3) {
    lm1 = loglik[(k-1)]
    lm = loglik[(k - 2)]
    lm_1 = loglik[(k - 3)]
    
    am = (lm1 - lm)/(lm - lm_1)
    lm1.Inf = lm + (lm1 - lm)/(1 - am)
    val = lm1.Inf - lm
    
    if (val < eps & val >= 0) 
      continue = FALSE
  }
  return(continue)
}


ab_update <- function(ab=NULL, S12=NULL) {
  # ab = (alpha, beta)
  # S12 = avg. of (log(x), log(1-x))
  
  g = -S12 +  digamma(ab) - digamma( sum(ab) )
  H    = diag(trigamma(ab)) -  trigamma( sum(ab) ) 
  # Hessian of NLL
  #if (ab[1] < 0.2) stop("here")
  # Solve Newton system: H * step = -g
  step <- tryCatch( solve(H, g),
                    error = function(e) {  -solve(H + 1e-8 * diag(2), g) }
  )
  ab - step
}





logdBBm <- function(x=NULL, n=NULL, alpha=NULL, beta=NULL) {
  ## We remove the choose(n,x) because we count the number of
  ## matrices with this arrangement differently.
  lbeta(x + alpha, n - x + beta) - lbeta(alpha, beta) 
}

dBBm <- function(x=NULL, n=NULL, alpha=NULL, beta=NULL) {
  ## We remove the choose(n,x) because we count the number of
  ## matrices with this arragement differently.
  exp( lbeta(x + alpha, n - x + beta) - lbeta(alpha, beta) )
}






rankComb <- function(z=NULL) {
  ## sorting the combination increasing order
  zSort = sort(z)
  zSeq = 0:(length(z)-1)
  ## add one as R does index 0
  sum( choose(zSort + zSeq, 1 + zSeq )) +1
}



getOrdComb <- function(cMax=NULL, nReps=NULL) {
  ## Ordered Combination
  combDat = comboGeneral(0:cMax, nReps, repetition = TRUE)[,nReps:1]
  OrdC = apply(combDat, 1, rankComb)
  combDat[order(OrdC), ]
}  







atiken <- function(loglik = NULL, k = NULL, eps = 1e-2) {
  continue = TRUE
  if (k > 3) {
    lm1 = loglik[(k-1)]
    lm = loglik[(k - 2)]
    lm_1 = loglik[(k - 3)]
    
    am = (lm1 - lm)/(lm - lm_1)
    lm1.Inf = lm + (lm1 - lm)/(1 - am)
    val = lm1.Inf - lm
    
    if (val < eps & val >= 0) 
      continue = FALSE
  }
  return(continue)
}



EM <- function(ipar=NULL, cdata=NULL, maxiter=100, atol=1e-2, maxN= 30) {
  #maxN = qpois(1-1e-4, mean(cdata) ) # for lambda = 30 this is 52
  df = data.frame( t(apply(cdata,1,sort)) )
  df = aggregate(list(num=rep(1,nrow(df))), df, length)
  wt = df$num
  countdata = as.matrix(df[,1:ncol(cdata)])
  
  mle  = ipar
  llik = numeric(maxiter)
  
  k = 1
  while (atiken(llik, k = k, eps = atol) & k < maxiter) {
    tem = EM.step(par=mle, countdata=countdata, wt=wt, maxN=maxN )
    # M-step
    mle = tem[-4]
    llik[k] = tem[4]
    k = k + 1
  }
  val = list( k = k, mle = mle,loglik = llik[(k-1)] )
  #val = list( loglik = llik[1:(k-1)], mle = mle)
  return(val)
}

denCTN3 <- function(kseq=2:3, N=2, p=0.95, theta=2) {
  xseq = 0:N
  den = outer(xseq, kseq, function(x,y) { 
    dbinom(x, size=N, prob=p)*dpois(y-x, lambda=theta) 
  })
  x.y = rep(kseq, each=N+1) - rep(xseq, times=length(kseq))  
  rbind( colSums(den), colSums(xseq*den), colSums( den*x.y) )
}


EM.step <- function(par=NULL, countdata=NULL, wt=NULL, maxN=30) {
  nSeq  = 0:maxN  
  cMin  = min(countdata) 
  cSeq  = cMin:max(countdata) 
  
  # (3 x maxN x maxC) array
  # for f( ci | n) b/c they are same when ci=k
  cxf.n = vapply(nSeq, function(k) { 
    denCTN3(cSeq, N=k, p=par[2], theta=par[3]) 
  }, FUN.VALUE = matrix(0, nrow=3, ncol=length(cSeq)) )
  
  tem  = vapply(1:nrow(countdata), function(k) { 
    z = countdata[k,]
    den.c = cxf.n[1, z + 1 - cMin,]  # ( 1 - cMin) converts the counts to row index 
    d.c = apply(den.c, 2, prod)   
    
    rbind(d.c, 
          rowSums(t(cxf.n[2, z + 1 - cMin,]/den.c)*d.c),
          rowSums(t(cxf.n[3, z + 1 - cMin,]/den.c)*d.c ) )
  }, FUN.VALUE = matrix(0, nrow=3, ncol=length(nSeq) ) )
  
  dn        = dpois(nSeq, lambda=par[1])
  dnCvec    = colSums( tem[1,,]*dn )
  wt.dnCvec = wt/dnCvec
  wt.sum    = sum(wt)
  
  N.c = sum( colSums( nSeq*tem[1,,]*dn )*wt.dnCvec )/sum(wt)
  F.c = sum( colSums( tem[3,,]* dn     )*wt.dnCvec )/(ncol(countdata)*wt.sum )
  X.c = sum( colSums( tem[2,,]* dn     )*wt.dnCvec )/(ncol(countdata)*wt.sum )
  
  val = c(N.c, X.c/N.c, F.c, sum(wt*log(dnCvec)))
  return(val)
}

