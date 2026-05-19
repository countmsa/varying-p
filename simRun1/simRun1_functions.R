library(RcppAlgos)
ab_update <- function(ab=NULL, S12=NULL) {
  # ab = (alpha, beta)
  # S12 = avg. of (log(x), log(1-x))
  alpha <- ab[1]
  beta  <- ab[2]
  
  # current (mu, phi)
  phi <- alpha + beta
  mu  <- alpha / phi
  
  # transformed parameters
  nu  <- qlogis(mu)   # log(mu / (1-mu))
  eta <- log(phi)
  
  g = -S12 +  digamma(ab) - digamma( sum(ab) )
  H    = diag(trigamma(ab)) -  trigamma( sum(ab) ) 
  H_theta <- matrix(0, 2, 2)
  
  A <- phi * mu * (1 - mu)
  
  # Jacobian J = d(alpha,beta)/d(nu,eta)
  J <- matrix(c(
    A, alpha,
    -A, beta
  ), nrow = 2, byrow = TRUE)
  
  # transformed gradient
  g_theta <- c(
    A * (g[1] - g[2]),
    alpha * g[1] + beta * g[2]
  )
  
  # ---- main Hessian term ----
  corr <- matrix(c(
    A * (1 - 2 * mu) * (g[1] - g[2]),
    A * (g[1] - g[2]),
    A * (g[1] - g[2]),
    alpha * g[1] + beta * g[2]
  ), nrow = 2, byrow = TRUE)
  
  H_theta <- t(J) %*% H %*% J + corr
  
  
  # Hessian of NLL
  #if (ab[1] < 0.2) stop("here")
  # Solve Newton system: H * step = -g
  step <- tryCatch( solve(H_theta, g_theta),
                    error = function(e) {  -solve(H_theta + 1e-4 * diag(2), g_theta) }
  )
  
  theta <- c(nu, eta)
  theta_new <- theta - step
  
  
  # map back to (alpha, beta)
  mu_new  <- plogis(theta_new[1])
  phi_new <- exp(theta_new[2])
  
  alpha_new <- mu_new * phi_new
  beta_new  <- (1 - mu_new) * phi_new
  
  c(alpha_new, beta_new)
}

alpha_V = function(par0 = NULL){
  l = par0[1]
  t = par0[2]
  a = par0[3]
  b = par0[4]
  mu = a / (a+b)
  l* mu^2 * (( a + 1 )*(a + b)/ (a*(a + b + 1) )) / (l * mu + t)
}


rPPdata <- function(nparts=20, nreps=3, par=NULL) {
  N = rpois( nparts, par[1])
  pdat = lapply(N, function(n){rbeta(n, par[3], par[4]) })
  xdat = matrix(0,nrow = nparts, ncol = nreps)
  zdat = matrix(0, nrow = nparts, ncol = (nreps + 1))
  for(m in 1:length(N)){
    p = pdat[[m]]
    n = N[m]
    
    if (n == 0) { 
      y = NA
      val = rep(0, nreps)
    } else if (n == 1) { 
      val = rbinom( nreps, size=1, prob=p)
      y = sum(val)
    } else {
      D = sapply(p, function(prob=NULL) { rbinom( nreps, size=1, prob=prob) })
      y = colSums( D )
      val = rowSums( D )
    }
    zdat[m,] =  as.numeric(table(factor(y, 0:nreps)))
    xdat[m,] = val
  }
  
  
  fdat = matrix( rpois( nparts*nreps, par[2]) , nrow=nparts, ncol=nreps)
  
  counts = list(zdat= zdat,xdat=xdat, fdat=fdat, ndat=N, pdat=pdat )
  return(counts)
}



## improved verision
getZN <- function(lcomb=NULL, M=NULL, nReps=NULL) {
  #znSet <- readRDS("znSet5.rds")
  #znSet <- readRDS("znSetx12r3M10.rds")
  lCombFile = rankComb( rep(15,nReps))
  
  if ( lcomb <= lCombFile & any(nReps == 4) & M <= 5 ) {
    #znSet <- readRDS(paste0("znSetx12r", nReps, "M5.rds", collpase=""))
    #if (nReps ==3) {znSet <- readRDS(paste0("znSetx30r", nReps, "M10.rds", collpase=""))
    #} else znSet <- readRDS(paste0("znSetx25r", nReps, "M6.rds", collpase=""))
    znSet = readRDS(paste0("znSetx12r4M5.rds", collpase=""))
    
    znSet = znSet[1:lcomb]
    znSet = lapply(znSet, function(z=NULL, M=NULL ) { 
      nSet = z$n <= min(z$n) + M
      z$znum = z$znum[nSet]
      z$n    = z$n[nSet]
      z
    }, M=M )
  } else {
    ## OR maybe return NULL instead? easier for simulations 
    ## just use try? stop()
    if ( lcomb > lCombFile ) warning("There is a combination is larger than (",paste0( rep(12, nReps-1), collpase=',' ), 12, ")" )    
    if ( all(nReps != c(3,4,5)) ) warning("nReps is not one of {3,4}")
    if ( M > 10 ) warning("M has to be <= 5")
    znSet = NULL
  }  
  
  znSet
}


### NEEDS to be improved 
#getRFXOLD <- function(rcombs=NULL) {
#  ## rcombs is vector of combinations to get
#  rfx <- readRDS("rfx20.rds")
#  rfx[rcombs]
#}
## improved verision
getRFX <- function(rcombs=NULL, nReps=NULL) {
  ## nReps can be equal to {3,4,5}
  ## The largeret combinations is rep(12,nReps)
  lComb = rankComb( rep(15,nReps))
  
  if ( (max(rcombs) <= lComb) & any(nReps == 4) ) {
    if (nReps == 4) {
      rfx <- sapply(rcombs, function(i=NULL) {
        readRDS( paste("rfx15r", nReps, "/", "rfx15r", nReps, "c", i, ".rds", sep='') )
      })
    }
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
    ZNden2[i, zn$n+1, ] = t(vapply( zn$znum, function(tem=NULL) { 
      den = exp(colSums(tem$z*logdBB))*tem$num
      denzn =  sum(den)
      con.den = den/denzn
      s1 = sum( colSums(tem$z*digam1)*con.den )
      s2 = sum( colSums(tem$z*digam2)*con.den )
      c(denzn, s1, s2)
    }, numeric(3)))
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


tolfun <- function(mab = NULL, k = NULL, eps = 1e-7) {
  continue = TRUE
  if (k > 3) {
    mab1 = mab[k-2]
    mab2 = mab[k-1]
    val = abs(abs(mab2) - abs(mab1))
    
    if (val < eps) 
      continue = FALSE
  }
  return(continue)
}











count_binary_fixed_margins <- function(row_sums, col_sums) {
  # ---- helpers ----
  is_zero_vec <- function(v) all(v == 0L)
  
  # Gale–Ryser feasibility (binary) with r (rows) and c (cols) sorted decreasing
  gale_ryser <- function(r, c) {
    if (sum(r) != sum(c)) return(FALSE)
    m <- length(r)
    # c_star(k) = sum_j min(k, c_j)
    c_star <- sapply(seq_len(m), function(k) sum(pmin(k, c)))
    all(cumsum(r) <= c_star)
  }
  
  # Convert column sums to histogram r: r[k] = # of columns with sum == k for k>=1
  # We omit k=0 since those columns can’t take a 1 in any future row.
  col_histogram <- function(c) {
    if (length(c) == 0) return(integer(0))
    maxc <- max(c)
    if (maxc == 0) return(integer(0))
    r <- integer(maxc)
    tab <- table(c)
    for (k in seq_len(maxc)) {
      r[k] <- as.integer(if (as.character(k) %in% names(tab)) tab[as.character(k)] else 0L)
    }
    r
  }
  
  # Reduce operator: r \ s = r - s + L s (L shifts s one to the left: s_{k+1} contributes to block k)
  reduce_hist <- function(r, s) {
    K <- max(length(r), length(s))
    r2 <- integer(K)
    # pad vectors
    rr <- integer(K); rr[seq_along(r)] <- r
    ss <- integer(K); ss[seq_along(s)] <- s
    # r2[k] = (r[k] - s[k]) + s[k+1]   (with s[K+1] := 0)
    for (k in seq_len(K)) {
      add_from_next <- if (k + 1L <= K) ss[k + 1L] else 0L
      r2[k] <- rr[k] - ss[k] + add_from_next
      if (r2[k] < 0L) return(NULL)  # infeasible
    }
    # trim trailing zeros (higher blocks that are empty)
    if (length(r2) > 0L) {
      last <- max(which(r2 != 0L), 0L)
      if (last == 0L) r2 <- integer(0) else r2 <- r2[seq_len(last)]
    }
    r2
  }
  
  # Enumerate all bounded compositions s of total 'need' across K blocks with upper bounds r
  bounded_compositions <- function(need, r) {
    K <- length(r)
    out <- list()
    cur <- integer(K)
    # quick bound: if total available positive columns < need, no compositions
    if (sum(r) < need) return(list())
    # recursion
    rec <- function(i, rem) {
      if (i > K) {
        if (rem == 0L) out[[length(out) + 1L]] <<- cur
        return()
      }
      # upper bound for this block
      ub <- min(r[i], rem)
      # pruning: if remaining capacity from i..K is < rem, stop
      if (sum(r[i:K]) < rem) return()
      for (x in 0:ub) {
        cur[i] <<- x
        rec(i + 1L, rem - x)
      }
    }
    rec(1L, need)
    out
  }
  
  # Safe choose (integer-ish). For typical small cases, base choose() is fine.
  choose_prod <- function(r, s) {
    # product over k of choose(r_k, s_k)
    if (length(r) == 0L) return(1)
    prod(mapply(function(ri, si) choose(ri, si), r, s))
  }
  
  # ---- preprocess & sanity checks ----
  rsum <- as.integer(row_sums)
  csum <- as.integer(col_sums)
  if (any(is.na(rsum)) || any(is.na(csum))) stop("row_sums and col_sums must be integer-like.")
  if (any(rsum < 0) || any(csum < 0)) stop("row_sums and col_sums must be nonnegative.")
  if (sum(rsum) != sum(csum)) return(0)
  
  m <- length(rsum); n <- length(csum)
  if (any(rsum > n) || any(csum > m)) return(0)
  
  # Sort decreasing (does not change the count)
  rsum <- sort(rsum, decreasing = TRUE)
  csum <- sort(csum, decreasing = TRUE)
  
  # Initial feasibility
  if (!gale_ryser(rsum, csum)) return(0)
  
  # Initial histogram of column sums
  r_hist <- col_histogram(csum)  # vector of length max(col_sums) with counts
  
  # Memoization: key is (remaining row sums | r_hist)
  memo <- new.env(parent = emptyenv(), hash = TRUE)
  
  key_of <- function(rs, rh) paste(c(rs, "|", rh), collapse = ",")
  
  # ---- main recursion ----
  f <- function(rs, rh) {
    # If no rows left, all remaining positive column demand must be zero
    if (length(rs) == 0L) {
      return(if (length(rh) == 0L) 1 else 0)
    }
    key <- key_of(rs, rh)
    ans <- get0(key, envir = memo, inherits = FALSE)
    if (!is.null(ans)) return(ans)
    
    need <- rs[1L]                # p1
    rs_rest <- rs[-1L]
    
    # Quick pruners:
    # (i) need cannot exceed number of columns with positive capacity
    if (sum(rh) < need) {
      memo[[key]] <- 0
      return(0)
    }
    # (ii) any column capacity > remaining rows (including current)? If so, impossible
    # Max positive capacity index present equals length(rh)
    if (length(rh) > (length(rs_rest) + 1L)) {
      memo[[key]] <- 0
      return(0)
    }
    
    total <- 0
    # enumerate s in C_r(need)
    comps <- bounded_compositions(need, rh)
    if (length(comps) == 0L) {
      memo[[key]] <- 0
      return(0)
    }
    for (s in comps) {
      # r \ s
      new_rh <- reduce_hist(rh, s)
      if (is.null(new_rh)) next
      # combinatorial factor ∏ choose(r_k, s_k)
      coeff <- choose_prod(rh, s)
      total <- total + coeff * f(rs_rest, new_rh)
    }
    memo[[key]] <- total
    total
  }
  
  f(rsum, r_hist)
}




partitions <- function(n, max.part = n, max.value = n) {
  out <- list()
  recurse <- function(n, max.value, prefix) {
    if (n == 0) {
      out[[length(out) + 1]] <<- prefix
    } else {
      kmax <- min(n, max.value)
      for (k in seq(kmax, 1)) {
        recurse(n - k, k, c(prefix, k))
      }
    }
  }
  recurse(n, max.value, numeric(0))
  # filter by max.part
  out[sapply(out, length) <= max.part]
}


gale_ryser <- function(r, lambda) {
  r <- sort(r, decreasing = TRUE)
  lambda <- sort(lambda, decreasing = TRUE)
  
  m <- length(r)
  n <- length(lambda)
  
  ## For each k, check sum of largest k rows <= sum min(lambda_j, k)
  for (k in seq_len(m)) {
    left  <- sum(r[1:k])
    right <- sum(pmin(lambda, k))
    if (left > right) return(FALSE)
  }
  
  ## Totals must match automatically if lambda is a partition of sum(r)
  TRUE
}



feasible_column_totals <- function(r, n) {
  m <- length(r)
  T <- sum(r)
  
  ## Step 1: generate all partitions of T with ≤ n parts, each ≤ m
  parts <- partitions(T, max.part = n, max.value = m)
  
  ## Step 2: keep only those that satisfy Gale–Ryser
  feas <- list()
  for (p in parts) {
    ## pad partition with zeros to length n
    lambda <- sort(c(p, rep(0, n - length(p))), decreasing = TRUE)
    if (gale_ryser(r, lambda)) {
      feas[[length(feas) + 1]] <- lambda
    }
  }
  
  do.call(rbind, feas)
}



# Generate all ordered integer tuples (x1,...,xk) with 0 <= xi <= ci.
# Returns an integer matrix with prod(ci + 1) rows and k columns.
# Row i is one tuple; columns correspond to x1..xk.
all_integer_tuples <- function(capacities) {
  capacities <- as.integer(capacities)
  if (any(is.na(capacities)) || any(capacities < 0L)) {
    stop("capacities must be nonnegative integers.")
  }
  k <- length(capacities)
  if (k == 0) return(matrix(integer(0), nrow = 1, ncol = 0))  # single empty tuple
  
  # Build a list of ranges 0:ci and use expand.grid efficiently
  ranges <- lapply(capacities, function(c) 0L:c)
  # expand.grid creates a data frame; convert to integer matrix with correct column order.
  df <- do.call(expand.grid, c(ranges, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  # By default, expand.grid puts the first list element as fastest-moving column.
  # If you want x1 to be the left-most column, the list order matches.
  M <- as.matrix(df)
  storage.mode(M) <- "integer"
  colnames(M) <- paste0("x", seq_len(k))
  M
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


qFn <- function(y=NULL) { 
  factorial(length(y))/prod(factorial(as.numeric(table(y)) ))
}

getCounts <- function(y=NULL, x=NULL, r=NULL) {
  num = apply(y,1, function(z) { count_binary_fixed_margins(x, z)*qFn(z)  } )
  z  =  apply(y, 1, function(z) {  as.numeric(table(factor(z, 0:r))) }) 
  rbind(z, num)
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

### THE EMS

logdBBy <- function(y=NULL, n=NULL, alpha=NULL, beta=NULL) {
  lchoose(n, y) + lbeta(y + alpha, n - y + beta) - lbeta(alpha, beta) 
}

dBBy <- function(y=NULL, n=NULL, alpha=NULL, beta=NULL) {
  exp( lchoose(n, y) + lbeta(y + alpha, n - y + beta) - lbeta(alpha, beta) )
}




create.Yloglik <- function(zdat=NULL, nreps=NULL){
  ztots = colSums(zdat)
  zseq = 0:nreps
  function(ab=NULL) {
    sum( ztots*logdBBy(zseq, nreps, ab[1], ab[2] ) )
  }
}



EMY.step <- function(ab=NULL, ztot=NULL, nreps=NULL) {
  ## ztot = the sum_i z_j 
  zseq = 0:(nreps)
  
  a.zseq = digamma( ab[1] + zseq ) - digamma( sum(ab) + nreps )
  b.zseq = digamma( ab[2] + nreps - zseq ) - digamma( sum(ab) + nreps )
  S12 = c( sum(a.zseq*ztot), sum(b.zseq*ztot) )
  
  ## Netwon-Raphson
  ab.new = ab_update(ab, S12)
  
  return(ab.new)
}


EMY.ab <- function(ab=NULL, zdat=NULL,ndat=NULL,fdat=NULL, nreps=NULL, maxiter=100, atol=1e-2) {
  ab = ab[3:4]
  llik = numeric(maxiter)
  loglik.fn = create.Yloglik(zdat, nreps)
  zAvg = colSums(zdat)/sum(zdat)
  zseq = 0:nreps
  mndat = mean(ndat)
  mfdat = mean(fdat)
  
  k = 1
  while (atiken(llik, k = k, eps = atol) & k < maxiter) {
  #while ( tolfun(mab,k,eps = atol) & k < maxiter) {
    # EM-step
    ab = EMY.step(ab, zAvg, nreps)
    
    llik[k] = loglik.fn( ab )
   # mab[k] = alpha_V(par0 = c(mndat,mfdat, ab) )
    k = k + 1
  }
  
  #val = list( loglik = llik[1:(k-1)], mle = c(mean(ndat),mean(fdat),ab) )
  val = list( k=k, mle = c(mndat,mfdat,ab), loglik = llik[k-1] )
  return(val)
}


##N


uniqueCatN <- function(cdat=NULL, ndat=NULL, wt=NULL) {
  nReps = ncol(cdat)
  
  # sort and find the unique pairs of (ci, ni)
  df = data.frame( t(apply(cdat,1,sort, decreasing=TRUE)), n=ndat)
  if (is.null(wt)) uCount = as.matrix(aggregate( data.frame(num=rep(1,nrow(cdat))), df, sum))
  else uCount = as.matrix(aggregate( data.frame(num=wt), df, sum))
  
  # cdat = the cdat  
  # uNum = number of times each values occurs
  # n    =  the value of n for 
  list( cdat   = uCount[,1:nReps], uNum = uCount[,nReps+2], n = uCount[,nReps+1] )
}

create.Nloglik <- function(cdat=NULL, ndat=NULL, wt=NULL) {
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  nReps = ncol(cdat); 
  samSize  = length(ndat);
  sum.ndat = sum(ndat); 
  
  uData  = uniqueCatN(cdat=cdat, ndat=ndat, wt=wt) 
  M      = max(uData$n - apply(uData$cdat,1,max))  
  rCombs = apply(uData$cdat, 1, rankComb) 
  rfx    = getRFX( rCombs, nReps )
  
  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs)
  xdat  = getOrdComb(cMax=max(uData$cdat), nReps=nReps)[1:lcomb,]
  znSet = getZN(lcomb, M, nReps )
  fdat  = xdat
  
  sum.ndat = sum(ndat); 
  zseq = 0:nReps;
  ZNden    = matrix(0, nrow=length(znSet), ncol= (max(uData$cdat)+1) + M )
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    #logN.den = sum.ndat*log( par1[1] ) - samSize*par1[1] 
    logN.den = sum.ndat*log( par1[1] ) - samSize*par1[1] - sum(lfactorial(ndat))
    
    F.den = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    for (i in 1:length(znSet) ) {
      zn = znSet[[i]]
      ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
    }
    
    logC.den = log(vapply(1:length(rfx), function(m=NULL) {
      fxnum = rfx[[m]]
      sum(ZNden[fxnum[,1], uData$n[m]+1]*F.den[fxnum[,2]]*fxnum[,3])
    }, numeric(1) ))   
    
    return( logN.den + sum(uData$uNum*logC.den ) )  
  }
}




create.getFAB <- function(cdat=NULL, ndat=NULL, wt=NULL) {
  nReps = ncol(cdat); 
  samSize  = length(ndat);
  if(is.null(wt)){
    sum.ndat = sum(ndat)
  } else {
    sum.ndat = sum(wt*ndat)
  }
  
  uData  = uniqueCatN(cdat=cdat, ndat=ndat, wt=wt) 
  M      = max(uData$n - apply(uData$cdat,1,max))  
  rCombs = apply(uData$cdat, 1, rankComb) 
  rfx    = getRFX( rCombs, nReps  ) ### list of matrices (fi, xi, num)
  
  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs)
  xdat  = getOrdComb(cMax=max(uData$cdat), nReps=nReps)[1:lcomb,]
  znSet = getZN(lcomb, M, nReps )
  fdat  = xdat
  
  zseq     = 0:nReps;
  maxN1M   = max(uData$cdat)+1 + M
  sum.fdat = rowSums(fdat)
  
  function(tab=NULL) {
    Fden = exp(rowSums(dpois(fdat, tab[1],log=TRUE)))
    
    s12.xn= getZNTable2(tab[2:3], znSet, nReps, maxN1M) 
    
    fs12 = sapply(1:length(rfx), function(m=NULL) {
      fxnum = rfx[[m]]
      if (nrow(fxnum) == 1){ xs12  = matrix(s12.xn[fxnum[,1], uData$n[m] +1,], 1, 3)
      } else {xs12  = s12.xn[fxnum[,1], uData$n[m] +1,]}
      
      xfden = xs12[,1]*Fden[fxnum[,2]]*fxnum[,3]
      cden = sum(xfden)
      
      xf.cden = xfden/cden
      sum.f = sum( sum.fdat[fxnum[,2]] * xf.cden ) 
      
      s12  = colSums(xs12[,-1,drop=FALSE]*xf.cden)
      #c(cden, sum.f, s12)
      c(sum.f, s12)
    } ) 
    sum.fs12 = colSums( t(fs12)*uData$uNum )
    
    #  S12 = sum.fs12[2:3]/sum(uData$n*num ) - digamma( sum(ab) + nReps )
    #  c( sum.fs12[1]/(nReps*sum(num) ), S12)
    S12 = sum.fs12[2:3] - sum.ndat*digamma( sum( sum(tab[2:3]) ) + nReps )
    c( sum.fs12[1], S12)
  }
}






getZNTable2 <- function(ab=NULL, znSet=NULL, nReps=NULL, maxN1M=NULL) {
  ZNden2   = array(0, dim=c(length(znSet), maxN1M, 3) ) 
  
  logdBB = logdBBm(0:nReps, nReps, ab[1], ab[2] )
  digam1 = digamma( ab[1] + 0:nReps )
  digam2 = digamma( ab[2] + nReps - 0:nReps )
  
  for (i in 1:length(znSet) ) {
    zn = znSet[[i]]
    # if (i== 1) stop('here')
    ZNden2[i, zn$n+1, ] = t(vapply( zn$znum, function(tem=NULL) { 
      den = exp(colSums(tem$z*logdBB))*tem$num
      denzn =  sum(den)
      con.den = den/denzn
      s1 = sum( colSums(tem$z*digam1)*con.den )
      s2 = sum( colSums(tem$z*digam2)*con.den )
      c(denzn, s1, s2)
    }, numeric(3)))
  }
  ZNden2
}


getZNTable <- function(ab=NULL, znSet=NULL, nReps=NULL, maxN1M=NULL) {
  nReps2 = nReps+2; zseq = 0:nReps ; nReps3=nReps+3;
  ZNden2   = array(0, dim=c(length(znSet), maxN1M, 3) ); 
  nReps23 = c(nReps2,nReps3);
  
  logdBB = logdBBm(zseq, nReps, ab[1], ab[2] )
  for (i in 1:length(znSet) ) {
    zn = znSet[[i]]
    if (i== 1) stop('here')
    
    denZN = exp( colSums(zn[-nReps23,,drop=FALSE]*logdBB))*zn[nReps2,]
    s1 = rowSums(t(zn[-nReps23,,drop=FALSE]*digamma( ab[1] + zseq ) )) *denZN 
    s2 = rowSums(t(zn[-nReps23,,drop=FALSE]*digamma( ab[2] + nReps - zseq ) )) *denZN 
    
    denZN = aggregate(data.frame(den=denZN, s1=s1, s2=s2), by= data.frame(nSet = zn[nReps3,]), sum)
    denZN[,3:4] = denZN[,3:4]/denZN$den
    
    ZNden2[i, denZN$nSet+1, ] = as.matrix(denZN[,-1])
  }
  ZNden2
}


tolfun <- function(mab = NULL, k = NULL, eps = 1e-2) {
  continue = TRUE
  if (k > 3) {
    mab1 = mab[k-1]
    mab2 = mab[k]
    val = mab2 - mab1
    
    if (val > eps) 
      continue = FALSE
  }
  return(continue)
}



EMN <- function(par=NULL, cdat=NULL, ndat=NULL, maxiter=100, atol=1e-2 ) {
  getFAB = create.getFAB(cdat, ndat)
  
  samSize  = nrow(cdat); 
  nReps    = ncol(cdat)
  sum.n    = sum(ndat)
  
  lam.mle   = mean(ndat); 
  llik      = numeric(maxiter)
  mab      = numeric(maxiter)
  loglik.fn = create.Nloglik(cdat, ndat)
  tab       = par[2:4]
  
  k = 1
  while (atiken(llik, k = k, eps = atol) & k < maxiter) {
    #while ( tolfun(mab,k,eps = atol) & k < maxiter) {
    # EM-step
    fab = getFAB(tab)
    
    tab[1]   = fab[1]/(nReps*samSize)
    tab[2:3] = ab_update(tab[2:3], fab[2:3]/sum.n )
    llik[k] = loglik.fn( c(lam.mle, tab) )
    mab[k] = alpha_V(par0 = c(lam.mle, tab))
    k = k + 1
  }
  
  #val = list( loglik = llik[1:(k-1)], mle = c(lam.mle, tab) )
  val = list(k= k, mle = c(lam.mle, tab), loglik = llik[k-1])
  return(val)
}



####X
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


create.Cloglik <- function(cdat=NULL, M=NULL, wt=NULL) {
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  nReps = ncol(cdat); 
  
  ## sort to find unique c rows
  uCount  = uniqueCdat(cdat=cdat, wt=wt) 
  
  ## The max possible C 
  cMax   = max(cdat); nReps = ncol(cdat); 
  rCombs = apply(uCount$cdat, 1, rankComb) 
  rfx   = getRFX( rCombs, nReps )
  
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
      ZNden[i, zn$n+1] = vapply( zn$zn, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1))
    }
    
    logC.den = log(sapply(rfx, function(fxnum=NULL) {
      sum( colSums(ZNden[fxnum[,1],,drop=FALSE]*F.den[fxnum[,2]]*fxnum[,3])*dpois( nSeq, par1[1])  )
    } ))   
    return( sum(uCount$cNum*logC.den) )  
  }
}


















create.getNFAB <- function(cdat=NULL, M=4,  wt=NULL) {
  nReps = ncol(cdat); 
  
  uData  = uniqueCdat(cdat=cdat, wt=wt) 
  rCombs = apply(uData$cdat, 1, rankComb) 
  rfx    = getRFX( rCombs, nReps ) ### list of matrices (fi, xi, num)
  
  ## The largest combination required is max(rCombs)
  lcomb = max(rCombs)
  xdat  = getOrdComb(cMax=max(uData$cdat), nReps=nReps)[1:lcomb,]
  znSet = getZN(lcomb, M, nReps)
  fdat  = xdat
  
  zseq     = 0:nReps;
  maxN1M   = max(uData$cdat)+1 + M 
  nSeq     = 0:(maxN1M -1)
  sum.fdat = rowSums(fdat)
  
  function(par1=NULL) {
    Fden    = exp(rowSums(dpois(fdat, par1[2],log=TRUE)))
    znTable = getZNTable2(par1[c(3,4)], znSet, nReps, maxN1M) 
    
    cfs12 = sapply(rfx, function(fxnum=NULL) {
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


EMC <- function(par1=NULL, cdat=NULL, maxiter=100, atol=1e-2, M=4) {
  
  getNFAB = create.getNFAB(cdat, M)
  samSize  = nrow(cdat); nReps    = ncol(cdat);
  
  loglik.fn = create.Cloglik(cdat, M)
  llik      = numeric(maxiter)
  mab      = numeric(maxiter)
  
  k = 1
  
  while (atiken(llik, k = k, eps = atol) & k < maxiter) {
    #while ( tolfun(mab,k,eps = atol) & k < maxiter) {
    # EM-step
    NFAB = getNFAB(par1)
    
    lambda.new = NFAB[1]/(samSize )
    theta.new  = NFAB[2]/(nReps*samSize )
    ab.new     = ab_update(par1[3:4], NFAB[3:4]/NFAB[1] )
    #ab.new = optim(par1[3:4],fn_ab_nu_eta, method = "BFGS", S12 = NFAB[3:4]/NFAB[1])$par
    par1 = c(lambda.new, theta.new, ab.new)
    
    llik[k] = loglik.fn( par1)
    #mab[k] = par1[3]/(par1[3] + par1[4])
    mab[k] = alpha_V(par0 = par1)
    k = k + 1
  }
  # par1[3:4] = c(exp(ab.new[1]), exp(ab.new[2]) +1 )
  
  
  val = list( k = k, mle = par1,loglik = llik[(k-1)] )
  #val = list( loglik = llik[1:(k-1)], mle = par1)
  return(val)
}

###X


Zn.N <- function(xdat=NULL, ndat=NULL,M=4) {
  nReps = ncol(xdat)
  #numX = nrow(xdat)
  val = list()
  for(i in 1:nrow(xdat)){
    n = ndat[i]
    xdatk = xdat[i,]
    tots = feasible_column_totals(xdatk, n)
    z   = apply(tots, 1, function(z) {  as.numeric(table(factor(z, 0:nReps ))) })
    num = apply(tots, 1, function(z) { count_binary_fixed_margins(xdatk, z)*qFn(z)  } )
    val[[i]] = list( znum= z, num =num)
  }
  val
}



create.Xloglik1 <- function(xdat=NULL, fdat=NULL, ndat=NULL) {
  # xdat is a n x r matrix  # fdat is a n x r matrix  # ndat is a n x 1 vector
  
  # some useful variables
  nReps = ncol(xdat); nReps2 = ncol(xdat)+2; zseq = 0:nReps ;
  samSize = length(ndat); sum.fdat = sum(fdat);  sum.ndat = sum(ndat);
  zseq = 0:nReps
  #convert xdat -> ydat -> zdat
  #zdat = sapply(1:samSize, function(m=NULL) { 
  #    y = feasible_column_totals(xdat[m,], ndat[m]) 
  #    getCounts(y, xdat[m,], r=nReps )
  # }  )
  zn.dat = Zn.N(xdat, ndat, M=4)
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    if(par1[2]==0) par1[2] = par1[2] + 1e-8
    logF.den =  sum.fdat*log( par1[2] ) - nReps*samSize*par1[2] - sum(log(factorial(fdat)))
    logN.den =  sum.ndat*log( par1[1] ) - samSize*par1[1] - sum(log(factorial(ndat)))
    
    # ( (maxN +1)  x nrow(x)   ) matrix for f( x | n) 
    #logX.den = log(unlist(lapply(zdat, function(z=NULL) {
    #  sum(exp( colSums(z[-nReps2,,drop=FALSE]*logdBBm(zseq, nReps, par1[3], par1[4] )))*z[nReps2,] )
    #} )) )   
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    logX.den = log(vapply( zn.dat, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1)) ) 
    
    return( sum(logN.den) + sum(logF.den ) + sum(logX.den ) )  
  }
}



EMX.step <- function(ab=NULL, zdat=NULL, ndat=NULL, nReps=NULL, samSize=NULL) {
  ## ztot = the sum_i z_j 
  nReps2 = nReps+2; zseq = 0:nReps ;
  
  sumS = sapply(1:samSize, function(k) {
    z=zdat[[k]]
    # print(k)
    if (ncol(z$znum) > 1 ) {
      #stop('here')
      z.cond.den = exp(colSums( z$znum*logdBBm(zseq, nReps, ab[1], ab[2] ) ))*z$num
      z.cond.den = z.cond.den/sum(z.cond.den)
      
      s1 = sum( colSums( t(z$znum*digamma( ab[1] + zseq ))*z.cond.den  ) )
      #s2 = sum( colSums( t(z[-nReps2,]*digamma( ab[2] + nReps - zseq )) )*z.cond.den  )
      s2 = sum( colSums( t(z$znum*digamma( ab[2] + nReps - zseq ))*z.cond.den ) )
    } else {
      s1 = sum( z$znum*digamma( ab[1] + zseq ) )
      s2 = sum( z$znum*digamma( ab[2] + nReps - zseq ) )
    }
    
    c(s1,s2 )
  } )  
  ## need add -digamma( sum(ab)) 
  S12 = rowSums(sumS)/sum(ndat) - digamma( sum(ab) + nReps )
  #S12 = rowSums(sumS) - ndat*digamma( sum(ab) + nReps )
  ## Netwon-Raphson
  ab.new = ab_update(ab, S12)
  
  return(ab.new)
}

EMX.ab <- function(par=NULL, xdat=NULL, fdat=NULL, ndat=NULL, maxiter=100, atol=1e-2) {
  
  # this is add hoc 
  maxN = max( qpois(1-1e-4, mean(ndat) ), max(xdat) +2 )
  nSeq   = 0:maxN  # for lambda = 30 this is 52
  nReps = ncol(xdat); nReps2 = ncol(xdat)+2; zseq = 0:nReps ;
  ab = par[3:4]
  samSize = length(ndat);
  
  #xdat -> ydat -> zdat
  # zdat = sapply(1:samSize, function(m=NULL) { 
  #  y = feasible_column_totals(xdat[m,], ndat[m]) 
  #  getCounts(y, xdat[m,], r=nReps )
  #}  )
  
  
  zdat = Zn.N(xdat, ndat, M=4)
  
  lam.mle   = mean(ndat); theta.mle = mean(fdat);
  llik = numeric(maxiter)
  loglik.fn = create.Xloglik1(xdat = xdat, fdat = fdat, ndat = ndat)
  
  k = 1
  while (atiken(llik, k = k, eps = atol) & k < maxiter) {
    #while ( tolfun(mab,k,eps = atol) & k < maxiter) {
    # EM-step
    ab = EMX.step(ab = ab, zdat=zdat, ndat=ndat, nReps=nReps, samSize= samSize) 
    
    llik[k] = loglik.fn( c(lam.mle, theta.mle, ab) )
    k = k + 1
  }
  
  #val = list( loglik = llik[1:(k-1)], mle = c(lam.mle, theta.mle, ab) )
  val =  list( k = k, mle = c(lam.mle, theta.mle, ab), loglik = llik[k-1] )
  return(val)
}










##R 

create.Rloglik <- function(xdat=NULL, fdat=NULL, M=5, wt=NULL) {
  # xdat is a n x r matrix  and fdat is a n x r matrix  
  # ni  goes from  max(xi1,...,xir) to max(xi1,...,xir) + M
  
  df = data.frame( t(apply(xdat,1,sort, decreasing=TRUE)) )
  if (is.null(wt)) {xCount = as.matrix(aggregate( data.frame(num=rep(1,nrow(xdat))), df, sum))
  } else {xCount = as.matrix(aggregate( data.frame(num=wt), df, sum))}
  
  # Some useful variables
  nReps = ncol(xdat)
  zseq = 0:nReps ; 
  samSize = nrow(fdat); 
  sum.fdat = sum(fdat);  
  
  
  #xdat -> ydat -> zn.dat 
  rXCombs = apply(xCount[,1:nReps], 1, rankComb) 
  zn.dat  = getZNSet(combSet=rXCombs, M=M,nReps )
  xNum    = xCount[,nReps+1]
  
  function(par1=NULL)  {
    # par1 is the canonical parameters and par1 = (lambda, theta, alpha, beta)
    logF.den =  sum.fdat*log( par1[2] ) - nReps*samSize*par1[2] - sum(lfactorial(fdat))
    
    # ( (maxN +1)  x nrow(x)   ) matrix for f( x | n) 
    logdBB = logdBBm(zseq, nReps, par1[3], par1[4] )
    logX.den = log(vapply(zn.dat, function(zn=NULL) {
      xden.n = vapply( zn$znum, function(tem=NULL) { sum( exp(colSums(tem$z*logdBB))*tem$num) }, numeric(1)) 
      sum( xden.n*dpois(zn$n, par1[1]) )
    }, numeric(1) ) ) 
    ind = which(logX.den!="-Inf")
    return( sum(logF.den ) + sum( xNum[ind]*logX.den[ind]) )  
  }
}

getZNSet <- function(combSet=NULL, M=NULL,nReps) {
  #znSet <- readRDS("znSet5.rds")
  #znSet <- readRDS("znSetx12r3M10.rds")
  #if (nReps ==3) {znSet <- readRDS(paste0("znSetx30r", nReps, "M10.rds", collpase=""))
  #} else znSet <- readRDS(paste0("znSetx25r", nReps, "M6.rds", collpase=""))
  
  znSet = readRDS(paste0("znSetx12r4M5.rds", collpase=""))
  
  
  znSet = znSet[combSet]
  znSet = lapply(znSet, function(z=NULL, M=NULL ) {
    nSet = z$n <= min(z$n) + M
    z$znum = z$znum[nSet]
    z$n    = z$n[nSet]
    z
  }, M=M )
  znSet
}



create.getNAB <- function(xdat=NULL, M=NULL, wt=NULL) {
  nReps = ncol(xdat)
  df = data.frame( t(apply(xdat,1,sort, decreasing=TRUE)) )
  if (is.null(wt)) {xCount = as.matrix(aggregate( data.frame(num=rep(1,nrow(xdat))), df, sum))
  } else {xCount = as.matrix(aggregate( data.frame(num=wt), df, sum))}
  
  rXCombs = apply(xCount[,1:nReps], 1, rankComb) 
  zn.dat  = getZNSet(rXCombs, M=M, nReps)
  #lcomb = max(rXCombs)
  #zn.dat = getZN(lcomb, M, nReps)
  xNum    = xCount[,nReps+1]
  zseq    = 0:nReps;
  
  function(lab=NULL) {
    ## (lambda, alpha, beta)
    
    logdBB = logdBBm(zseq, nReps, lab[2], lab[3] )
    digam1 = digamma( lab[2] + zseq )
    digam2 = digamma( lab[3] + nReps - zseq )
    
    ns12 = vapply(zn.dat, function(zn=NULL) {
      ds12.n = t(vapply( zn$znum, function(tem=NULL) { 
        den = exp(colSums(tem$z*logdBB))*tem$num
        cden = den/sum(den)
        s1 = sum( colSums(tem$z*digam1)*cden )
        s2 = sum( colSums(tem$z*digam2)*cden )
        c(sum(den), s1, s2)
      }, numeric(3)))
      ## ds12.n is ( length(zn$n) x 3) matrix
      denzn = ds12.n[,1]*dpois(zn$n, lab[1])
      colSums(cbind(zn$n, ds12.n[,-1])*denzn)/sum(denzn)
    }, numeric(3))
    
    sum.ns12 = colSums( t(ns12)*xNum,na.rm=TRUE )
    
    S12 = sum.ns12[2:3] - sum.ns12[1]*digamma( sum( lab[2:3] ) + nReps )
    c( sum.ns12[1], S12)
  }
}




EMR <- function(par=NULL, xdat=NULL, fdat=NULL, M=4, maxiter=100, atol=1e-2) {
  
  getNAB = create.getNAB(xdat, M)
  
  samSize   = nrow(fdat); 
  theta.mle = mean(fdat);
  llik      = numeric(maxiter)
  loglik.fn = create.Rloglik(xdat, fdat, M)
  lab       = c(par[1], par[3:4])
  
  k = 1
  while (atiken(llik, k = k, eps = atol) & k < maxiter) {
    #while ( tolfun(mab,k,eps = atol) & k < maxiter) {
    # EM-step
    NAB = getNAB(lab)
    lab = c(NAB[1]/samSize, ab_update(lab[2:3], NAB[2:3]/NAB[1]) )
    
    llik[k] = loglik.fn( c(lab[1], theta.mle, lab[2:3]) )
    #mab[k] = alpha_V(par0 = c(lab[1], theta.mle, lab[2:3]) )
    k = k + 1
    #print(k)
  }
  
  #val = list( loglik = llik[1:(k-1)], mle = c(lab[1], theta.mle, lab[2:3]) )
  val = list(k= k, mle = c(lab[1], theta.mle, lab[2:3]), loglik = llik[k-1] )
  return(val)
}

###FOR SP
## ALPHA ###############
transformation.alphaM <- function(par=NULL) {
  lambda = par[1] 
  p      = par[2]
  theta  = par[3]
  
  mu = lambda*p + theta
  val = numeric(3)
  
  val[1] = (p^2*mu - lambda*p^3)/mu^2
  val[2] = ( 2*lambda*p*mu - (lambda*p)^2 )/mu^2
  val[3] = - lambda*p^2/mu^2
  
  val 
}

var.alphaM <- function(fisherInfpar=NULL, par=NULL) {
  dpar = transformation.alphaM(par)
  
  as.numeric( t(dpar) %*% solve(fisherInfpar) %*% (dpar) )
}

### MU
getmu <- function(par=NULL) { par[1]*par[2]+par[3] }


####### FOR EMs ########
atiken <- function(loglik = NULL, k = NULL, eps = 1e-2) {
  continue = TRUE
  if (k > 3) {
    lm1 = loglik[(k-1)]
    lm = loglik[(k - 2)]
    lm_1 = loglik[(k - 3)]
    
    am = (lm1 - lm)/(lm - lm_1)
    lm1.Inf = lm + (lm1 - lm)/(1 - am)
    val = lm1.Inf - lm
    
    #if (is.na(val)) continue = FALSE
    
    if (val < eps & val >= 0 | is.na(val)) 
      continue = FALSE
  }
  return(continue)
}

##SINGLE PHASE EM ##########################
denCTN00 <- function(kseq=2:3, N=2, p=0.95, theta=2) {
  ## calculates f(c=kseq| n=N)
  den = outer(0:N, kseq, function(x,y) {dbinom(x, size=N, prob=p)*dpois(y-x, lambda=theta) })
  return( colSums(den) )
}

create.wtloglikFIXEDp <- function(cdat=NULL, maxN=30) {
  df = data.frame( t(apply(cdat,1,sort)) )
  df = aggregate(list(num=rep(1,nrow(df))), df, length)
  wt = df$num
  countdata = df[,1:ncol(cdat)]
  
  #maxN = qpois(1-1e-4, mean(cdat) ) # for lambda = 30 this is 52
  nSeq   = 0:maxN  # for lambda = 30 this is 52
  cMin   = min(countdata)
  cSeq   = cMin:max(countdata) 
  
  function(par1=NULL)  {
    # convert the canonical parameters 
    par =par1
    
    # (maxN x cMax) matrix for f( ci | n) is the same for all ci=k
    dci.n = sapply(nSeq, function(k) { denCTN00(cSeq, N=k, p=par[2], theta=par[3]) } )
    
    # (maxN x nrow(cdat)) matrix of joint densities f(c1,...ck| n ) for each row 
    dcVec.n = apply(countdata, 1, function(z) { 
      # ( 1 - cMin) converts the counts to row index 
      apply(dci.n[ z + (1 - cMin),], 2, prod)
    } )
    
    dCvec = colSums( dcVec.n*dpois(nSeq, lambda=par[1]) )
    
    return( sum(wt*log(dCvec)) )  
  }
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


###fixed functions

## Create reinspection loglik
create.rloglik1 <- function(xdat=NULL, fdat=NULL) {
  # xdat is a n x r matrix
  # fdat is a n x r matrix
  
  # this is add hoc 
  maxN = max( qpois(1-1e-4, mean(xdat + fdat) ), max(xdat) +2 )
  nSeq   = 0:maxN  # for lambda = 30 this is 52
  
  num.fdat = length(fdat)
  sum.fdat = sum(fdat)
  
  function(par=NULL)  {
    # par is the canonical parameters 
    if(par[3]==0) par[3] = par[3] + 1e-8
    
    logf.den =  sum.fdat*log( par[3] ) - num.fdat*par[3] - sum(log(factorial(fdat)))
    #logf.den = sum( dpois(fdat, lambda=par[3], log=TRUE) )
    
    # ( (maxN +1)  x nrow(x)   ) matrix for f( x | n) 
    dx.n = t(sapply(nSeq, function(k) { rowSums(dbinom( xdat, size=k, prob=par[2], log=TRUE)) } ))
    
    dx = colSums( exp(dx.n +dpois(nSeq, lambda=par[1], log=TRUE)) )
    
    return( sum(log(dx)) + logf.den )  
  }
}




EMr.step <- function(par=NULL, xdat=NULL, maxN=30) {
  ## EM step for reinspection plan
  nSeq  = 0:maxN  
  
  
  dx.n = t(sapply(nSeq, function(k) { rowSums(dbinom( xdat, size=k, prob=par[2], log=TRUE)) } ))
  dxn  = exp(dx.n + dpois(nSeq, lambda=par[1], log=TRUE))
  dx = colSums( dxn )
  
  expn.x = rowSums( t(nSeq * dxn)/dx )
  
  lambda.hat = mean(expn.x)
  p.hat = sum(xdat)/( ncol(xdat)*sum(expn.x))
  
  return(c(lambda.hat, p.hat))
}

## log likelihood for MLEs for TC

create.loglikn2 <- function(cdat=NULL,ndat=NULL) {
  # par = (beta, eta)
  #tkval = table(kdata)
  #val = as.numeric( names(tkval) )
  #freq = as.numeric( tkval )
  cdat=cdat
  ndat=ndat
  
  function(par=NULL) {
    #-1*(sum( freq*log(denCTN(kseq= val, N=N, p=par[2], theta=par[1]))	) + sum(log(dpois(ni,lambda=par[3]))))
    #-1*sum(sapply(ni, function(N) {sum( freq*log(denCTN(kseq= val, N=N, p=par[2], theta=par[1]))	) } )) +
    # -1*sum(log(dpois(ni,lambda=par[3])))
    -1*sum( log( sapply(1:nrow(cdat), function(i){
      denCTN2(cdat =cdat[i,],n=ndat[i], p=par[2], theta=par[3])*dpois(ndat[i],par[1])
    })) )
    
  }
}



###send in one obs (r measurements) and its particular n
denCTN2 <- function(cdat=NULL, n=NULL, p=NULL, theta=NULL) {
  den = prod(sapply(cdat, function(k) {
    if (n - k >= 0) {
      xseq = 0:min(n,k)
    } else {
      xseq = 0:n
    }
    sum(dbinom(xseq, n,p)*dpois(k-xseq,theta))				
  }))
  
  den
}



######### LOG LIK GOLD STANDARD
create.loglikGS <- function(xdat=NULL,fdat=NULL, ndat = NULL) {
  
  cdat=cdat
  ndat=ndat
  fdat=fdat
  
  function(par1=NULL) {
    #-1*(sum( freq*log(denCTN(kseq= val, N=N, p=par[2], theta=par[1]))	) + sum(log(dpois(ni,lambda=par[3]))))
    #-1*sum(sapply(ni, function(N) {sum( freq*log(denCTN(kseq= val, N=N, p=par[2], theta=par[1]))	) } )) +
    # -1*sum(log(dpois(ni,lambda=par[3])))
    sum( log( sapply(1:nrow(xdat), function(i){
      dpois(ndat[i], par1[1])*prod(dbinom(xdat[i,], size = ndat[i], prob = par1[2])*dpois(fdat[i,], par1[3]))
    })) )
    
  }
}

create.loglikGSp <- function(xdat=NULL,fdat=NULL, ndat = NULL) {
  
  cdat=cdat
  ndat=ndat
  fdat=fdat
  
  function(par1=NULL) {
    sum( log( sapply(1:nrow(xdat), function(i){
      prod(dbinom(xdat[i,], size = ndat[i], prob = par1[2]))
    })) )
    
  }
}






#transformation so we dont need to use constrained optimization
par21 <- function(par=NULL) {
  ## go from par2 to par1 
  par[1] = exp(par[1])
  par[2]   = exp(par[2])/(1+exp(par[2]))
  par[3] = exp(par[3])
  return(par)
}
par12 <- function(par=NULL) {
  ## go from par1 to par2
  par[1] = log(par[1])
  par[2]   = log(par[2]/(1-par[2]))
  par[3] = log(par[3])
  return(par)
}



EM.r <- function(ipar=NULL, xdat=NULL, fdat=NULL, maxiter=100, atol=1e-2) {
  ## EM for reinspection plan
  
  # this is add hoc 
  maxN = max( qpois(1-1e-4, mean(xdat + fdat) ), max(xdat) +2 )
  mle.theta = mean(fdat)
  mle12  = ipar[1:2]
  
  llik = numeric(maxiter)
  loglik.fn = create.rloglik1(xdat, fdat)
  
  k = 1
  while (atiken(llik, k = k, eps = atol) & k < maxiter) {
    # EM-step
    mle12 = EMr.step(par=mle12, xdat=xdat, maxN=maxN )
    
    mle = c(mle12, mle.theta)
    llik[k] = loglik.fn( mle )
    k = k + 1
  }
  
  val = list( loglik = llik[(k-1)], mle = mle)
  return(val)
}




