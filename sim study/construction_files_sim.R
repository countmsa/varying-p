



construct.ZN <- function(xdat=NULL, M=NULL) {
  nReps = ncol(xdat)
  #numX = nrow(xdat)
  zndat = apply(xdat, 1, function(xdatk) {
    nSeq = max( xdatk):( max( xdatk ) +M )
    val = lapply(nSeq, function(n=NULL) {
      #cbind(t(getCounts(feasible_column_totals(xdat[k,], n) , xdat[k,], r=nReps )),n)
      tots = feasible_column_totals(xdatk, n)
      list( z   = apply(tots, 1, function(z) {  as.numeric(table(factor(z, 0:nReps ))) }), 
            num = apply(tots, 1, function(z) { count_binary_fixed_margins(xdatk, z)*qFn(z)  } ) )
    })
    list( znum= val, n=nSeq)
  }) 
  zndat
}


construct.RFX <- function(cSet=NULL) {
  ## Given a set of c find (x,f) such that x+f=c
  ## output the indices for x and f
  
  nReps = ncol(cSet)
  rfx = sapply(1:nrow(cSet), function(m=NULL) { 
    xdat = all_integer_tuples(cSet[m,]) 
    fdat = matrix(cSet[m,], nrow=nrow(xdat), ncol=ncol(xdat), byrow=TRUE) - xdat
    
    rfx = aggregate(list(num=rep(1, nrow(fdat))), 
                    by= data.frame(rx=apply(xdat, 1, rankComb), rf=apply(fdat, 1, rankComb)), sum)
    as.matrix(rfx)
  }  )
  rfx
}

#znSetx15r5M5 = construct.ZN(xdat=getOrdComb(cMax=15, nReps=5), M=5)

###30 
#cset12 = all_integer_tuples(c(12,12,12)) 
#uCount12  = uniqueCdat(cdat=cset12) 

#nReps = ncol(cset12); 
#rCombs12 = apply(uCount30$cdat, 1, rankComb) 

#test.results = construct.RFX(uCount12$cdat)

#for (i in 1:length(test.results)){
#  subs = test.results[[i]]
#  j = rCombs12[i]
#  saveRDS(subs, file= paste("rfx12r", nReps, "/", "rfx12r", nReps, "c", j, ".rds", sep=''))
#}

