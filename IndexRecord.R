
#-------------------------------------------------------------------------------
# functions from R packages 'fda' and 'fdaoutlier'. These are references on 
# the main paper

combinat=function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}
#BD2
fBD2=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=apply(rmat,1,min)-1
  up=n-apply(rmat,1,max)
  (up*down+n-1)/combinat(n,2)
}
#MBD
fMBD=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  (rowSums(up*down)/p+n-1)/combinat(n,2)
}

extremal_depth <- function(dts){
  if(is.data.frame(dts)){
    dts <- as.matrix(dts)
  }
  
  if(!is.array(dts) || !is.numeric(dts))
    stop("Argument \"dts\" must be a numeric matrix or dataframe.")
  
  if (any(!is.finite(dts))){
    stop("Missing or infinite values are not allowed in argument \"data\"")
  }
  if(nrow(dts) < 3) stop("The number of curves must be greater than 2")
  
  ddim <- dim(dts)
  n <- ddim[1]
  p <- ddim[2]
  pwdepth <- pwise_depth(dt = dts, n = n) # matrix of n by p
  pmfs <- apply(pwdepth, 1, function(x){
    pmf <- table(x)/p
    return(c(as.numeric(names(pmf[1])), # depth level and mass
             pmf[1]))
  })
  
  depth_levels <- pmfs[1,]
  masses <- pmfs[2, ]
  # order functions according to depth_levels and mass
  ordered_functions <- sapply(sort(unique(depth_levels),
                                   method = "quick"),
                              function(x){
                                fns_depth_level <- which(depth_levels == x)
                                if(length(fns_depth_level) > 1){
                                  fns_depth_level[order(masses[fns_depth_level], decreasing = T)]
                                }else{
                                  fns_depth_level
                                }
                              })
  ordered_functions <- unlist(ordered_functions)
  depth_values <- ((1:n)/n)[order(ordered_functions)]
}

pwise_depth <- function(dt, n) {
  pdepth <- apply(dt, 2, function(i){
    (1 - abs(2*rank(i) - n -1)/n) # for rank r, d = ((r - 1) + (n - r))/n
  })
  return(pdepth)
}

#computing time records --------------------------------------------------------

IndexRecord <- function(data, depth=c('MBD','ED'), plot=FALSE){
  # It computes the index times where record curves are observed
  # args:
  #   data: a matrix with dim m x n, where m is the total num of points observed 
  #         for a function, and n is the total number of functions
  #   depth: depth to be use to order curves and estimate the records
  #
  # values: list 
  #   Est: a estatistic value 
  #   Records: Index/time where records are estimated
  #   UpperR: Index/time where upper record curves are estimated
  #   LowerR: Index/time where lower record curves  are estimated
  method <- match.arg(depth)
  N <- dim(data)[2]
  p <- dim(data)[1]
  Indx.r.U <- Indx.r.L <-  rep(0,N) # indicators where records are observed
  Yu <- Yl <-  rep(0,N)  # indx process of upper/lower recrods (sumsum of Indx.r.U)
  indexminmax <- matrix(NA,N,2)  # a matrix containing the two most extreme curves at each time
  indexminmax[1,] <- indexminmax[2,] <- c(1,2)   # a matrix containing the two most extreme curves at each time
  # classification of the two first records
  mean0 <-  rowMeans(data[,1:2])
  # is the first curve a lower/upper record?
  if( sum(data[,1]< mean0)>floor(p/2) ){
    Indx.r.L[1] <- 1
    Yl[1:2] <- 1
  }else{
    Indx.r.U[1] <- 1
    Yu[1:2] <- 1
  }
  # is the second function a lower/upper record?
  if( sum(data[,2]< mean0) > floor(p/2) ){
    Indx.r.L[2] <- 1
    Yl[1:2] <- 2
  }else{
    Indx.r.U[2] <- 1
    Yu[1:2] <- 2
  }
  
  if(method=='MBD'){
    for (k in 3:N) {
      new.data=data[,1:k]
      depth=fMBD(new.data)
      index=order(depth) # from the smallest to the biggest
      extreme.index <- index[1:2]
      
      if(k%in%extreme.index){
        mean0 <-  rowMeans(new.data)
        if(sum(data[,k]< mean0) > floor(p/2)){
          # k is candidate for lower record
          Indx.r.L[k] <- 1*(depth[k] < depth[Yl[k-1]])
          Yl[k] <- k*Indx.r.L[k] + Yl[k-1]*1*(depth[k] >= depth[Yl[k-1]])
          Yu[k] <- Yu[k-1]
        }else{
          # k is candidate for upper record
          Indx.r.U[k] <- 1*(depth[k] < depth[Yu[k-1]])
          Yu[k] <- k*Indx.r.U[k] + Yu[k-1]*1*(depth[k] >= depth[Yu[k-1]])
          Yl[k] <- Yl[k-1]
        }
      }else{
        Yu[k] <- Yu[k-1]
        Yl[k] <- Yl[k-1]
      }
      indexminmax[k,] <- extreme.index[order(extreme.index)]
    }
  }else{
    # Continuar aqui................................
    for (k in 3:N) {
      new.data=data[,1:k]
      depth= extremal_depth(t(new.data))
      index=order(depth) # from the smallest to the biggest
      extreme.index <- index[1:2]
      
      if(k%in%extreme.index){
        mean0 <-  rowMeans(new.data)
        if(sum(data[,k]< mean0) > floor(p/2)){
          # k is candidate for lower record
          Indx.r.L[k] <- 1*(depth[k] < depth[Yl[k-1]])
          Yl[k] <- k*Indx.r.L[k] + Yl[k-1]*1*(depth[k] >= depth[Yl[k-1]])
          Yu[k] <- Yu[k-1]
        }else{
          # k is candidate for upper record
          Indx.r.U[k] <- 1*(depth[k] < depth[Yu[k-1]])
          Yu[k] <- k*Indx.r.U[k] + Yu[k-1]*1*(depth[k] >= depth[Yu[k-1]])
          Yl[k] <- Yl[k-1]
        }
      }else{
        Yu[k] <- Yu[k-1]
        Yl[k] <- Yl[k-1]
      }
      indexminmax[k,] <- extreme.index[order(extreme.index)]
    }
  }
  Ind.Record <- Indx.r.L+Indx.r.U # upper and lower record 
  timeC <- 1/sqrt(1:N)
  Ind.S=timeC*cumsum(Ind.Record)
  
  if(plot){
    
    ## modificar esto
    par(mfrow=c(1,2))
    plot(1:length(Ind.Record), cumsum(Ind.Record), type = 's', col=1,
         ylab = 'num of records', xlab = 'index of time', main=paste('Record times with', method))
    lines(1:length(Indx.r.L), cumsum(Indx.r.L), type = 's', col=4 )
    lines(1:length(Indx.r.U), cumsum(Indx.r.U), type = 's', col=2 )
    legend("topleft", legend = c('total records','upper records', "lower records" ), col = c(col=1,2,4),
           ncol = 1, cex = 1, lwd = 2, bty='n')
    matplot(data, type = 'l', col = grey(.7,.4), main='Functional Data', xlab = 's')
    matplot(data[,Indx.r.U==1], type = 'l', col = 2, add = TRUE )
    matplot(data[,Indx.r.L==1], type = 'l', col = 4, add = TRUE )
    par(mfrow=c(1,1))
  }
  
  return(list(Est=Ind.S, Records=Ind.Record, UpperR=Indx.r.U, 
              LowerR=Indx.r.L, IndexSupInf=indexminmax) )
}

