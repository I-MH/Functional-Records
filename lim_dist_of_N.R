

lim_dist_Up_records <- function(x){
  # limit dist of functional upper records for a functional random walk
  sqrt(1/pi)*exp(- (x^2)/4 )
}


find_quantile_lim_dist_Up <- function(quantile = 0.95){
  # quantile of lim_dist_Up_records
  xx <- seq(.01,4, by=0.005 )
  quant <- c()
  for (j in 1:length(xx)) {
    quant[j] <- integrate(lim_dist_Up_records, lower=0, upper=xx[j])$value
  }
  return(xx[which(quant >= quantile)[1]])
}


lim_dist_records <- function(x){
  # limit dist of functional records for a functional random walk
  sqrt(2/pi)*( x^2 )*exp(- x^2 / 2 )
}

find_quantile_lim_dist_r <- function(quantile = 0.95){
  # quantile of lim_dist_records
  xx <- seq(.1,5, by=0.005 )
  quant <- c()
  for (j in 1:length(xx)) {
    quant[j] <- integrate(lim_dist_records, lower=0, upper=xx[j])$value
  }
  return(xx[which(quant >= quantile)[1]])
}

if(FALSE){
  # 
  find_quantile_lim_dist_Up( quantile = 0.005)  # 0.01
  find_quantile_lim_dist_Up( quantile = 0.01)   # 0.02
  find_quantile_lim_dist_Up( quantile = 0.05)   # 0.09
  find_quantile_lim_dist_Up( quantile = 0.1)    # 0.18
  find_quantile_lim_dist_Up( quantile = 0.5)    # 0.955
  # quantiles
  quantilesUp <- c(0.01, 0.02, 0.09, 0.18, 0.955) # c(.5,1,5,10,50) %
  # 
  find_quantile_lim_dist_r( quantile = 0.005)  # 0.27
  find_quantile_lim_dist_r( quantile = 0.01)   # 0.34
  find_quantile_lim_dist_r( quantile = 0.05)   # 0.595
  find_quantile_lim_dist_r( quantile = 0.1)    # 0.765
  find_quantile_lim_dist_r( quantile = 0.5)    # 1.54
  # quantiles
  quantiles <- c(0.27, 0.34, 0.595, 0.765, 1.54) # c(.5,1,5,10,50) %
}




