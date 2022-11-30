# it uses 
#source('/R/IndexRecord.R')
#source('/R/lim_dist_of_N.R')

RB.unit.root.test <- function(data, total.r=NULL, depth='MBD'){
  # Test for a functional unit root
  # args:
  #   data: a matrix with dim m x n, where m is the total num of points observed 
  #         for a function, and n is the total number of functions
  #   total.r: total number of functional records of the data
  #   depth: type of depth to be used if total.r=NULL to compute the records
  # values: result of the test
  
  N <- dim(data)[2]
  if(is.null(total.r)){
    est.records <- IndexRecord(data, depth= depth)
    total.r <- sum(est.records$Records)
  }
  Tstat <- round(total.r/sqrt(N),4)
  p_values <- round(integrate(lim_dist_records, lower=0, upper=Tstat)$value,
                    4)
  
  cat("\n")
  cat("Unit root test based on records for a functional time series")
  cat("\n")
  cat("Null hypothesis: the series has a functional random walk component\n")
  cat("\n")
  cat(paste("Test statistic value = ", Tstat, sep = ""))
  cat("\n")
  cat(paste("p-value = ", p_values, sep = ""), "\n")
  
  cat(paste("Total number of recrods = ", total.r, sep = ""), 
      "\n")
}
