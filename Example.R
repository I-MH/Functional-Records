
# load IndexRecord.R
require(demography)
require(ggplot2)

years <- fr.mort$year
ages <- fr.mort$age[1:101]
france.fit <- fdm(fr.mort) # use model fitted with fdm
Data <- france.fit$fitted$y[1:101,]
N <- dim(Data)[2] # num of curves

# compute the index/time where records are estimated
Idx.records <- IndexRecord(data=Data, depth = 'MBD')

# optional plot
matplot(Data, type = 'l', col = grey(.7,.4))
matplot(Data[,Idx.records$UpperR==1], type = 'l', col = 2, add = TRUE )
matplot(Data[,Idx.records$LowerR==1], type = 'l', col = 4, add = TRUE )
