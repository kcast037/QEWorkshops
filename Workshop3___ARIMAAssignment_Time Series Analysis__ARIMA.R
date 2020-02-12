load("C:/Users/katca/Desktop/Reproducible_Science/ARIMA_Workshop.RData")
library(zoo)
library(tseries)
library(forecast)
library(xts)
nee <- ts( mangroves$nee, start= 1, frequency=30)
par(mfrow=c(1,1), mai=c(0.25,0.8,0.1, 0.1)) 
plot( nee, typ="l", ylab= "NEE",
xlab="")
plot(nee)
lines(tsclean(nee),col="red")
nee <-tsclean(nee)

nee.d <- decompose(nee, 'multiplicative')
plot(nee.d)

adf.test(nee)

acf(nee, lag.max=45)
pacf(nee, lag.max=45)

arima.nee1 <-auto.arima(nee, trace=TRUE)

tsdisplay(residuals(arima.nee1), lag.max=45)

arima.nee2 <-arima(nee , order=c(10,1,3), seasonal= list(order=c(2,0,2))) 

tsdisplay(residuals(arima.nee2), lag.max= 30)

AIC(arima.nee1, arima.nee2) 

par(mfrow=c(1,1))

plot(nee , typ="l"); lines(fitted(arima.nee2),col="red")
checkresiduals(arima.nee2, lag=36)
par(mfrow=c(1,1)) 
plot(nee , typ="l"); lines(fitted(arima.nee2),col="red")
plot(forecast(arima.nee2, h=30))



sal <- ts(mangroves$salinity.max, start= 1, frequency=30)
par(mfrow=c(1,1), mai=c(0.25,0.8,0.1, 0.1)) 

plot(sal , typ="l", ylab= "Salinity", xlab="")

plot(sal , typ="l", ylab= "Salinity", xlab="") 
lines(tsclean(sal) , col="red")

sal <- tsclean(sal)

sal.d <- decompose(sal, 'multiplicative') 
plot(sal.d)

adf.test(sal) 

adf.test(diff(sal)) 

ccf( diff(sal),nee, na.action = na.pass, lag.max=40, plot=TRUE)
arima.nee3 <-auto.arima(nee, xreg=c(diff(sal),0), trace=TRUE)

AIC(arima.nee2, arima.nee3 ) 
sal.i <- sal 
sal.i[sal.i < 25 ]<- 0 
sal.i[sal.i >= 25 ]<- 1
plot(sal.i)

arima.nee4 <-auto.arima(nee, xreg=sal.i, trace=TRUE)
AIC(arima.nee2,arima.nee4 ) 
checkresiduals(arima.nee4, lag=36)

par(mfrow=c(1,1)) 
plot(nee , typ="l"); lines(fitted(arima.nee4),col="red")

#________________________________________________________________

# Creating better model begins here
# 1. Create a Timeseries object

everair <-ts(mangroves$tair, start= 1, frequency=30)

# visualize data
par(mfrow=c(1,1), mai=c(0.25, 0.8, 0.1, 0.1))
plot(everair , typ="l", ylab= "tair", xlab="")

#       removing any outliers that could bias the model by skewing statistical summaries
plot(everair , typ="l", ylab= "tair", xlab="") 
lines(tsclean(everair) , col="orange")
everair <- tsclean(everair)

# 2. Decompose the time series

everair.d <- decompose(everair, 'multiplicative')
plot(everair.d)

# 3a. Test for Stationarity 

adf.test(everair) 
adf.test(diff(everair)) 

#4 Explore correlations

ccf( diff(everair), nee, na.action =na.pass, lag.max=40, plot=TRUE)

# 5 Explore Models of NEE

arima.nee3 <-auto.arima(nee, xreg=c(diff(everair),0), trace=TRUE)

AIC(arima.nee2, arima.nee3 )



checkresiduals(arima.nee3, lag=36)


par(mfrow=c(1,1)) 
plot(nee , typ="l"); lines(fitted(arima.nee3),col="red")
#AIC 660 arima.nee3

