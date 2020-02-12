load("C:/Users/katca/Desktop/Reproducible_Science/NLM_Workshop.RData")
library(nlstools)
par(mai=c(1,1,0.1,0.1))
plot(harv$TIMESTAMP, harv$NEE, 
     ylab=expression(paste("NEE (",mu,"mol m"^{-2} ~ s^{-1} ~ ")" )), xlab="")
#creating a plot of NEE. harv is whole dataset, dollar sign seperates, timestamp is that data column you want to work with

plot( NEE ~ PAR, data= day)
y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, data=day[which(day$MONTH == 07),], 
        start=list(a1= -1 , ax= -1, r= 1), 
        na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
summary(y)
lrcModel <- function(PAR, a1, ax, r) { 
  NEE <- (a1 * PAR * ax)/(a1 * PAR + ax) + r 
  return(NEE) 
  }
lrc.int <- function (mCall, LHS, data){
  x <- data$PAR
  y <- data$NEE
  r <- max(na.omit(y), na.rm=T) #maximum NEE
  ax <- min(na.omit (y), na.rm=T) #minimum NEE
  a1 <- (r+ax)/2 #midway between r and a1
  
  a1[a1 > 0] <- -0.1
  r[r > 50] <- ax* -1
  r[r < 0] <- 1
  
  value = list(a1, ax, r) #Must be included for selfStart function
  names(value) <-mCall[c("a1", "ax", "r")] #must be included for selfStart funct
  return(value)
}
# random note: telling what you want the program to do as you would want to present the data visually
SS.lrc <- selfStart(model=lrcModel,initial= lrc.int)
iv <- getInitial(NEE ~ SS.lrc('PAR', "a1", "ax", "r"),
                 data = day[which(day$MONTH == 07),])
iv
y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, day[which(day$MONTH == 07),],
         start=list(a1= iv$a1 , ax= iv$ax, r= iv$r),
         na.action=na.exclude, trace=F, control=nls.control (warnOnly=T))
summary(y)
res.lrc <- nlsResiduals(y)
par(mfrow=c(2,2))
plot(res.lrc, which=1)
plot(res.lrc, which=3)
plot(res.lrc, which=4)
plot(res.lrc, which=5)
results <- nlsBoot(y, niter=100 )
summary(results)
plot(results, type = "boxplot")




#Exercise: How variable are NEE rates over an annual cycle in Harvard Forest?


parms.Month <- data.frame( 
  MONTH=numeric(), 
  a1=numeric(), 
  ax=numeric(), 
  r=numeric(), 
  a1.pvalue=numeric(), 
  ax.pvalue=numeric(), 
  r.pvalue=numeric(), stringsAsFactors=FALSE, row.names=NULL)

parms.Month[1:12, 1] <- seq(1,12,1) # Adds months to the file
nee.day <- function(dataframe) { y = nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, dataframe, 
                                         start=list(a1= iv$a1 , ax= iv$ax, r= iv$r),
                                         na.action=na.exclude, trace=F, 
                                         control=nls.control(warnOnly=T))

y.df <- as.data.frame(cbind(t(coef(summary(y)) [1:3, 1]), t(coef(summary(y)) [1:3, 4])))
names(y.df) <-c ("a1","ax", "r", "a1.pvalue", "ax.pvalue", "r.pvalue")
return (y.df )
}

try(for(j in unique(day$MONTH)){
  iv <- getInitial(NEE ~ SS.lrc('PAR', "a1", "ax", "r"), data = day[which(day$MONTH == j),])
  y3 <- try(nee.day(day[which(day$MONTH == j),]), silent=T)
  try(parms.Month[c(parms.Month$MONTH == j ), 2:7 ] <- cbind(y3), silent=T)
  rm(y3) }, silent=T)
parms.Month

boot.NEE <- data.frame(parms.Month[, c("MONTH")]); names (boot.NEE) <- "MONTH" 
boot.NEE$a1.est <- 0 
boot.NEE$ax.est<- 0 
boot.NEE$r.est<- 0 
boot.NEE$a1.se<- 0 
boot.NEE$ax.se<- 0 
boot.NEE$r.se<- 0

for ( j in unique(boot.NEE$Month)){
y1 <-day[which(day$MONTH == j),]
iv <- getInitial( NEE - SS.lrc('PAR', "a1", "ax", "r"), data=y1)
day.fit <- nls( NEE ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r, data=y1, 
                start=list(a1= iv$a1 , ax= iv$ax, r= iv$r), 
                na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
try(results <- nlsBoot(day.fit, niter=100 ), silent=T)
try(a <- t(results$estiboot)[1, 1:3], silent=T)
try(names(a) <- c('a1.est', 'ax.est', 'r.est'), silent=T)
try( b <- t(results$estiboot)[2, 1:3], silent=T)
try(names(b) <- c('a1.se', 'ax.se', 'r.se'), silent=T)
try(c <- t(data.frame(c(a,b))), silent=T)
try(boot.NEE[c(boot.NEE$MONTH == j), 2:7] <- c[1, 1:6], silent=T)
try(rm(day.fit, a, b, c, results, y1), silent=T)
}

lrc <-merge( parms.Month, boot.NEE, by.x="MONTH", by.y="MONTH")
lrc












#================================================================================================================
#Challenge  Fit monthly temperature respponse curves using a similar approach with the night data from harv(night)

load("C:/Users/katca/Desktop/Reproducible_Science/NLM_Workshop.RData")
library(nlstools)
rm(list=ls())

trcModel <- function(TA, a, b) {
  y=a * exp(b*TA)
  return(y)
}

trc.int <- function (mCall, LHS, data){
  x <- data$TA
  y <- data$NEE
  
  a <-1.00703982 + -0.08089044* (min(na.omit(y)))
  b <- 0.051654 + 0.001400 * (min(na.omit(y))) 
  
  value = list(a, b)
  names(value) <- mCall[c("a", "b")]
  return(value)
}



SS.trc <- selfStart(model=trcModel,initial= trc.int)

trc <- getInitial(y ~ SS.trc('TA', "a", "b"), 
                 data = night[which(night$MONTH == 12),]) 
trc 

y = nls( NEE ~ (a * exp(b * TA)), night[which(night$MONTH == 12),], 
         start=list(a= trc$a , b= trc$b), 
         na.action=na.exclude, trace=F, control=nls.control(warnOnly=T))
summary(y)

res.trc <- nlsResiduals(y) 
par(mfrow=c(2,2)) 
plot(res.trc, which=1)# Residulas vs fitted values (Constant Variance) 
plot(res.trc, which=3) # Standardized residuals 
plot(res.trc, which=4) # Autocorrelation 
plot(res.trc, which=5) # Histogram (Normality)




parms.Month <- data.frame(
  MONTH=numeric(),
  a=numeric(),
  b=numeric(), 
  a.pvalue=numeric(),
  b.pvalue=numeric(), stringsAsFactors=FALSE, row.names=NULL)
parms.Month[1:12, 1] <- seq(1,12,1)


nee.night <- function(dataframe){y.df = nls(NEE ~ a * exp(b*TA), 
                                            dataframe, start=list(a=iv$a , b=iv$b ),
                                            na.action=na.exclude, trace=F,
                                            control=nls.control(warnOnly=T))

y.df <- as.data.frame(cbind(t(coef(summary(y.df))[1:2, 1]), t(coef(summary(y.df)) [1:2, 4])))

names(y.df) <- c("a", "b", "a.pvalue", "b.pvalue")                      
return(y.df)}

try(for(j in unique(night$MONTH)){
  print(j)
  
  trc <- getInitial(NEE ~ SS.trc('TA', "a", "b"), data = night[which(night$MONTH == j),]) 
  
  y4 <- try(nee.night(night[which(night$MONTH == j),]), silent=T) # Fit night model
  
  try(parms.Month[c(parms.Month$MONTH == j ), 2:5 ] <- cbind(y4), silent=T)
  
  rm(y4)
}, silent=T)


boot.NEE <- data.frame(parms.Month[, c("MONTH")]); names (boot.NEE) <- "MONTH"
boot.NEE$a.est<- 0
boot.NEE$b.est<- 0
boot.NEE$a.se<- 0
boot.NEE$b.se<- 0


for ( j in unique(boot.NEE$MONTH)){
  print(j)
  y1 <-night[which(night$MONTH == j),]
  
  iv <- getInitial(NEE ~ SS.trc('TA',"a", "b"), data = y1) 
  
  night.fit <- nls(NEE ~ a * exp(b*TA), 
                   data=y1, start=list(a= iv$a , b=iv$b ),
                   na.action=na.exclude, trace=F,
                   control=nls.control(warnOnly=T))
  
  results <- nlsBoot(night.fit, niter=100 )
  a <- t(results$estiboot)[1, 1:2]
  names(a) <- c('a.est', 'b.est')
  b <- t(results$estiboot)[2, 1:2]
  names(b) <- c('a.se', 'b.se')
  c <- t(data.frame(c(a,b)))
  boot.NEE[c(boot.NEE$MONTH == j), 2:5] <- c[1, 1:4]
  rm(night.fit, a, b, c, results, y1)
}

trc <- merge( parms.Month, boot.NEE)
trc

results <- nlsBoot (y4, niter=100 ) #bootstrap, way of subsampling data set, number of times it doesnt run well is part of the error estimation, we had it run 100 times. How well model will predict results that our data.
summary(results)
plot(results, type = "boxplot")



summary(y)
# show images of month, how they are similar. include table, and image for a month, august worst month bc more variability, compare with other month that has less variablity as the plot wont be so scattered



warning()
