

UTgold <- 'lightgoldenrod' #saves 'nice' colors
UTblue <- rgb(5/255, 110/255, 167/255, 1)



library(readxl)
CPI2 <- read_excel("CPI2.xlsx")
#####   view(CPI2) 




################ read data##############################################################
library(tseries)
cpi <- ts(CPI2$CPI, start = c(1975, 01), end = c(2005, 12), frequency = 12)

cpi



##calculate inflation#################################################################

inf =  (cpi - lag (cpi, -1))/lag(cpi, -1)*100

inf

### inf = newdata ######### both are the same
library(forecast)  ### autoplot, 

Inflation <- inf

##############plot of inflation #########################################


attach(CPI2)

newdata <- data.frame(CPI2$Date[2:372],inf) #### = inflation

newdata

plot(newdata,type="l")
plot(inf)


inf
newdata




################# normality plot and descriptive statistics ##############################
par(mfrow=c(1,1))

plot(inf, ylab = "Inflation rate", xlab = "year",main = "", col = UTblue, lwd = 1.5, ylim=c(-5,25),xaxt='n',yaxt='n')

xticks_irf <- seq(1975,2005, 5)

axis(1, at = xticks_irf, labels = xticks_irf)

axis(2, at = seq(-5,25,5), labels = seq(-5,25,5))

abline(h=seq(-5,25,5), v=seq(1975,2005,5), col="grey", lty = "solid")

## abline(h=c(-2,0,2,4), v=c(2,4,6,8,10,12,14,16,18,20,22,24), col="grey", lty = "solid")

abline(h=0, col="green", lwd = 2)
legend("topright",legend = c("Inflation rate"), col = c(UTblue), lty = 1:1, lwd = 2.5:2.5, box.col = "black", text.font = 10)



library(forecast)

?autoplot

autoplot(Inflation, col = "black", title = "")

par(mfrow=c(1,1))

plot(Inflation, type= "h")
plot(Inflation, type= "b")

library(ggplot2)


##############################    checking if the data is seasonally adjusted###########

library(seasonal)

sea <- seas(inf)
plot(sea)




summary(inf)


### histogram of inflation ###############################################################
UTgold <- 'lightgoldenrod' #saves 'nice' colors
UTblue <- rgb(5/255, 110/255, 167/255, 1)

hist(inf, breaks=15, freq=FALSE, col=UTgold, 
     xlab= expression(paste('Estimated Error Terms ', widehat(u[t]))), ylab='Density in %', main='', las = 1) 

xfit <- seq(min(inf),max(inf),length=10000) # saves 10,000 x-values to ..
yfit <- dnorm(xfit,mean=mean(inf),sd=sd(inf)) # .. calculate the norm. distr.
lines(xfit, yfit, col="blue", lwd=3) # inserts it in the plot
lines(density(inf), lwd=2) # inserts a non-parametric distr.




## Jarque Bera and Shapiro-Wilk Tests, 3rd and 4th moments
library(moments)
library(tseries)

??jarque.bera.test
jarque.bera.test(inf)
skewness(inf)
kurtosis(inf)


??shapiro.test
library(stats)

shapiro.test(inf)

library(fBasics)
descstats <- basicStats(Inflation)

print(descstats)


######################## QQ plot of inflation ########################
qqnorm(inf)
qqline(inf)


#####   Ljung Box Test ### arch.test() function also included autocorrelation test###################
library(urca)
library(urca)

res <- ur.df(inf,lags=11,type="trend")@testreg$residuals

lag.max <- round(log(length(res))) # m = 7

k <- dim(ur.df(inf, lags = 11, type = "trend")@testreg$coef)[1]

library(FitAR)

x <- LjungBoxTest(res, k=k,lag.max = lag.max, StartLag = 1)
plot(x[,3], main = "Ljung-Box Q Test", ylab = "P-value", xlab = "Lag", ylim= c(0,1))
abline(h=0.05, col = "3")
#It is free of autocorrelation
??LjungBoxTest
LBQPlot(res, lag.max = lag.max, StartLag = 1, k = k, SquaredQ = FALSE)

LjungBoxTest(res)

Box.test(res)

###########################Unit Root Tests ################################

??ur.df

library(urca)

summary(ur.df(na.omit(inf), selectlags = c("AIC"), type='trend'))
summary(ur.df(na.omit(inf), selectlags = c("AIC"), type='drift'))
summary(ur.df(na.omit(inf), selectlags = c("AIC"), type='none'))




??PP.test
library(stats)



PP.test(inf)


#############################################################################

### We conclude that  the null hypothesis is rejected, as statistic value is smaller than its critical value at 5% significant level,
### so we do not have a unit root and now the inflation series are stationary.


### The inflation rate is the order of integration 0

#########################################################################################################################################################################


########## model selection############### PACF breaks up after 2 lags ########## AR(2) ### 

library(FinTS)
library(FitAR)
par(mfrow=c(2,1))

Acf(inf,lag.max = 12)

Pacf(inf)
acf(inf)
par(mfrow=c(1,1))

Acf(inf, lag.max = 12, plot = TRUE, na.action = na.contiguous, demean = TRUE)
Acf(inf, lag.max = 12, type = c("correlation", "covariance", "partial"), plot = TRUE, na.action = na.contiguous, demean = TRUE)

acf2AR(inf)
??acf

library(FinTS)
library(polynom)

## plotArmaTrueacf(inf) #############################################################
library(forecast)

par(mfrow=c(1,1))

library(astsa)

acf2(inf, max.lag=12, plot=TRUE, main="Inflation", ylim=NULL, pacf=FALSE,
     na.action = na.pass)

acf1(inf, max.lag=NULL, plot=TRUE, main="Inflation", ylim=c(-0.1,0.5),xlim = c(0,10), pacf=TRUE,
     na.action = na.pass)

is.acf(inf)

### CPI2 %>%
### plot_acf_diagnostics(Date, CPI, .interactive = interactive)
inf

#### from acf and pacf we MA(2) model is appropriate #################################

model_arma=arima(inf, order=c(2,0,0))

model_arma

########################## for myself to check other stuffs trivial ###############
library(stats)

??phillips
??test

library(urca)

??Ljung

??arch

##### check for autocorrelation and arch effect #############################################################

library(aTSA)

arch.test(model_arma,output = TRUE)



############ GARCH model(1,1)#########################################################
library(rmgarch)
library(rugarch)

ug_spec <-  ugarchspec(mean.model = list(armaOrder=c(2,0)))


ugfit=ugarchfit(spec=ug_spec,data=inf)

ugfit


###########################################################################

library(fGarch)

## Sum of ARCH Coefficients 
sum(ugfit@fit$matcoef[5:6,1])

######################################################################################
names(ugfit@fit)
names(ugfit@model)



### save the estimated conditional variances
ug_var=ugfit@fit$var

ug_var





####save the estimated residuals
ug_res2=ugfit@fit$residuals

####################################################################################
################### residuals squared and conditional variance #####################
ug_res2
ug_var

newdata1 <- data.frame(CPI2$Date[2:372],ug_res2^2)
newdata1


newdata2 <-  data.frame(CPI2$Date[2:372],ug_var)
newdata2



#################################################
newdata
newdata1
newdata2

#########  inflation uncertainty variable - unc - #######################################
newdata2$ug_var

unc <- ts(newdata2$ug_var, start = c(1975,02), end = c(2005,12), frequency = 12)

unc

variance <- unc

ug_res2.2 <- ts(newdata1$ug_res2.2, start = c(1975,02), end = c(2005,12), frequency = 12)
ug_res2.2


plot(ug_res2.2,type='l',ylim=c(-5,300))
lines(unc, col='green')
lines(inf)


########## descriptive statistics of inflation uncertainty and Normality tests

UTgold <- 'lightgoldenrod' #saves 'nice' colors
UTblue <- rgb(5/255, 110/255, 167/255, 1)

library(moments)
skewness(variance)
kurtosis(variance)

library(fBasics)
descstats <- basicStats(variance)
print(descstats)

library(moments)
library(tseries)

??jarque.bera.test
jarque.bera.test(variance)

shapiro.test(variance)

plot(inf, type = "l")
plot(unc, type = "l")


### Unit root test for uncertainty


summary(ur.df(na.omit(unc), selectlags = c("AIC"), type='drift'))



##########################histogram of uncertainty #########

## CairoPNG("inf_infuncertainty.png", width=2400, height=2400, pointsize = 36)

plot(newdata2, type = "l")
par(mfrow=c(1,1))

hist(variance, breaks=15, freq=FALSE, col=UTgold, 
     xlab= expression(paste('Estimated Error Terms ', widehat(u[t]))), ylab='Density in %', main='', las = 1) 

xfit <- seq(min(variance),max(variance),length=10000) # saves 10,000 x-values to ..
yfit <- dnorm(xfit,mean=mean(variance),sd=sd(variance)) # .. calculate the norm. distr.
lines(xfit, yfit, col="blue", lwd=3) # inserts it in the plot
lines(density(variance), lwd=2) # inserts a non-parametric distr.

## dev.off()

par(mfrow=c(1,1))
#############################################################################################

library(plotly)

########################## inflation and inflation uncertainty plots #################################

??plot

par(mfrow=c(2,1))

plot(unc,type='l',tck = 1,col = "red", ylim=c(-5,180),xlab = "year",ylab = "Conditional variance",
                                                  main = "Inflation uncertainty")
plot(inf,type="l",tck = 1,col = UTblue,xlab = "year", main = "Inflation", ylab = "Inflation rate")


par(mfrow=c(1,1))

##################################### VAR MODEL and  ####################################



##finding the optimal lags  ###################################

library(vars)


Var_lag <- cbind(inf,unc)

colnames(Var_lag) <- cbind("Inflation","Inlfation Uncertainty")
lag_select <- VARselect(Var_lag,lag.max = 10,type= "const")
lag_select

## unc <- newdata2$ug_var

model2 <- VAR(data.frame(inf,unc),type = "const", lag.max = 15, ic = "AIC")


summary(model2)

##################################################################################################
??dynlm



library(dynlm)


################################## F test to check the overall impact significance ###############

###### impact of inflation on inflation uncertainty
library(dynlm)



###########################################
mod <- dynlm(unc ~ L(inf, 1:9) + L(unc,1:9))



summary(mod)

## 4b) F-Test for inflation variables: 

summary(mod)$coef[,0]

# Define coeffcients 

library(car)
library(carData)

linearHypothesis(mod,c("L(inf, 1:9)1 +
                        L(inf, 1:9)2 +
                        L(inf, 1:9)3 +
                        L(inf, 1:9)4 +
                        L(inf, 1:9)5 +
                        L(inf, 1:9)6 +
                        L(inf, 1:9)7 +
                        L(inf, 1:9)8 +
                        L(inf, 1:9)9 = 0"))

# inflation variable is jointly significant



###################################################
###### impact of inflation uncertainty on inflation


mod1 <- dynlm(inf ~ L(unc,1:9) + L(inf, 1:9))

summary(mod1)
## 4b) F-Test for inflation uncertainty variables: 

summary(mod1)$coef[,0]
# Define coefficients 

linearHypothesis(mod1,c(
                       "L(unc, 1:9)1  +
                        L(unc, 1:9)2  +
                        L(unc, 1:9)3  +
                        L(unc, 1:9)4  +
                        L(unc, 1:9)5  +
                        L(unc, 1:9)6  +
                        L(unc, 1:9)7  +
                        L(unc, 1:9)8  +
                        L(unc, 1:9)9 = 0"))


############################################################################################################
??linearHypothesis



######################impulse response function  ##############################################################

#######  1st way to plot ################################################

## plot(irf(model2,n.ahead = 25))


######  2nd way to plot ################################################
library(vars)

inflation_infuncertainty <- irf(model2,impulse = "inf", response = "unc"
                                ,boot = TRUE, n.ahead = 25,ortho = TRUE, main = "Shock from Inflation")

infuncertainty_inflation <- irf(model2,impulse = "unc", response = "inf"
                                
                                ,boot = TRUE ,n.ahead = 25,ortho = TRUE, main = "Shock from Inflation Uncertainty")


inflation_inflation <- irf(model2,impulse = "inf", response = "inf"
                           ,boot = TRUE, n.ahead = 25,ortho = TRUE, main = "Shock from Inflation")


infuncertainty_infuncertainty <- irf(model2,impulse = "unc", response = "unc"
                                     ,boot = TRUE ,ortho = TRUE,n.ahead = 25, main = "Shock from Inflation Uncertainty")


par(mfrow=c(1,1))

plot(inflation_infuncertainty, ylab = "Inflation Uncertainty", main = "Shock from Inflation")

plot(inflation_inflation, ylab = "Inflation", main = "Shock from Inflation")

plot(infuncertainty_inflation, ylab = "Inflation", main = "Shock from Inflation Uncertainty")

plot(infuncertainty_infuncertainty, ylab = "Inflation Uncertainty", main = "Shock from Inflation Uncertainty")


library(lpirfs)
library(panelvar)
??girf



######## granger causality to test whether inflation cause ug_var (inflation uncertainity) with order 1

library(lmtest)

grangertest(unc~inf, order = 9)

################# granger causality to test ug_var cause inf with order 1

grangertest(inf~unc, order = 9)



####   Variance Decomposition ##################################################

vd <- fevd(model2, n.ahead = 24)
plot(vd)

??fevd

########################  VAR DIagnostics ######################################

######################  Serial Correlation #####################################
library(vars)

Serial <- serial.test(model2,lags.pt = 12, type = "PT.asymptotic")
Serial

library(sandwich)



##################### Heteroskedasticity Test ##################################

archtest2 <- arch.test(model2,lags.multi = 12,multivariate.only = TRUE)
archtest2


###################  normality test ############################################

normality <- normality.test(model2)
normality


#####   testing for structural breaks in the residuals #########################
### please close R, open and run this code first again. 

library(vars)
stability <- stability(model2, type = "OLS-CUSUM")
plot(stability)

?stability

##### there is no structural breaks
################################################################################




##############################################################################################

##  ###  3rd way to plot IRFs

irf_inf_infuncertainty_95  = irf(VAR(data.frame(inf,unc), p=4), impulse="inf", response="unc",     
                                 n.ahead = 24, runs=500, ci=0.95)


irf_infuncertainty_inf_95   = irf(VAR(data.frame(inf,unc), p=4), impulse="unc", response="inf",     
                                  n.ahead = 24, runs=500, ci=0.95)

irf_inf_inf_95 = irf(VAR(data.frame(inf,unc), p=4), impulse="inf", response="inf",     
                     n.ahead = 24, runs=500, ci=0.95)

irf_infuncertainty_infuncertainty_95   = irf(VAR(data.frame(inf,unc), p=4), impulse="unc", response="unc",     
                                             n.ahead = 24, runs=500, ci=0.95)

??irf


plot(irf_inf_infuncertainty_95, ylab = "Inflation Uncertainty", main = "Shock from Inflation")

plot(irf_infuncertainty_inf_95, ylab = "Inflation", main = "Shock from Inflation Uncertainty")

plot(irf_inf_inf_95, ylab = "Inflation", main = "Shock from Inflation")

plot(irf_infuncertainty_infuncertainty_95, ylab = "Inflation Uncertainty", main = "Shock from Inflation Uncertainty")

##############################################################################################
#################################################

## Define xticks_irf
## E.g., 1,7,13,19,25

library(Cairo)

library(graphics)

par(mfrow=c(1,1))

################# 4th way to plot ############################

######### IRFs - shock from inflation uncertainty on inflation


CairoPNG("response of inflation to uncertainty.png", width=1600, height=1000, pointsize = 36)

xticks_irf <- seq(1,25, 2)

plot(irf_infuncertainty_inf_95[[1]]$unc, ylab="Inflation", xlab="", main="Shock from inflation uncertainty", type="l", lwd=2, ylim=c(-0.7,0.4), xaxt='n',yaxt='n')


polygon(c(1:25, rev(1:25)),c(irf_infuncertainty_inf_95[[2]]$unc,rev(irf_infuncertainty_inf_95[[3]]$unc)), col="lightgrey", border=NA)

axis(1, at = xticks_irf, labels = xticks_irf-1)

axis(2, at = seq(-0.4,0.2,0.2), labels = seq(-0.4,0.2,0.2))

lines(irf_infuncertainty_inf_95[[1]]$unc, col="blue", lwd = 2, pch=2)

abline(h=c(-2,0,2,4), v=c(1,3,5,7,9,11,13,15,17,19,21,23,25), col="grey", lty = "solid")

## abline(h=c(-2,0,2,4), v=c(2,4,6,8,10,12,14,16,18,20,22,24), col="grey", lty = "solid")

abline(h=0, col="red")
dev.off()

######### IRFs - shock from inflation uncertainty on inflation uncertainty


CairoPNG("response of uncertainty to uncertainty.png", width=1600, height=1000, pointsize = 36)

xticks_irf <- seq(1,25, 2)

plot(irf_infuncertainty_infuncertainty_95[[1]]$unc, ylab="Inflation Uncertainty", xlab="", main="Shock from inflation uncertainty", type="l", lwd=2, ylim=c(-1,12), xaxt='n',yaxt='n')


polygon(c(1:25, rev(1:25)),c(irf_infuncertainty_infuncertainty_95[[2]]$unc,rev(irf_infuncertainty_infuncertainty_95[[3]]$unc)), col="lightgrey", border=NA)

axis(1, at = xticks_irf, labels = xticks_irf-1)

axis(2, at = seq(-10,12,2), labels = seq(-10,12,2))

lines(irf_infuncertainty_infuncertainty_95[[1]]$unc, col="blue", lwd = 2, pch=2)

abline(h=c(0), v=c(1,3,5,7,9,11,13,15,17,19,21,23,25), col="grey", lty = "solid")

## abline(h=c(-2,0,2,4), v=c(2,4,6,8,10,12,14,16,18,20,22,24), col="grey", lty = "solid")

abline(h=0, col="red")
dev.off()

######## IRFs - shock from inflation  on inflation 

CairoPNG("response of inflation to inflation.png", width=1600, height=1000, pointsize = 36)

xticks_irf <- seq(1,25, 2)

plot(irf_inf_inf_95[[1]]$inf, ylab="Inflation", xlab="", main="Shock from inflation", type="l", lwd=2, ylim=c(0,2.1), xaxt='n',yaxt='n')


polygon(c(1:25, rev(1:25)),c(irf_inf_inf_95[[2]]$inf,rev(irf_inf_inf_95[[3]]$inf)), col="lightgrey", border=NA)

axis(1, at = xticks_irf, labels = xticks_irf-1)

axis(2, at = seq(0,2,0.5), labels = seq(0,2,0.5))

lines(irf_inf_inf_95[[1]]$inf, col="blue", lwd = 2, pch=2)

abline(h=c(-2,0,2,4), v=c(1,3,5,7,9,11,13,15,17,19,21,23,25), col="grey", lty = "solid")

##abline(h=c(-2,0,2,4), v=c(2,4,6,8,10,12,14,16,18,20,22,24), col="grey", lty = "solid")

abline(h=0, col="red")

dev.off()


######## IRFs - shock from inflation  on inflation uncertainty

CairoPNG("response of uncertainty to inflation.png", width=1600, height=1000, pointsize = 36)

xticks_irf <- seq(1,25, 2)

plot(irf_inf_infuncertainty_95[[1]]$inf, ylab="Inflation Uncertainty", xlab="", main="Shock from inflation", type="l", lwd=2, ylim=c(-3,9), xaxt='n',yaxt='n')


polygon(c(1:25, rev(1:25)),c(irf_inf_infuncertainty_95[[2]]$inf,rev(irf_inf_infuncertainty_95[[3]]$inf)), col="lightgrey", border=NA)

axis(1, at = xticks_irf, labels = xticks_irf-1)

axis(2, at = seq(-2,8,2), labels = seq(-2,8,2))

lines(irf_inf_infuncertainty_95[[1]]$inf, col="blue", lwd = 2, pch=2)

abline(h=c(0), v=c(1,3,5,7,9,11,13,15,17,19,21,23,25), col="grey", lty = "solid")

## abline(h=c(-2,0,2,4), v=c(2,4,6,8,10,12,14,16,18,20,22,24), col="grey", lty = "solid")

abline(h=0, col="red")

dev.off()


