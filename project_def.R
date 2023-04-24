rm(list = ls())
graphics.off()

library(readxl)
library(tseries)
library(rugarch)
library(readr)
library(fDMA)
library(fGarch)
library(rumidas)
library(xts)
library(zoo)
library(FinTS)
library(MCS)
library(segMGarch)
library(pracma)
library(moments)
library(ggplot2)
library(lmtest)

Energy_data_ <- read_excel("Energy data .xls")

Data <- na.omit(Energy_data_)
Data1 <- Data[, c(1,2)]
sum(is.na(Data1))
Data2 <- Data1[-36, ]

dataset <- Data2

rm(Data, Data1, Data2, Energy_data_)

### WTI ###

WTI <- dataset
WTI_ts <- ts(WTI, start = c(2020,2), frequency = 213)
ts.WTI <- WTI_ts[ ,2]
plot(ts.WTI, ylab = "WTI")
acf(ts.WTI, main = "WTI") #no stationary
pacf(ts.WTI)

logWTI_ts <- log(ts.WTI)
diffWTI_ts <- diff(logWTI_ts)

plot(diffWTI_ts)
M <- mean(diffWTI_ts)
abline(h = M, col = "red")
acf(diffWTI_ts) #stationary
pacf(diffWTI_ts)

pdf("WTI.pdf")
plot(ts.WTI, ylab = "WTI")
acf(ts.WTI, main = "") 
pacf(ts.WTI, main = "")
acf(diffWTI_ts, main = "")
pacf(diffWTI_ts, main = "")
dev.off()

# DICKEY-FULLER TEST
# Here below, we are verifying if our TS is stationary or not through the Dickey-Fuller test: we adopt this test two times;
# in the first case, we will test the first time series and then we will test the log-diff1 TS
# in the last case, we expect to reject the null Ho due to the stationarity

adf.test(ts.WTI, alternative ="stationary")
# We cannot reject the null hypothesis because the process is non stationary
adf.test(diffWTI_ts, alternative = "stationary") 
# Computing the first differences, we have a strong relevance that the process is stationary (p-values < 0.05)

logWTI_ts <- log(ts.WTI) #logprices
diffWTI_ts <- diff(logWTI_ts) #difference of order 1 (log) = the daily returns series
plot.ts(diffWTI_ts) #stationarity and daily returns
retWTI <- diffWTI_ts
sq_retWTI <- retWTI^2
plot(retWTI, main = "WTI returns", ylab = "Returns")
plot(sq_retWTI, main = "WTI squared returns", ylab = "Returns")

# ACF and PACF
par(mfrow = c(4,2))
acf(retWTI)
acf(retWTI, lag.max = 70)
pacf(retWTI)
pacf(retWTI, lag.max = 70)

acf(sq_retWTI)
acf(sq_retWTI, lag.max = 70)
pacf(sq_retWTI)
pacf(sq_retWTI, lag.max = 70)

par(mfrow = c(1,1))

#QQ - PLOT
qqnorm(y = retWTI, main = "Normal QQ plot WTI", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(retWTI, distribution = qnorm, col = "red")

# STATISTICAL ANALYSIS (CI = 95%)
basicStats(retWTI) # negative skewness and leptokurtic distribution. Our TS does not have ~N distribution
# This means that high and low returns in absolute value are more frequent than those expected for a Gaussian distribution
# In statistical aspect, it means that high return in absolute value represent huge gain or losses

# KURTOSIS TEST
kurtosis(sq_retWTI) #coherently with the previous result, we note that the values exceed positively 0 and confirm the leptokurtic distribution

# BOX-PIERCE TEST Ho: no correlation
Box.test(retWTI, lag = 35, type = "Box-Pierce") #We need to reject the null. For this reason we can conclude that within the 35th lag there's correlation
Box.test(retWTI, lag = 100, type = "Box-Pierce") #We cannot reject the null because crossing the 100th lag, there's no correlation among residuals/returns

Box.test(sq_retWTI, lag = 35, type = "Box-Pierce") #We reject the null
Box.test(sq_retWTI, lag = 100, type = "Box-Pierce") #After the 100th lag, there's no correlation 

# SHAPIRO-WILK TEST (We can also test the normality of residuals/returns: HO = normally distributed)
shapiro.test(retWTI) #As we said before, the data are not normally distributed, therefore we reject the null

# HISTOGRAMS
pdf("histWTI.pdf")

hist(retWTI, xlab = "Daily stock prices", prob = T, main = "WTI daily returns", freq = F, breaks = 20) #shifted distribution different from the Normal one

xfit <- seq(min(retWTI), max(retWTI), length = 1000)
yfit <- dnorm(xfit, mean = M, sd = sd(retWTI)) #normal distribution
lines(xfit, yfit, type = "l", col = "blue", lwd = 1.5) 

y1fit <- dstd(xfit, mean = M, sd = sd(retWTI)) #student-t distribution
lines(xfit, y1fit, type = "l", col = "red", lwd = 1.5) 

y2fit <- dged(xfit, mean = M, sd = sd(retWTI)) #generalized error distribution
lines(xfit, y2fit, type = "l", col = "green", lwd = 1.5)

y3fit <- dsstd(xfit, mean = M, sd = sd(retWTI)) #skew student-t distribution
lines(xfit, y3fit, type = "l", col = "black", lwd = 1.5)

y4fit <- dsnorm(xfit, mean = M, sd = sd(retWTI)) #skew normal distribution
lines(xfit, y4fit, type = "l", col = "purple", lwd = 1.5)

legend("topright", legend = c("std", "ged", "sstd", "snorm"), col = c("red", "green", "black", "purple"), 
       lty = 1, cex = 0.9)
dev.off()

# FITTING MODELS
models <- c("sGARCH", "eGARCH")
distributions <- c("norm", "std", "snorm", "sstd", "ged")
m.d_bind <- list()
m <- c()
d <- c()
for (m in models) {
  for (d in distributions) {
    m.d_bind[[paste(m, d, sep = "-" )]] <-  
      ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                 variance.model = list(model = m, garchOrder = c(1, 1)),
                 distribution.model = d)
  }
}

spec.model <- names(m.d_bind); spec.model

fitmod <- list()
for (s in spec.model) {
  fitmod[[s]] <- ugarchfit(spec = m.d_bind[[s]], data = retWTI)
}
fitmod

model_1 <- fitmod$`eGARCH-sstd`
plot(model_1, which = "all")
plot(model_1, which = 2) 

res1 <- model_1@fit$residuals
res1

plot(res1)
abline(h=0, col = "red")

acf(res1)
acf(res1, lag.max = 50) # after lag 29, the residuals seems uncorrelated
acf(res1, lag.max = 100) 
pacf(res1, lag.max = 50)
pacf(res1, lag.max = 100)

jarque.bera.test(res1) # Ho: normality
Box.test(res1, lag = 30, type = "Ljung-Box") # Ho: no correlations

dates <- as.Date(WTI$...1)
dates_1 <- dates[-1]
dates <- dates_1 
rm(dates_1)

retWTI_xts <- xts(retWTI, order.by = dates)

# PERFORMING VaR (GARCH models)
VaR_test <- list()
for (s in spec.model) {
  VaR_test[[s]] <- ugarchroll(spec = m.d_bind[[s]], data = retWTI_xts, forecast.length = 150, refit.every = 5,
                            refit.window = "moving", calculate.VaR = T)
}

VaR_test

plot(VaR_test$`eGARCH-sstd`, which = 4, VaR.alpha = 0.05)
plot(VaR_test$`eGARCH-sstd`, which = 4, VaR.alpha = 0.01)

# BACKTESTING PROCEDURE GARCH (VaR)

#Report gives information about the VaR backtest with VaR 5%
BT95 <- list()
for (s in spec.model) {
  BT95[[s]] <- as.data.frame(report(VaR_test[[s]], type = "VaR", VaR.alpha = 0.05, conf.level = 0.95))
}

#Report gives information about the VaR backtest with VaR 1%
BT99 <- list()
for (s in spec.model) {
BT99[[s]] <- as.data.frame(report(VaR_test[[s]], type = "VaR", VaR.alpha = 0.01, conf.level = 0.99))
}

# DQ test at confidence level 95%
DQ5 <- list()
for (s in spec.model) {
  DQ5[[s]] <- DQtest(retWTI, VaR = VaR_test[[s]]@forecast$VaR$`alpha(5%)`, VaR_level = 0.95)
}

# DQ test at confidence level 99%
DQ <- list()
for (s in spec.model) {
  DQ[[s]] <- DQtest(retWTI, VaR = VaR_test[[s]]@forecast$VaR$`alpha(1%)`, VaR_level = 0.99)
}

# MCS PROCEDURE
# Loss function with VaR confidence level 95%

Loss5 <- do.call(cbind, lapply(spec.model, function(s) { 
  LossVaR(tau = 0.05, realized = VaR_test[[s]]@forecast[["VaR"]][["realized"]],
          evaluated = VaR_test[[s]]@forecast[["VaR"]][["alpha(5%)"]])
}))

colnames(Loss5) <- spec.model

# Superior Set of Models 95%

SSM5 <- MCSprocedure(Loss = Loss5, alpha = 0.2, B = 5000, statistic = "Tmax")
SSM5

# Loss function with VaR confidence level 99%

Loss <- do.call(cbind, lapply(spec.model, function(s) { 
       LossVaR(tau = 0.01, realized = VaR_test[[s]]@forecast[["VaR"]][["realized"]],
                                evaluated = VaR_test[[s]]@forecast[["VaR"]][["alpha(1%)"]])
}))

colnames(Loss) <- spec.model                        

# Superior Set of Models 99%

SSM <- MCSprocedure(Loss = Loss, alpha = 0.2, B = 5000, statistic = "Tmax")
SSM

#################################################################
#################################################################

# Now, we add the weekly data of Covid deaths in US
data_Covid <- read_excel("Covid deaths.xlsx")

dates1 <- as.Date(data_Covid$...1)
Covid_ts <- ts(data_Covid$`Deaths USA`, start = c(2020,3), frequency = 44)
plot(Covid_ts)
Covid_xts <- xts(Covid_ts, order.by = dates1)
diffCovid.ts <-diff(Covid_ts)
acf(diffCovid.ts)
adf.test(diffCovid.ts, alternative = "stationary") 

Date <- as.Date(dates)
retWTI_TS <- ts(retWTI)
WTI_weekly_sum <- apply.weekly(retWTI_TS, sum)
WTI_weekly <- as.xts(coredata(WTI_weekly_sum), order.by = Date, by = "week", length.out = length(WTI_weekly_sum))

Deaths <- ts(data_Covid$`Deaths USA`) 
Date1 <- as.Date(data_Covid$...1)
death_weekly_sum <- apply.weekly(Deaths, sum)
death_weekly <- as.xts(coredata(death_weekly_sum), order.by = Date1, by = "week", length.out = length(death_weekly_sum))

mat <- mv_into_mat(WTI_weekly["2020-04-06/"], diff(death_weekly), K = 4, type = "weekly")
mat

# ESTIMATION STRATEGY FOR GARCH-MIDAS MODELs
GM_model <- list()
distr1 <- c("norm","std")
skewn <- c("YES", "NO")
di <- c()
sk <- c()
for (di in distr1) {
  for(sk in skewn) {
    GM_model[[paste(di, sk,  sep = "-")]] <-
        ugmfit(model = "GM", skew = sk, distribution = di,
           daily_ret = retWTI_xts["2020-04-06/"], mv_m = mat,
           K = 4, out_of_sample = 150)
  }
}

summary.rumidas(GM_model$`norm-YES`)
GM_model$`norm-YES`$inf_criteria

summary.rumidas(GM_model$`std-YES`)
GM_model$`std-YES`$inf_criteria

summary.rumidas(GM_model$`norm-NO`)
GM_model$`norm-NO`$inf_criteria

summary.rumidas(GM_model$`std-NO`)
GM_model$`std-NO`$inf_criteria

# VaR
var_snorm <- qnorm(0.05) * GM_model$`norm-YES`$est_vol_oos
var_sstd <- qstd(0.05) * GM_model$`std-YES`$est_vol_oos
var_norm_no <- qnorm(0.05) * GM_model$`norm-NO`$est_vol_oos
var_std_no <- qstd(0.05) * GM_model$`std-NO`$est_vol_oos

var_snorm99 <- qnorm(0.01) * GM_model$`norm-YES`$est_vol_oos
var_sstd99 <- qstd(0.01) * GM_model$`std-YES`$est_vol_oos
var_norm_no99 <- qnorm(0.01) * GM_model$`norm-NO`$est_vol_oos
var_std_no99 <- qstd(0.01) * GM_model$`std-NO`$est_vol_oos

# Backtesting VaR
VaRTest(alpha = 0.05, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_snorm), conf.level = 0.95)
VaRTest(alpha = 0.05, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_sstd), conf.level = 0.95)
VaRTest(alpha = 0.05, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_norm_no), conf.level = 0.95)
VaRTest(alpha = 0.05, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_std_no), conf.level = 0.95)
VaRTest(alpha = 0.01, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_snorm99), conf.level = 0.99)
VaRTest(alpha = 0.01, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_sstd99), conf.level = 0.99)
VaRTest(alpha = 0.01, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_norm_no99), conf.level = 0.99)
VaRTest(alpha = 0.01, actual = as.numeric(retWTI_xts[166:315]), VaR = as.numeric(var_std_no99), conf.level = 0.99)

# DQ test
DQtest(retWTI_xts, as.numeric(var_snorm), VaR_level = 0.95)
DQtest(retWTI_xts, as.numeric(var_sstd), VaR_level = 0.95)
DQtest(retWTI_xts, as.numeric(var_norm_no), VaR_level = 0.95)
DQtest(retWTI_xts, as.numeric(var_std_no), VaR_level = 0.95)
DQtest(retWTI_xts, as.numeric(var_snorm99), VaR_level = 0.99)
DQtest(retWTI_xts, as.numeric(var_sstd99), VaR_level = 0.99)
DQtest(retWTI_xts, as.numeric(var_norm_no99), VaR_level = 0.99)
DQtest(retWTI_xts, as.numeric(var_std_no99), VaR_level = 0.99)

specifications5 <- c("midas_snorm", "midas_sstd", "midas_norm_no", "midas_std_no")
var.comp5 <- list(midas_snorm = var_snorm, midas_sstd = var_sstd, midas_norm_no = var_norm_no, midas_std_no = var_std_no)

# Loss function with VaR confidence level 95%
midas_Loss5 <- do.call(cbind, lapply(specifications5,
                                     function(s) LossVaR(tau = 0.05, realized = retWTI_xts[166:315], evaluated = var.comp5[[s]])))
colnames(midas_Loss5) <- specifications5

# Superior Set of Models 95%
SSM_midas5 <- MCSprocedure(Loss = midas_Loss5, alpha = 0.2, B = 5000, statistic = "Tmax")

# Loss function with VaR confidence level 99%   
specifications <- c("midas_snorm", "midas_norm_no", "midas_std_no") # creiamo una lista per eliminare manualmente dall'mcs i modelli che non hanno superato 2/3 test del Backtesting
midas_succ.test <- list(midas_snorm = var_snorm99, midas_norm_no = var_norm_no99, midas_std_no = var_std_no99)

midas_Loss <- do.call(cbind, lapply(specifications,
                                     function(s) LossVaR(tau = 0.01, realized = retWTI_xts[166:315], evaluated = midas_succ.test[[s]])))
colnames(midas_Loss) <- specifications

# Superior Set of Models 99%
SSM_midas <- MCSprocedure(Loss = midas_Loss, alpha = 0.2, B = 5000, statistic = "Tmax")

x <- as.Date(dataset$...1)

bestmodels <- qplot(y = VaR_test[["eGARCH-snorm"]]@forecast[["VaR"]][["alpha(1%)"]], x = x[166:315], geom = 'line', color = 'blue') +
  geom_point(aes(x = x[166:315] , y = VaR_test[["eGARCH-snorm"]]@forecast[["VaR"]][["realized"]], color = "gray")) + 
  geom_line(aes(y= VaR_test[["eGARCH-norm"]]@forecast[["VaR"]][["alpha(5%)"]], x = x[166:315], color='black')) +
  geom_line(aes(y = var_sstd, x = x[166:315], color = "red")) +
  geom_line(aes(y = var_std_no99, x = x[166:315], color = "535"))+
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  scale_color_manual(values = c("black", "blue", "red", "grey", 535), labels = c("eGARCH-norm VaR 5%", "eGARCH-snorm VaR 1%", "GM-sstd VaR 5%", "Realized", "GM-std VaR 1%")) +
  labs(y = 'Daily Returns' , x = '', title = 'WTI: comparison among all the best models') + theme_light() + 
  theme(legend.position = "bottom",
        legend.title = element_text(size = 3), 
        legend.text = element_text(size = 7)) 
bestmodels

Losses95 <- cbind(Loss5, midas_Loss5)
SSM_midas95 <- MCSprocedure(Loss = Losses95, alpha = 0.2, B = 5000, statistic = "Tmax")

Losses99 <- cbind(Loss, midas_Loss)
SSM_midas99 <- MCSprocedure(Loss = Losses99, alpha = 0.2, B = 5000, statistic = "Tmax")

