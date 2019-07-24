## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)


## ------------------------------------------------------------------------
library(forecast)
library(MTS)
setwd("C:/Graduate_School_at_Cornell/Logistics/Crses/555/555project")
load ("Anomaly.RData")


## ----warning=FALSE, results='hide', fig1,fig.width=5.5,fig.height=3.5,fig.cap="\\label{fig:fig1}ARIMA (left) and VARMA (right) model fitting to the Meterological station data points"----
anom <- Anomaly[,1]
d_anom <- diff(anom)
dd_anom <- diff(d_anom)
plot(dd_anom, col = 'red', lty = 3, ylab = "Meterological station")
lines(d_anom, col = 'blue', lty = 3)
lines(anom)
legend("bottomright", col = c("black", "blue", "red"),
legend = c("Meterological stations measurement",
           "First differenced meterological stations measurements",
           "Second differenced meterological stations measurements"),
lty = c(1,3,3), cex = 0.53)


## ----warning=FALSE, results='hide'---------------------------------------
ndiffs(anom, alpha = 0.05, test = "adf", type = "trend")
ndiffs(anom, alpha = 0.05, test = "pp", type = "trend")


## ----warning=FALSE, results='hide'---------------------------------------
ps <- 0:11
qs <- 0:1
aic <- matrix(nrow = length(ps), ncol = length(qs))
n <- length(anom)

for (i in 1:length(ps)) {
  # cat("i=", i, "\n")
  for (j in 1:length(qs)) {
    p <- ps[i]
    q <- qs[j]
    uml <- try(arima(anom, order = c(p, 1, q), method = "ML", xreg = time(anom)))
    
    if (inherits(uml, "try-error") | !uml$code == 0) {
      aic[i, j] <- NA
    } else {
    sig.sq <- uml$sigma2
    aic[i, j] <- log(sig.sq) + (n + 2*(p + q + 1))/n
    }
  }
}

aic.pq <- which(aic == min(aic), arr.ind = TRUE)
bestp <- ps[aic.pq[1]]
bestq <- qs[aic.pq[2]]


## ----warning=FALSE, results='hide'---------------------------------------
cutoff <- 1500
to.use <- 1:cutoff
not.to.use <- (cutoff + 1):length(anom)
anom_trim <- window(anom, end=time(anom)[cutoff]) 
n.ahead <- length(not.to.use)

y <- anom
y[not.to.use] <- NA
t <- time(anom)


## ----warning=FALSE,results='hide'----------------------------------------
# ARIMA(bestp,1,0)
fit_arima_pd0 <- arima(anom_trim, order = c(bestp, 1, 0), method = "ML", xreg = (time(anom_trim)))
pred_pd0 <- predict(fit_arima_pd0, n.ahead = n.ahead, se = TRUE,
                newxreg = (max(time(anom_trim)) + 1:n.ahead/12))

# ARIMA(0,1,bestq)
fit_arima_0dq <- arima(anom_trim, order = c(0, 1, bestq), method = "ML", xreg = (time(anom_trim)))
pred_0dq <- predict(fit_arima_0dq, n.ahead = n.ahead, se = TRUE,
                newxreg = (max(time(anom_trim)) + 1:n.ahead/12))

# ARIMA(bestp,1,bestq)
fit_arima_pdq <- arima(anom_trim, order = c(bestp, 1, bestq), method = "ML", xreg = (time(anom_trim)))
pred_pdq <- predict(fit_arima_pdq, n.ahead = n.ahead, se = TRUE,
                newxreg = (max(time(anom_trim)) + 1:n.ahead/12))

AIC(fit_arima_pd0)
AIC(fit_arima_0dq)
AIC(fit_arima_pdq)


## ----eval = FALSE, echo = FALSE------------------------------------------
## ocean = Anomaly[,2]
## d_ocean = diff(ocean)
## dd_ocean = diff(d_ocean)
## par(mfrow = c(1,3))
## plot(ocean, main = "Anomaly", ylab = "Ocean")
## plot(d_ocean, main = "d_Anomaly", ylab = "Ocean")
## plot(dd_ocean, main = "Anomaly", ylab = "Ocean")
## 
## par(mfrow = c(1,2))
## #par(mar = rep(2, 4))
## acf_docean <- acf(d_ocean, lag.max = 15, plot = FALSE)
## pacf_docean <- acf(d_ocean, lag.max = 15, type = "partial", plot = FALSE)
## plot(acf_docean, main = "ACF", ylab = "ACF")
## plot(pacf_docean, main = "PACF", ylab = "PACF")


## ----warning=FALSE, results='hide'---------------------------------------
ocean <- Anomaly[,2]
ocean_trim <- window(ocean, end=time(ocean)[cutoff]) 
n.ahead <- length(not.to.use)

anom.demean <- anom - mean(anom)
ocean.demean <- ocean - mean(ocean)
#X.pro <- cbind(anom_trim, ocean_trim)
X.pro <- cbind(anom.demean, ocean.demean)
X.pro <- X.pro[1:cutoff, ]
Z.pro <- 1:nrow(X.pro)

linfit <- lm(X.pro~Z.pro)
ps <- 0:11
qs <- 0:1
#aic <- matrix(nrow = length(ps), ncol = length(qs))
#bic <- matrix(nrow = length(ps), ncol = length(qs))
aic.pro <- matrix(c(-8.795,-8.795,-8.895,-8.912,-8.930,-8.930,-8.930,-8.929,-8.948,-8.205,-8.933,-8.952,-8.960,-8.959,-8.958,-8.955,-8.952,-8.962),nrow=9,ncol=2,byrow = FALSE)
bic.pro <- matrix(c(-8.780,-8.780,-8.867,-8.870,-8.874,-8.860,-8.845,-8.830,-8.835,-8.191,-8.904,-8.909,-8.903,-8.888,-8.873,-8.856,-8.839,-8.834), nrow=9,ncol=2,byrow = FALSE)

n <- length(anom)


## ----eval = FALSE, echo = FALSE------------------------------------------
## # Find the AIC/BIC minimizing VARMA(p,q) model
## for (i in 1:length(ps)) {
##   cat("i=", i, "\n")
##   for (j in 1:length(qs)) {
##     pv <- ps[i]
##     qv <- qs[j]
##     armaMTS.pro <- VARMA(linfit$residuals, p = pv, q = qv, include.mean = FALSE)
##     aic[i, j] <- armaMTS.pro$aic
##     bic[i, j] <- armaMTS.pro$bic
##   }
## }


## ----warning=FALSE, results='hide'---------------------------------------
aic.noNA <- aic.pro[1:9,]
bic.noNA <- bic.pro[1:9,]
aic.pq <- which(aic.noNA == min(aic.noNA), arr.ind = TRUE)
bic.pq <- which(bic.noNA == min(bic.noNA), arr.ind = TRUE)
bestp <- ps[bic.pq[1]]
bestq <- qs[bic.pq[2]]

armaMTS.best <- VARMA(linfit$residuals, p = bestp, q = bestq, include.mean = FALSE)
varma.pred <- VARMApred(armaMTS.best, h = length(ocean) - nrow(X.pro))

# plot the forecasts with VARMA(2,1) model
# par(mfrow = c(2, 1))


## ----warning=FALSE, results='hide', fig2,fig.width=7,fig.height=3,fig.cap="\\label{fig:fig2}ARIMA (left) and VARMA (right) model fitting to the Meterological station data points"----
fig2start <- 1400
t_anom <- (max(time(anom_trim)) + 1:n.ahead/12)

par(mfrow = c(1, 2))
par(mar = c(3,4,1,1))
plot(time(anom)[fig2start:length(anom)],anom[fig2start:length(anom)], col = "gray", type = "l",
     xlab = "Time", ylab = "Meterological Station",
     ylim = range(-2,3))
lines(time(anom_trim)[fig2start:length(anom_trim)], anom_trim[fig2start:length(anom_trim)])
lines(t_anom, pred_pd0$pred, col = "blue")
lines(t_anom, pred_pd0$pred + qnorm(0.975)*pred_pd0$se, col = "blue", lty = 2)
lines(t_anom, pred_pd0$pred - qnorm(0.975)*pred_pd0$se, col = "blue", lty = 2)

lines(t_anom, pred_0dq$pred, col = "red")
lines(t_anom, pred_0dq$pred + qnorm(0.975)*pred_0dq$se, col = "red", lty = 2)
lines(t_anom, pred_0dq$pred - qnorm(0.975)*pred_0dq$se, col = "red", lty = 2)

lines(t_anom, pred_pdq$pred, col = "purple")
lines(t_anom, pred_pdq$pred + qnorm(0.975)*pred_pdq$se, col = "purple", lty = 2)
lines(t_anom, pred_pdq$pred - qnorm(0.975)*pred_pdq$se, col = "purple", lty = 2)

legend("topleft", col = c("blue", "blue"),
legend = c("ARIMA(8,1,0) Model predictions", 
           "95% CIs of ARIMA(8,1,0) predictions"),lty = c(1,2,1,2,1,2), cex = 0.5)
legend("bottomleft", col = c("red","red","purple","purple"),
legend = c("ARIMA(0,1,1) Model predictions", 
           "95% CIs of ARIMA(0,1,1) predictions", 
           "ARIMA(8,1,1) Model predictions", 
           "95% CIs of ARIMA(8,1,1) predictions"),lty = c(1,2,1,2,1,2), cex = 0.5)

par(mar = c(3,4,1,1))
lineartrend <- linfit$coefficients[, 1][1] + linfit$coefficients[, 1][2]*(1:length(anom))
plot(time(anom)[fig2start:length(anom)],lineartrend[fig2start:length(anom)], col = "red", type = "l", xlab = "Time", ylab = "Meterological station", ylim = range(-2,3))
lines(time(anom)[fig2start:length(anom)],anom[fig2start:length(anom)], col = "gray", type = "l")
lines(anom.demean[fig2start:length(anom)], col = "gray", type = "l")
lines(time(anom)[fig2start:nrow(X.pro)], X.pro[, 1][fig2start:nrow(X.pro)])
lines(time(anom)[(nrow(X.pro) + 1):length(anom)], varma.pred$pred[, 1] + linfit$coefficients[, 1][1] + linfit$coefficients[, 1][2]*time(anom)[nrow(X.pro) + 1:nrow(varma.pred$pred)], col = "blue")
lines(time(anom)[(nrow(X.pro) + 1):length(anom)], varma.pred$pred[, 1] + linfit$coefficients[, 1][1] + linfit$coefficients[, 1][2]*time(anom)[nrow(X.pro) + 1:nrow(varma.pred$pred)] + qnorm(0.975)*varma.pred$se.err[, 1], col = "blue",
      lty = 2)
lines(time(anom)[(nrow(X.pro) + 1):length(anom)], varma.pred$pred[, 1] + linfit$coefficients[, 1][1] + linfit$coefficients[, 1][2]*time(anom)[nrow(X.pro) + 1:nrow(varma.pred$pred)] - qnorm(0.975)*varma.pred$se.err[, 1], col = "blue",
      lty = 2)
legend("topleft", col = c("blue", "blue", "red"),
legend = c("VARMA(2,1) Model predictions", 
           "95% CIs of VARMA(2,1) predictions", 
           "Linear fit over 90% of the time series data"), 
           lty = c(1,2,1), cex = 0.6)


## ----eval=FALSE, echo=FALSE----------------------------------------------
## # plot Ocean and forecasts
## par(mar = c(3,4,1,1))
## plot(time(ocean),linfit$coefficients[, 2][1] + linfit$coefficients[, 2][2]*(1:length(ocean)), col = "red", type = "l", xlab = "Time", ylab = "Ocean", ylim = range(ocean.demean, linfit$coefficients[, 2][1] + linfit$coefficients[, 2][2]*time(ocean)))
## lines(ocean.demean, col = "gray", type = "l")
## lines(time(ocean)[1:nrow(X.pro)], X.pro[, 2])
## lines(time(ocean)[(nrow(X.pro) + 1):length(ocean)], varma.pred$pred[, 2] + linfit$coefficients[, 2][1] + linfit$coefficients[, 2][2]*time(ocean)[nrow(X.pro) + 1:nrow(varma.pred$pred)], col = "blue")
## lines(time(ocean)[(nrow(X.pro) + 1):length(ocean)], varma.pred$pred[, 2] + linfit$coefficients[, 2][1] + linfit$coefficients[, 2][2]*time(ocean)[nrow(X.pro) + 1:nrow(varma.pred$pred)] + qnorm(0.975)*varma.pred$se.err[, 2], col = "blue",
##       lty = 2)
## lines(time(ocean)[(nrow(X.pro) + 1):length(ocean)], varma.pred$pred[, 2] + linfit$coefficients[, 2][1] + linfit$coefficients[, 2][2]*time(ocean)[nrow(X.pro) + 1:nrow(varma.pred$pred)] - qnorm(0.975)*varma.pred$se.err[, 2], col = "blue",
##       lty = 2)
## legend("topleft", col = c("blue", "blue", "red"),
## legend = c("VARMA(2,1) Model predictions",
##            "95% CIs of VARMA(2,1) predictions",
##            "Linear fit over 90% of the time series data"),
##            lty = c(1,2,1), cex = 0.6)
## #abline(a = linfit$coefficients[, 2][1], b = linfit$coefficients[, 2][2])


## ----warning=FALSE,results='hide'----------------------------------------
library(MARSS)
state.anom <- anom
t.state.anom <- time(state.anom)
n <- length(anom)
state.anom <- state.anom - mean(state.anom)
state.anom[(cutoff+1):n] <- NA
anom.demean <- anom - mean(anom)

month <- factor(round((round(time(state.anom), 3) - floor(time(state.anom)))*12))
covariates <- t(model.matrix(~month-1)[, -12])

# Basic model
model.basic <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                    Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                    x0=matrix("mu"), tinitx=0)
fit.basic <- MARSS(c(state.anom), model=model.basic, method = "kem")
fit.basic <- MARSS(c(state.anom), model=model.basic, method = "BFGS",
                   inits = fit.basic)

# put the covariates in the observation part of the equation
model.covy <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   D="unconstrained", d=covariates,
                   x0=matrix("mu"), tinitx=0 )
fit.covy <- MARSS(c(state.anom), model=model.covy, method = "kem")
fit.covy <- MARSS(c(state.anom), model=model.covy, method = "BFGS",
                  inits = fit.covy)

# put the covariates in the states part of the equation
model.covx <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   C="unconstrained", c=covariates,
                   x0=matrix("mu"), tinitx=0 )
fit.covx <- MARSS(c(state.anom), model=model.covx, method = "kem")
fit.covx <- MARSS(c(state.anom), model=model.covx,  method = "BFGS",
                  inits = fit.covx)

# put the covariates in the observation and states parts of the equation
model.covyx <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   D="unconstrained", d=covariates,
                   C="unconstrained", c=covariates,
                   x0=matrix("mu"), tinitx=0 )
fit.covyx <- MARSS(c(state.anom), model=model.covyx, method = "kem")
fit.covyx <- MARSS(c(state.anom), model=model.covyx,  method = "BFGS",
                  inits = fit.covyx)


## ----warning=FALSE, results='hide'---------------------------------------
AIC(fit.basic)
AIC(fit.covy)
AIC(fit.covx)
AIC(fit.covyx)


## ----warning=FALSE,results='hide'----------------------------------------

season <- factor(round((round(time(state.anom), 3) - floor(time(state.anom)))*3))
covariates.ss <- t(model.matrix(~season-1)[, -3])

# put the covariates in the observation part of the equation
model.covy.ss <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   D="unconstrained", d=covariates.ss,
                   x0=matrix("mu"), tinitx=0 )
fit.covy.ss <- MARSS(c(state.anom), model=model.covy.ss, method = "kem")
fit.covy.ss <- MARSS(c(state.anom), model=model.covy.ss, method = "BFGS",
                  inits = fit.covy.ss)

# put the covariates in the states part of the equation
model.covx.ss <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   C="unconstrained", c=covariates.ss,
                   x0=matrix("mu"), tinitx=0 )
fit.covx.ss <- MARSS(c(state.anom), model=model.covx.ss, method = "kem")
fit.covx.ss <- MARSS(c(state.anom), model=model.covx.ss, method = "BFGS",
                  inits = fit.covx.ss)

# put the covariates in the observation and states parts of the equation
model.covyx.ss <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   D="unconstrained", d=covariates.ss,
                   C="unconstrained", c=covariates.ss,
                   x0=matrix("mu"), tinitx=0 )
fit.covyx.ss <- MARSS(c(state.anom), model=model.covyx.ss, method = "kem")
fit.covyx.ss <- MARSS(c(state.anom), model=model.covyx.ss,  method = "BFGS",
                  inits = fit.covyx.ss)


## ----warning=FALSE, results='hide'---------------------------------------
AIC(fit.covy.ss)
AIC(fit.covx.ss)
AIC(fit.covyx.ss)


## ----warning=FALSE, fig3,fig.width=8,fig.height=4.5,fig.cap="\\label{fig:fig3}Univariate state space model fitting to the Meterological station data points"----
par(mfrow = c(1, 2))
par(mar = c(3,4,1,1))
plot(t.state.anom[(fig2start+1):n], anom.demean[(fig2start+1):n], col = "gray", type = "l", 
     xlab = "Time", ylab = "Meterological Station",
     main = "Yearly Seasonality",
     ylim = range(c(-1.8,1.5), na.rm = TRUE))
lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.basic$ytT)[is.na(state.anom)], col = "red")
lines(t[is.na(state.anom)], c(fit.basic$ytT + qnorm(0.975)*fit.basic$ytT.se)[is.na(state.anom)],
      col = "red", lty = 2)
lines(t[is.na(state.anom)], c(fit.basic$ytT - qnorm(0.975)*fit.basic$ytT.se)[is.na(state.anom)],
      col = "red", lty = 2)

lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.covy$ytT)[is.na(state.anom)], col = "blue")
lines(t[is.na(state.anom)], c(fit.covy$ytT + qnorm(0.975)*fit.covy$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)
lines(t[is.na(state.anom)], c(fit.covy$ytT - qnorm(0.975)*fit.covy$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)

lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.covx$ytT)[is.na(state.anom)], col = "green")
lines(t[is.na(state.anom)], c(fit.covx$ytT + qnorm(0.975)*fit.covx$ytT.se)[is.na(state.anom)],
      col = "green", lty = 2)
lines(t[is.na(state.anom)], c(fit.covx$ytT - qnorm(0.975)*fit.covx$ytT.se)[is.na(state.anom)],
      col = "green", lty = 2)

lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.covyx$ytT)[is.na(state.anom)], col = "purple")
lines(t[is.na(state.anom)], c(fit.covyx$ytT + qnorm(0.975)*fit.covyx$ytT.se)[is.na(state.anom)],
      col = "purple", lty = 2)
lines(t[is.na(state.anom)], c(fit.covyx$ytT - qnorm(0.975)*fit.covyx$ytT.se)[is.na(state.anom)],
      col = "purple", lty = 2)

legend("bottomleft", col = c("red","red", "blue","blue", "green","green", "purple","purple"),
legend = c("Predictions-Basic Model", 
           "95% CIs of predictions", 
           "Predictions-Monthly seasonality in observations", 
           "95% CIs of predictions",
           "Predictions-Monthly seasonality in states", 
           "95% CIs of predictions",
           "Predictions-Monthly seasonality in both", 
           "95% CIs of predictions"),lty = c(1,2,1,2,1,2,1,2), cex = 0.63)

par(mar = c(3,4,1,1))
plot(t.state.anom[(fig2start+1):n], anom.demean[(fig2start+1):n], col = "gray", type = "l", 
     xlab = "Time", ylab = "Meterological Station",
     main = "Quarterly Seasonality",
     ylim = range(c(-1.8,1.5), na.rm = TRUE))
lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.basic$ytT)[is.na(state.anom)], col = "red")
lines(t[is.na(state.anom)], c(fit.basic$ytT + qnorm(0.975)*fit.basic$ytT.se)[is.na(state.anom)],
      col = "red", lty = 2)
lines(t[is.na(state.anom)], c(fit.basic$ytT - qnorm(0.975)*fit.basic$ytT.se)[is.na(state.anom)],
      col = "red", lty = 2)

lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.covy.ss$ytT)[is.na(state.anom)], col = "blue")
lines(t[is.na(state.anom)], c(fit.covy.ss$ytT + qnorm(0.975)*fit.covy.ss$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)
lines(t[is.na(state.anom)], c(fit.covy.ss$ytT - qnorm(0.975)*fit.covy.ss$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)

lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.covx.ss$ytT)[is.na(state.anom)], col = "green")
lines(t[is.na(state.anom)], c(fit.covx.ss$ytT + qnorm(0.975)*fit.covx.ss$ytT.se)[is.na(state.anom)],
      col = "green", lty = 2)
lines(t[is.na(state.anom)], c(fit.covx.ss$ytT - qnorm(0.975)*fit.covx.ss$ytT.se)[is.na(state.anom)],
      col = "green", lty = 2)

lines(t.state.anom[!is.na(state.anom)], state.anom[!is.na(state.anom)])
lines(t.state.anom[is.na(state.anom)], c(fit.covyx.ss$ytT)[is.na(state.anom)], col = "purple")
lines(t[is.na(state.anom)], c(fit.covyx.ss$ytT + qnorm(0.975)*fit.covyx.ss$ytT.se)[is.na(state.anom)],
      col = "purple", lty = 2)
lines(t[is.na(state.anom)], c(fit.covyx.ss$ytT - qnorm(0.975)*fit.covyx.ss$ytT.se)[is.na(state.anom)],
      col = "purple", lty = 2)

legend("bottomleft", col = c("red","red", "blue","blue", "green","green", "purple","purple"),
legend = c("Predictions-Basic Model", 
           "95% CIs of predictions", 
           "Predictions-Quarterly seasonality in observations", 
           "95% CIs of predictions",
           "Predictions-Quarterly seasonality in states", 
           "95% CIs of predictions",
           "Predictions-Quarterly seasonality in both", 
           "95% CIs of predictions"),lty = c(1,2,1,2,1,2,1,2), cex = 0.63)


## ----warning=FALSE,results='hide'----------------------------------------
state.anom.resi <- anom
state.anom.resi <- state.anom.resi - mean(state.anom.resi)
t.state.anom <- time(state.anom)
n <- length(anom)
state.anom.resi[(cutoff+1):n] <- NA
X.resi <- state.anom.resi[1:cutoff]
Z.resi <- 1:cutoff
linfit.state <- lm(X.resi~Z.resi)
linfit.state.residuals <- state.anom.resi
linfit.state.residuals[1:cutoff] <- linfit.state$residuals

state.ocean.resi <- ocean
state.ocean.resi <- state.ocean.resi - mean(state.ocean.resi)
state.ocean.resi[(cutoff+1):n] <- NA
X.ocean.resi <- state.ocean.resi[1:cutoff]
linfit.state.ocean <- lm(X.ocean.resi~Z.resi)
linfit.state.ocean.residuals <- state.ocean.resi
linfit.state.ocean.residuals[1:cutoff] <- linfit.state.ocean$residuals

#month <- factor(round((round(time(state.anom), 3) - floor(time(state.anom)))*12))
#covariates.resi <- matrix(t(model.matrix(~month-1)[, -12]))

# Model the residuals
model.basic.resi <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                    Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                    x0=matrix("mu"), tinitx=0)
fit.basic.resi <- MARSS(c(linfit.state.residuals), model=model.basic.resi, method = "kem")
fit.basic.resi <- MARSS(c(linfit.state.residuals), model=model.basic.resi, method = "BFGS",
                   inits = fit.basic.resi)


## ----eval=FALSE, echo=FALSE----------------------------------------------
## # Model the residuals, put the covariates in the observation part of the equation
## # However, this model takes too long to fit
## model.covy.resi <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"),
##                    Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
##                    D="unconstrained", d=covariates.resi,
##                    x0=matrix("mu"), tinitx=0 )
## fit.covy.resi <- MARSS(c(linfit.state.residuals), model=model.covy.resi, method = "kem")
## fit.covy.resi <- MARSS(c(linfit.state.residuals), model=model.covy.resi, method = "BFGS",
##                   inits = fit.covy.resi)
## 


## ----warning=FALSE,results='hide'----------------------------------------
# covariates be the entire time vector
state.anom.cov <- anom
state.anom.cov <- state.anom.cov - mean(state.anom.cov)
n <- length(anom)
state.anom.cov[(cutoff+1):n] <- NA

covariates.cov <- matrix(1:length(anom), nrow = 1, ncol = length(anom))
# Had tried: covariates.cov <- matrix(1:length(anom), nrow = length(anom), ncol = 1), did not work

model.basic.cov <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                    Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                    D="unconstrained", d=covariates.cov,
                    x0=matrix("mu"), tinitx=0)
fit.basic.cov <- MARSS(c(state.anom.cov), model=model.basic.cov, method = "kem")
fit.basic.cov <- MARSS(c(state.anom.cov), model=model.basic.cov, method = "BFGS",
                   inits = fit.basic.cov)


## ----warning=FALSE, fig4,fig.width=8,fig.height=3,fig.cap="\\label{fig:fig4}Univariate state space model fitting to the Meterological station data points"----
par(mfrow = c(1, 2))
par(mar = c(3,4,1,1))
plot(time(anom)[(fig2start+1):n],(linfit.state$coefficients[1] + linfit$coefficients[2]*(1:length(anom)))[(fig2start+1):n], col = "red", type = "l", xlab = "Time", ylab = "Ocean", main = "Model the residuals of linear fit", ylim = range(-0.5,1.6), cex.main=0.8)
lines(time(anom)[(fig2start+1):n],anom.demean[(fig2start+1):n], col = "gray", type = "l")
lines(time(anom)[(fig2start+1):cutoff],anom.demean[(fig2start+1):cutoff], type = "l")

lines(time(anom)[is.na(state.anom)], c(fit.basic.resi$ytT)[is.na(state.anom)] + linfit.state$coefficients[1] + linfit.state$coefficients[2]*(1:length(anom))[is.na(state.anom)], col = "blue")
lines(time(anom)[is.na(state.anom)], linfit.state$coefficients[1] + linfit.state$coefficients[2]*(1:length(anom))[is.na(state.anom)] + c(fit.basic.resi$ytT)[is.na(state.anom)] + qnorm(0.975)*(fit.basic.resi$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)
lines(time(anom)[is.na(state.anom)], linfit.state$coefficients[1] + linfit.state$coefficients[2]*(1:length(anom))[is.na(state.anom)] + c(fit.basic.resi$ytT)[is.na(state.anom)] - qnorm(0.975)*(fit.basic.resi$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)

legend("topleft", col = c("blue","blue","red"),
legend = c("Univariate State Space Model predictions", 
           "95% CIs of predictions", 
           "Linear fit over 90% of the time series data"),lty = c(1,2,1), cex = 0.60)

par(mar = c(3,4,1,1))
plot(time(anom)[(fig2start+1):n], anom.demean[(fig2start+1):n], col = "gray", type = "l", 
     xlab = "Time", ylab = "Meterological Station",
     main = "Model with the covariates being entire time vector",
     ylim = range(c(-0.5,1.6), na.rm = TRUE),
     cex.main=0.8)
lines(time(anom)[(fig2start+1):cutoff], anom.demean[(fig2start+1):cutoff], type = "l")
lines(time(anom)[is.na(state.anom)], c(fit.basic.cov$ytT)[is.na(state.anom)]+ fit.basic.cov$coef["D.D"]*time(anom)[is.na(state.anom)], col = "blue")

lines(time(anom)[is.na(state.anom)], c(fit.basic.cov$ytT)[is.na(state.anom)]+ fit.basic.cov$coef["D.D"]*time(anom)[is.na(state.anom)] + qnorm(0.975)*c(fit.basic.cov$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)

lines(time(anom)[is.na(state.anom)], c(fit.basic.cov$ytT)[is.na(state.anom)]+ fit.basic.cov$coef["D.D"]*time(anom)[is.na(state.anom)] - qnorm(0.975)*c(fit.basic.cov$ytT.se)[is.na(state.anom)],
      col = "blue", lty = 2)

#lines(t[is.na(state.anom)], c(fit.basic.cov$ytT - qnorm(0.975)*fit.basic.cov$ytT.se)[is.na(state.anom)],
#      col = "blue", lty = 2)

legend("topleft", col = c("blue","blue"),
legend = c("Univariate State Space Model predictions", 
           "95% CIs of predictions"),lty = c(1,2), cex = 0.57)


## ----warning=FALSE, results='hide'---------------------------------------
# Multivariate state space model with both Meterological station and Ocean data
anom.demean <- anom - mean(anom)
ocean.demean <- ocean - mean(ocean)

Y.pro <- cbind(c(anom.demean), c(ocean.demean))

Y.pro[(cutoff+1):nrow(Y.pro), ] <- NA
covariates.multi <- matrix(1:length(anom), nrow = 1, ncol = length(anom))

model.multi <- list(
  ### Outcome Equation
  Z=matrix(c("a11", "a21", "a12", "a22"),2,2),
  A=matrix(0, 2, 1), 
  R="unconstrained", 
  ### State Equation
  B=matrix(c("phi11", "phi21", "phi12", "phi22"),2,2), 
  U=matrix(0, 2, 1),
  Q="unconstrained",
  D="unconstrained", d = covariates.multi,
  ### Initial value(s)
  x0=matrix("mu", 2, 1), 
  tinitx=1)

fit.multi <- MARSS(t(Y.pro), model=model.multi, method = "kem")
fit.multi <- MARSS(t(Y.pro), model=model.multi, method = "BFGS", inits = fit.multi)

AIC(fit.multi)


## ----warning=FALSE, results='hide'---------------------------------------
# Multivariate state space model with seasonality on the residuals of the linear fit
# Model the residuals
Y.mul <- cbind(c(linfit.state.residuals), c(linfit.state.ocean.residuals))
Y.mul[(cutoff+1):nrow(Y.mul), ] <- NA

Z.mul <- model.matrix(~factor(round(time(anom), 3) - floor(time(anom))))
Z.mul <- Z.mul[, -ncol(Z.mul)]


## ----eval=FALSE, echo=FALSE----------------------------------------------
## # Observation multivariate
## model.mul1 <-  list(
##   ### Outcome Equation
##   Z=matrix(c("a11", "a12"),2,1),
##   A=matrix(0, 2, 1),
##   R="unconstrained",
##   ### State Equation
##   B=matrix("phi"),
##   U=matrix(0, 1, 1),
##   Q=matrix("q"),
##   ### Initial value(s)
##   x0=matrix("mu"),
##   tinitx=1)
## fit.mul1 <- MARSS(t(Y.mul), model=model.mul1, method = "kem")
## fit.mul1 <- MARSS(t(Y.mul), model=model.mul1, method = "BFGS", inits = fit.mul1)
## 
## 
## # Both multivariate
## model.mul2 <-  list(
##   ### Outcome Equation
##   Z=matrix(c("a11", "a21", "a12", "a22"),2,2),
##   A=matrix(0, 2, 1),
##   R="unconstrained",
##   ### State Equation
##   B=matrix(c("phi11", "phi21", "phi12", "phi22"),2,2),
##   U=matrix(0, 2, 1),
##   Q="unconstrained",
##   ### Initial value(s)
##   x0=matrix("mu", 2, 1),
##   tinitx=1)
## fit.mul2 <- MARSS(t(Y.mul), model=model.mul2, method = "kem")
## fit.mul2 <- MARSS(t(Y.mul), model=model.mul2, method = "BFGS", inits = fit.mul2)
## 
## #three states
## model.mul3 <-  list(
##   ### Outcome Equation
##   Z=matrix(c("a11", "a12", "a13", "a21", "a22", "a23"),2,3, byrow = TRUE),
##   A=matrix(0, 2, 1),
##   R="unconstrained",
##   ### State Equation
##   B=matrix(c("phi11", "phi12", "phi13", "phi21", "phi22", "phi23", "phi31", "phi32", "phi33"),3,3, byrow=TRUE),
##   U=matrix(0, 3, 1),
##   Q="unconstrained",
##   ### Initial value(s)
##   x0=matrix("mu", 3, 1),
##   tinitx=1)
## fit.mul3 <- MARSS(t(Y.mul), model=model.mul3, method = "kem")
## fit.mul3 <- MARSS(t(Y.mul), model=model.mul3, method = "BFGS", inits = fit.mul3)
## 
## AIC(fit.mul1)
## AIC(fit.mul2)
## AIC(fit.mul3)


## ----warning=FALSE, results='hide'---------------------------------------
# for knitting
model.mul3 <-  list(
  ### Outcome Equation
  Z=matrix(c("a11", "a12", "a13", "a21", "a22", "a23"),2,3, byrow = TRUE),
  A=matrix(0, 2, 1), 
  R="unconstrained", 
  ### State Equation
  B=matrix(c("phi11", "phi12", "phi13", "phi21", "phi22", "phi23", "phi31", "phi32", "phi33"),3,3, byrow=TRUE),
  U=matrix(0, 3, 1), 
  Q="unconstrained", 
  ### Initial value(s)
  x0=matrix("mu", 3, 1), 
  tinitx=1)
fit.mul3 <- MARSS(t(Y.mul), model=model.mul3, method = "kem")
fit.mul3 <- MARSS(t(Y.mul), model=model.mul3, method = "BFGS", inits = fit.mul3)


## ----fig5,fig.width=8,fig.height=3,fig.cap="\\label{fig:fig5}Multivariate state space model fitting to the Meterological station data points"----
par(mfrow = c(1, 2))
par(mar = c(3,4,1,1))
plot(time(anom)[(fig2start+1):n],(linfit.state$coefficients[1] + linfit$coefficients[2]*(1:length(anom)))[(fig2start+1):n], col = "red", type = "l", xlab = "Time", ylab = "Meterological station", main = "Model the residuals of linear fit", ylim = range(-0.5,1.6), cex.main=0.8)

lines(time(anom)[(fig2start+1):n],anom.demean[(fig2start+1):n], col = "gray", type = "l")
lines(time(anom)[(fig2start+1):cutoff],anom.demean[(fig2start+1):cutoff], type = "l")

linear.mul <- linfit.state$coefficients[1] + linfit.state$coefficients[2]*(1:length(anom))[is.na(state.anom)]

lines(time(anom)[is.na(state.anom)], linear.mul + fit.mul3$ytT[1, is.na(Y.mul[, 1])], col = "blue")
lines(time(anom)[is.na(state.anom)], linear.mul + fit.mul3$ytT[1, is.na(Y.mul[, 1])] + qnorm(0.975)*fit.mul3$ytT.se[1, is.na(Y.mul[, 1])], col = "blue", lty = 2)
lines(time(anom)[is.na(state.anom)], linear.mul + fit.mul3$ytT[1, is.na(Y.mul[, 1])] - qnorm(0.975)*fit.mul3$ytT.se[1, is.na(Y.mul[, 1])], col = "blue", lty = 2)
legend("topleft", col = c("blue","blue","red"),
legend = c("Multivariate State Space Model predictions", 
           "95% CIs of predictions", 
           "Linear fit over 90% of the time series data"),lty = c(1,2,1), cex = 0.58)

par(mar = c(3,4,1,1))
plot(time(anom)[(fig2start+1):n], anom.demean[(fig2start+1):n], col = "gray", type = "l", 
     xlab = "Time", ylab = "Meterological Station",
     main = "Model with the covariates being entire time vector",
     ylim = range(c(-0.5,1.6), na.rm = TRUE),
     cex.main=0.8)
lines(time(anom)[(fig2start+1):cutoff],anom.demean[(fig2start+1):cutoff], type = "l")

lines(time(anom)[is.na(state.anom)], fit.multi$coef["D.Y1"]*time(anom)[is.na(state.anom)] + fit.multi$ytT[1, is.na(Y.mul[, 1])], col = "blue")
lines(time(anom)[is.na(state.anom)], fit.multi$coef["D.Y1"]*time(anom)[is.na(state.anom)] + fit.multi$ytT[1, is.na(Y.mul[, 1])] + qnorm(0.975)*fit.multi$ytT.se[1, is.na(Y.mul[, 1])], col = "blue", lty = 2)
lines(time(anom)[is.na(state.anom)], fit.multi$coef["D.Y1"]*time(anom)[is.na(state.anom)] + fit.multi$ytT[1, is.na(Y.mul[, 1])] - qnorm(0.975)*fit.multi$ytT.se[1, is.na(Y.mul[, 1])], col = "blue", lty = 2)

legend("topleft", col = c("blue","blue"),
legend = c("Mutltivariate State Space Model predictions", 
           "95% CIs of predictions"),lty = c(1,2), cex = 0.57)


## ----echo=FALSE, eval=FALSE----------------------------------------------
## 
## fig4.left <- c(fit.basic.resi$ytT)[is.na(state.anom)] + linfit.state$coefficients[1] + linfit.state$coefficients[2]*(1:length(anom))[is.na(state.anom)]
## 
## fig4.right <- c(fit.basic.cov$ytT)[is.na(state.anom)]+ fit.basic.cov$coef["D.D"]*time(anom)[is.na(state.anom)]
## 
## fig5.left <- linear.mul + fit.mul3$ytT[1, is.na(Y.mul[, 1])]
## fig5.right <- fit.multi$coef["D.Y1"]*time(anom)[is.na(state.anom)] + fit.multi$ytT[1, is.na(Y.mul[, 1])]
## 
## mean((fig4.left- anom.demean[(cutoff+1):n])^2)
## mean((fig4.right- anom.demean[(cutoff+1):n])^2)
## mean((fig5.left- anom.demean[(cutoff+1):n])^2)
## mean((fig5.right- anom.demean[(cutoff+1):n])^2)
## 

