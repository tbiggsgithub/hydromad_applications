#  Read in PEQ zoo file, try hydromad for Guejito
library(hydromad)

outdir = "G:/mydocuments/SDSU_classes/2017F/780/w12_IHACRES/Rfiles/"
fname = "Guejito.PEQout.csv"

xin = read.zoo(paste0(outdir,fname))
xyplot(xin,with.P=TRUE)

# Choose one 10-year period to calibrate the model.  which year range should you use?

x = window(xin,start=as.Date(c("1972-10-01"),end=as.Date("1981-09-30")))
# Check to make sure data subsetting worked
xyplot(x,with.P=TRUE)

modx <- hydromad(x, sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))
fitx <- fitByOptim(modx)
fitx
summary(fitx)
xyplot(fitx, with.P = TRUE)
qqmath(fitx,type = c("l", "g"),scales=list(y=list(log=TRUE)))

# Extract fitted values and plot in separate time series
fitx.values = fitx$fitted.values
fitx.values[fitx.values==0]=0.0001
plot(fitx.values,log="y")
lines(x$Q,col="blue")

#  Run calibrated model on validation period
xval = window(xin,start=as.Date(c("1955-10-01"),end=as.Date("1965-09-30")))
fitx.val = update(fitx,newdata=xval)
fitx.val
summary(fitx.val)
xyplot(fitx.val,log="y")
xval$Q[xval$Q==0]=0.0001
plot(xval$Q,log="y")
xval.fitted = fitx.val$fitted.values
xval.fitted[fitx.values==0]=0.0001
lines(xval.fitted,col="blue")
qqmath(fitx.val,type = c("l", "g"),scales=list(y=list(log=TRUE)))
#  Estimate delay
xd=x
names(xd) = c("U","Q","E")
estimateDelay(DATA=x)

#  Add "l" and "p" for semi-arid watersheds
modx.w.l <- hydromad(x, sma = "cwi", l=c(0.1,99),routing = "expuh",
                     tau_s = c(2,100))
modx.w.l
fitx.w.l <- fitByOptim(modx.w.l)  # ONly works if take out v_s range.
fitx.w.l
summary(fitx.w.l)
xyplot(fitx.w.l, with.P = TRUE, type = c("l", "g"))
qqmath(fitx.w.l,type = c("l", "g"),scales=list(y=list(log=TRUE)))


# Fit by Dream
library(dream)
fitx.dream <- fitByDream(modx, control = list(ndraw = 500))
summary(fitx.dream)
xyplot(fitx.dream,with.P=TRUE)
qqmath(fitx.dream,type = c("l", "g"),scales=list(y=list(log=TRUE)))

#  Determine how parameters are trading off
cov2cor(vcov(fitx.dream))
#  Which parameters are highly correlated?

symnum(cov2cor(vcov(fitx.dream)))

str(fitx.dream$fit.result)

## plot log-likelihood value convergence over time
xyplot(window(optimtrace(fitx.dream, raw = TRUE), start = 50),
       superpose = TRUE, auto.key = FALSE,
       xlab = "function evaluations", ylab = "neg. log likelihood")

# Find "behavioral" models.  Doesn't work
fitglu <- defineFeasibleSet(fitx.dream, threshold = 0.05,
                            glue.quantiles = c(0.05, 0.95))

# Try gridded sampling
fitx.sim = simulate(modx,nsim=144,sampletype = "latin")
coef(fitx.sim)
summary(fitx.sim)

fitx.sim1 = simulate(modx,nsim=144,sampletype = "all")
levelplot(result ~ tw + f, fitx.sim1, cex = 2,
          panel = panel.levelplot.points,
          main = "R Squared (of sq.rt. data) over parameter space") +
  layer_(panel.2dsmoother(...))


mod1 <- update(modx, routing = "armax", rfit = list("ls", order = c(2,1)))
               
sim1 <- simulate(mod1, 144, sampletype = "all", FUN = objFunVal,
                                objective = ~ nseStat(Q, X, trans = sqrt))               
levelplot(result ~ tw + f, sim1, cex = 2,
          panel = panel.levelplot.points,
          main = "R Squared (of sq.rt. data) over parameter space") +
  layer_(panel.2dsmoother(...))



