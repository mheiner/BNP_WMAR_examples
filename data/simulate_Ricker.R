rm(list=ls())

f = function(x, theta) x * exp( theta - x)
curve(f(x, 2.6), from=-0.2, to=7)

L = 10

TT = 10e3 
nburn = 100
ninit = 5
n = ninit + nburn + TT + L

sig = 0.09
theta = 2.6







### lognormal

set.seed(1)

y0 = numeric(n)
y0[1:ninit] = exp(rnorm(ninit, 0.0, sd=sig))

for (t in (ninit+1):n) {
	y0[t] = f(y0[t-2], theta) * exp(rnorm(1, 0.0, sd=sig))
}

y1 = y0[-c(1:(ninit+nburn))]

length(y1) == (TT + L)

yX = embed(y1, L+1)

y = yX[,1]
X = yX[,-1]

str(y)
str(X)

# nplot = 200
# plot.ts(y[1:nplot])

# par(mfrow=c(3,2))
# for (i in 1:5) plot(X[1:nplot,i], y[1:nplot], main=paste("Lag", i))

# library("rgl")
# plot3d(X[1:nplot,3], X[1:nplot,2], X[1:nplot,1], type='p', xlab="lag 2", ylab="lag 1", zlab="y")
# dev.off()


## for validation

nprime = 1000
tprime = sort(sample(1001:10e3, nprime, replace=FALSE))

y_valid = y[tprime]
X_valid = X[tprime,]

nsim_valid = 2000

Y_valid = matrix(NA, nrow=nsim_valid, ncol=nprime)
for (i in 1:nprime) {
	Y_valid[,i] = f(X_valid[i,2], theta) * exp(rnorm(nsim_valid, 0.0, sd=sig))
}

i = 2
hist(Y_valid[,i], breaks=30, main=paste("y[t-2] =", X_valid[i,2]))

save(file="Ricker_sigle_lognormal.rda", y, X, y_valid, X_valid, Y_valid, tprime)



### lognormal on two lags

set.seed(2)

y0 = numeric(n)
y0[1:ninit] = exp(rnorm(ninit, 0.0, sd=sig))

for (t in (ninit+1):n) {
	y0[t] = f(y0[t-2], theta) * exp(rnorm(1, 0.0, sd=sig*(y0[t-1])))
}


y1 = y0[-c(1:(ninit+nburn))]

length(y1) == (TT + L)

yX = embed(y1, L+1)

y = yX[,1]
X = yX[,-1]

str(y)
str(X)

# nplot = 70
# plot.ts(y[1:nplot])

# par(mfrow=c(3,2))
# for (i in 1:5) plot(X[1:nplot,i], y[1:nplot], main=paste("Lag", i))

# library("rgl")
# plot3d(X[1:nplot,3], X[1:nplot,2], X[1:nplot,1], type='p', xlab="lag 2", ylab="lag 1", zlab="y")
# dev.off()


# check empirical distributions
# x1val = 0.5
# x2val = 1.0
# indx = which( sqrt((X[,1] - 0.5)^2 + (X[,2] - 1.0)^2) < 0.5 )
# length(indx)
# hist(y[indx], freq=F)
# curve(dlnorm(x, log(x2val) + theta - x2val, sig*x1val), 0, 10, add=T, col="blue")
# dlnorm(0:10, log(x2val) + theta - x2val, sig*x1val)

## for validation

nprime = 1000
tprime = sort(sample(1001:10e3, nprime, replace=FALSE))

y_valid = y[tprime]
X_valid = X[tprime,]

nsim_valid = 2000

Y_valid = matrix(NA, nrow=nsim_valid, ncol=nprime)
for (i in 1:nprime) {
	Y_valid[,i] = f(X_valid[i,2], theta) * exp(rnorm(nsim_valid, 0.0, sd=sig*X_valid[i,1]))
}

# i = 50
# (i = which.max(y_valid))
# hist(Y_valid[,i], breaks=30, main=paste("y[t-1] =", X_valid[i,1], "y[t-2] =", X_valid[i,2]))

# library("rgl")
# plot3d(X_valid[,2], X_valid[,1], y_valid, type='p', xlab="lag 2", ylab="lag 1", zlab="y")
# dev.off()



save(file="Ricker_twolag_lognormal.rda", y, X, y_valid, X_valid, Y_valid, tprime)







### original

set.seed(3)

y0 = numeric(n)
y0[1:ninit] = exp(rnorm(ninit, 0.0, sd=sig))

for (t in (ninit+1):n) {
	y0[t] = f(y0[t-2], theta) + rnorm(1, 0.0, sd=sig)
}

y1 = y0[-c(1:(ninit+nburn))]

length(y1) == (TT + L)

yX = embed(y1, L+1)

y = yX[,1]
X = yX[,-1]

str(y)
str(X)

# nplot = 200
# plot.ts(y[1:nplot])

# par(mfrow=c(3,2))
# for (i in 1:5) plot(X[1:nplot,i], y[1:nplot], main=paste("Lag", i))

# library("rgl")
# plot3d(X[1:nplot,3], X[1:nplot,2], X[1:nplot,1], type='p', xlab="lag 2", ylab="lag 1", zlab="y")
# dev.off()


## for validation

nprime = 1000
tprime = sort(sample(1001:10e3, nprime, replace=FALSE))

y_valid = y[tprime]
X_valid = X[tprime,]

nsim_valid = 2000

Y_valid = matrix(NA, nrow=nsim_valid, ncol=nprime)
for (i in 1:nprime) {
	Y_valid[,i] = f(X_valid[i,2], theta) + rnorm(nsim_valid, 0.0, sd=sig)
}

# i = 2
# hist(Y_valid[,i], breaks=30, main=paste("y[t-2] =", X_valid[i,2]))

save(file="Ricker_single.rda", y, X, y_valid, X_valid, Y_valid, tprime)


 
rm(list=ls())
