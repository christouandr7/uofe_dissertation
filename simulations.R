library(fastICA)
library(ProDenICA)

#####################FIRST STEPS - 2 VARIABLES UNIFORM DISTRIBUTION################
#2 variables following uniform distribution with 200 samples
x <- runif(200) - 1/2
y <- runif(200) - 1/2
xy <- cbind(x,y)
par(mfrow = c(1, 1))
plot(xy, xlab = "signal1", ylab = "signal2")

A <- matrix(c(runif(4)), ncol = 2, nrow = 2)	#random mixing matrix
mixture <- xy %*% A	#mixing procedure (A*s)
plot(mixture, xlab = "mixture1", ylab = "mixture2")

#use fastICA to unmix the mixture signals
unmix <- fastICA(mixture, 2, alg.typ = "parallel", fun = "logcosh", 
								 alpha = 1, method = "C", row.norm = FALSE, maxit = 200, 
								 tol = 0.0001, verbose = TRUE)

#plot the ics
par(mfrow = c(1,1))
plot(unmix$S[,1], unmix$S[,2], xlab = "", main = "Independent components after ICA")


######################TIME SERIES SIGNALS##################
#two time series signals randomly generated
par(mfrow = c(2,1))
r1 <- runif(100, -100, 100)
r1 <- c(r1, runif(100, -10, 10))
plot(r1, type = "l", main = "Signal1", xlab = "Time", ylab = "")

r2 <- runif(50, -10, 10)
r2 <- c(r2, runif(100, -100, 100))
r2 <- c(r2, runif(50, -10, 10))
plot(r2, type = "l", main = "Signal2", xlab = "Time", ylab = "")

A <- matrix(c(runif(4)), ncol = 2, nrow = 2) #random mixing matrix
s <- cbind(r1, r2)
mixx <- t(A %*% t(s))	#mixing procedure - multiply A with source signals to get mixtures

#plots
par(mfrow = c(1, 2),
		oma = c(2, 2, 1, 1), # two rows of text at the outer left and bottom margin
		mar = c(1, 1, 1.5, 1))
plot(mixx[,1], xlab = "Time",type = "l")
title("Mixture1", line = 0.5)
plot(mixx[,2], xlab = "Time",type = "l")
title("Mixture2", line = 0.5)

#use fastICA to unmix the mixture signals and get the independent components
unmix <- fastICA(mixx, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1, 
								 method = "C", row.norm = FALSE, maxit = 200, 
								 tol = 0.000001, verbose = TRUE)
#plot the independent components (time series signals)
par(mfrow = c(1, 2),
		oma = c(2, 2, 1, 1), # two rows of text at the outer left and bottom margin
		mar = c(1, 1, 1.5, 1))
plot(unmix$S[,1], xlab = "Time",type = "l")
title("IC1", line = 0.5)
plot(unmix$S[,2], xlab = "Time",type = "l")
title("IC2", line = 0.5)

#use the correlation matrix
abs(cor(s, unmix$S))


######################MORE TIME SERIES SIGNALS - 4 SIGNALS##################
par(mfrow = c(2,2))
r1 <- runif(100, -100, 100)
r1 <- c(r1, runif(100, -10, 10))
plot(r1, type = "l", main = "Signal 1", xlab = "Time", ylab = "")

r2 <- runif(25, -10, 10)
r2 <- c(r2, runif(100, -100, 100))
r2 <- c(r2, runif(75, -10, 10))
plot(r2, type = "l", main = "Signal 2", xlab = "Time", ylab = "")

r3 <- runif(50, -10, 10)
r3 <- c(r3, runif(100, -100, 100))
r3 <- c(r3, runif(50, -10, 10))
plot(r3, type = "l", main = "Signal 3", xlab = "Time", ylab = "")

r4 <- runif(75, -10, 10)
r4 <- c(r4, runif(100, -100, 100))
r4 <- c(r4, runif(25, -10, 10))
plot(r4, type = "l", main = "Signal 4", xlab = "Time", ylab = "")


A <- matrix(c(runif(16)), ncol = 4, nrow = 4) #random mixing matrix

s <- cbind(r1, r2)
s <- cbind(s, r3)
s <- cbind(s, r4)
mix4 <- t(A %*% t(s))	#multiply source signals with mixing matrix (mixign procedure)

#plot the mixtures
par(mfrow = c(2, 2),
		oma = c(2, 2, 1, 1), # two rows of text at the outer left and bottom margin
		mar = c(1, 1, 2.5, 1))
for (i in c(1:4)){
	plot(mix4[,i], xlab = "Time", ylab = "", type = "l")
	titleee <- paste("Mixture ",i, sep = "")
	title(titleee, line = 0.5)
}

#run fastICA to get the independent components - unmix the singals
unmix <- fastICA(mix4, 4, alg.typ = "parallel", fun = "logcosh", alpha = 1, 
								 method = "C", row.norm = FALSE, maxit = 200, 
								 tol = 0.000001, verbose = TRUE)

#plot the independent components
for (i in c(1:4)){
	plot(unmix$S[,i], xlab = "Time", ylab = "", type = "l")
	titleee <- paste("IC ",i, sep = "")
	title(titleee, line = 0.5)
}

#correlation matrix
round(abs(cor(s, unmix$S)),3)

#use ProDenICA to unmix the mixtures
unmixProDenICA <- ProDenICA(mix4, k = 4, whiten = TRUE, maxit = 200)

for (i in c(1:4)){
	plot(unmixProDenICA$s[,i], xlab = "Time", ylab = "", type = "l")
	titleee <- paste("IC ",i, sep = "")
	title(titleee, line = 0.5)
}

#get the correlation matrix of source signals with the independent components
#estimated by the ProDenIDA algorithm
abs(cor(s, unmixProDenICA$s))



################################Count time passed for the two algorithms################
#compare the computational times for fastICA and ProDenICA
time_start_fastICA <- Sys.time()
unmix <- fastICA(mix4, 4, alg.typ = "parallel", fun = "logcosh", alpha = 1, 
								 method = "C", row.norm = FALSE, maxit = 200, 
								 tol = 0.000001, verbose = TRUE)
time_end_fastICA <- Sys.time()


time_start_ProDenICA <- Sys.time()
unmixProDenICA <- ProDenICA(mix4, k = 4, whiten = TRUE, maxit = 200)
time_end_ProDenICA <- Sys.time()


print(time_end_fastICA - time_start_fastICA)
print(time_end_ProDenICA - time_start_ProDenICA)

