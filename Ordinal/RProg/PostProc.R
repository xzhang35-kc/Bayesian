#---------------------------------------------------
# Post Process of Results
#---------------------------------------------------

#---------------------------------------------------
# load library
library(boa)

#---------------------------------------------------
# set up working directory
setwd("C:/CloudDrive/Dropbox/ZhangXiao/JASA/SubmittedCode/DataResults")

#---------------------------------------------------
# setup the parameters
Inum <- 10000 # number of iterations
N <- 50 # number of samples
K <- 2  # number of covariates
Rep <- 5 # dimension of outcome (number of repeated measures)
Ordcut <- 4# number of levels for ordinal variable
 
#---------------------------------------------------
# obtain results
Beta <- matrix(scan("CorrOrdBeta.dat"), ncol = K, byrow = T)
b.m <- apply(X = Beta, MARGIN = 2, FUN = mean)
b.m
b.sd <- apply(X = Beta, MARGIN = 2, FUN = sd)
b.sd

R <- matrix(scan("CorrOrdR.dat"), ncol = Rep * Rep, byrow = T)
R.m <- apply(X = R, MARGIN = 2, FUN = mean)
R.mean <-  matrix(R.m, Rep, Rep)
R.mean
R.sd <- apply(X = R, MARGIN = 2, FUN = sd)
R.sd <- matrix(R.sd, Rep, Rep)
R.sd

WR <- matrix(scan("CorrOrdSig.dat"), ncol = Rep * Rep, byrow = T)
WR.m <- apply(X = WR, MARGIN = 2, FUN = mean)
WR.mean <- matrix(WR.m, Rep, Rep)
WR.mean

WR.sd <- apply(X = WR, MARGIN = 2, FUN = sd)
WR.sd <- matrix(WR.sd, Rep, Rep)
WR.sd

Gama <- matrix(scan("CorrOrdGama.dat"), ncol = Rep * Ordcut, byrow = T)
Gama.m <- apply(X = Gama, MARGIN = 2, FUN = mean)
Gama.sd <- apply(X = Gama, MARGIN = 2, FUN = sd)
G1 <- round(Gama.m, 2)
G1
G2 <- round(Gama.sd, 2)
G2

#---------------------------------------------------
# Process results
# correlation
rS12 <- as.matrix(R[, 2])
rS13 <- as.matrix(R[, 3])
rS14 <- as.matrix(R[, 4])
rS15 <- as.matrix(R[, 5])
rS23 <- as.matrix(R[, 8])
rS24 <- as.matrix(R[, 9])
rS25 <- as.matrix(R[,10])
rS34 <- as.matrix(R[,14])
rS35 <- as.matrix(R[,15])
rS45 <- as.matrix(R[,20])

# conefficients
bS1 <- rep(0, Inum)
bS2 <- rep(0, Inum)
gS11 <- rep(0, Inum)
gS12 <- rep(0, Inum)
gS21 <- rep(0, Inum)
gS22 <- rep(0, Inum)
gS31 <- rep(0, Inum)
gS32 <- rep(0, Inum)
gS41 <- rep(0, Inum)
gS42 <- rep(0, Inum)
gS51 <- rep(0,Inum)
gS52 <- rep(0,Inum)

for(i in 1:Inum) 
{
	bS1[i]  <- Beta[i, 1] / sqrt(WR[i, 1])
	bS2[i]  <- Beta[i, 2] / sqrt(WR[i, 1])	
	gS11[i] <- Gama[i, 2] / sqrt(WR[i, 1])
	gS12[i] <- Gama[i, 3] / sqrt(WR[i, 1])
	gS21[i] <- Gama[i, 6] / sqrt(WR[i, 7])
	gS22[i] <- Gama[i, 7] / sqrt(WR[i, 7])
	gS31[i] <- Gama[i, 10] / sqrt(WR[i, 13])
	gS32[i] <- Gama[i, 11] / sqrt(WR[i, 13])
	gS41[i] <- Gama[i, 14] / sqrt(WR[i, 19])
	gS42[i] <- Gama[i, 15] / sqrt(WR[i, 19])
	gS51[i] <- Gama[i, 18] / sqrt(WR[i, 25])
	gS52[i] <- Gama[i, 19] / sqrt(WR[i, 25])
}

BetaM <- c(mean(bS1), mean(bS2))
BetaSD <- c(sd(bS1), sd(bS2))
BetaSOut <- round(data.frame(mean = BetaM, sd = BetaSD), 2)

GamaM <- c(mean(gS11), mean(gS12), mean(gS21), mean(gS22), mean(gS31), mean(gS32), 
mean(gS41), mean(gS42), mean(gS51), mean(gS52))
GamaSD <- c(sd(gS11), sd(gS12), sd(gS21), sd(gS22), sd(gS31), sd(gS32), 
sd(gS41), sd(gS42), sd(gS51), sd(gS52))
GamaSOut <- round(rbind(GamaM, GamaSD), 2)

RSMout <- round(R.mean, 2)
RSMout
RSSDout <- round(R.sd, 2)

WRSMout <- round(WR.mean, 2)
WRSMout
WRSSDout <- round(WR.sd, 2)

#---------------------------------------------------
# Output Table 1 for MH-ID with N-ID
BetaSOut

#---------------------------------------------------
# Output Figure 1 for MH-ID with N-ID
par(mfrow = c(3, 4))
# coefficent of intercept, beta_0
a2 <- acf(bS1, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(beta[0]), ylab = "ACF", lty = 1, lwd = 3)
# coefficent of intercept, beta_0
a2 <-  acf(bS2, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(beta[1]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r12
a2 <- acf(rS12, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[12]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r13
a2 <- acf(rS13, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[13]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r14
a2 <- acf(rS14, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[14]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r15
a2 <- acf(rS15, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[15]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r23
a2 <- acf(rS23, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[23]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r24
a2 <- acf(rS24, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[24]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r25
a2 <- acf(rS25, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[25]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r34
a2 <- acf(rS34, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[34]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r35
a2 <- acf(rS35, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[35]), ylab = "ACF", lty = 1, lwd = 3)
# correlation r45
a2 <- acf(rS45, plot = FALSE)
plot(a2$acf, type = "l", frame = FALSE, pch = 19, col = "black", xlab = bquote(r[45]), ylab = "ACF", lty = 1, lwd = 3)
