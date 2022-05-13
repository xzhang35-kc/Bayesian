
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#------------------Start the analysis--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

setwd("C:/ZXiao/MVP-Ordinal/Submission/JCGS/R-prog/RLMS")

wd <- "C:/ZXiao/MVP-Ordinal/Submission/JCGS/R-prog/RLMS/" 

N <- 2439

#----------------------------------------------------------------------
#----------------------------------------------------------------------


burn.in = 0
Output="PID-B0-"
Plot="PID-B0-"


#----------------------------------------------------------------------
#----------------------------------------------------------------------
K <- 7 # number of repeats
P <- 5 # number of covariates
G = 20000 # interation number
CP = 4 # number of cutpoints
S = 1
s = 1
alpha = 0.05 # quantiles
NC = c(4,4,4,4,4,4,4) #cutpoints
cut = sum(NC) #cutpoint number
#----------------------------------------------------------------------
#----------------------------------------------------------------------

Nbeta = c("b1", "b2", "b3", "b4", "b5")
Ngama = c("g11", "g12", "g21", "g22", "g31", "g32", "g41", "g42", "g51", "g52", "g61", "g62", "g71", "g72")
Nr <- c("r12", "r13", "r14", "r15", "r16", "r17", 
"r23", "r24", "r25", "r26", "r27", 
"r34", "r35", "r36", "r37", 
"r45", "r46", "r47", 
"r56", "r57", 
"r67")

Nw <- c("w11","w22", "w33", "w44", "w55", "w66", "w77")
Lgama = c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27)
Lr = c(2:7, 10:14, 18:21, 26:28, 34:35, 42)
Lw = diag(t(matrix(c(1:49), 7, 7)))
#----------------------------------------------------------------------
#----------------------------------------------------------------------

a = length(Ngama)
b = length(Nr)
c = length(Nw)
indexG = c(11, 12, 21, 22, 31, 32, 41, 42, 51, 52, 61, 62, 71, 72)
indexR = c(12:17, 23:27, 34:37, 45:47, 56:57, 67)
indexW = c(11, 22, 33, 44, 55, 66, 77)
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# load functions
source("simu-funs.R")
#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#---PX-GSM-------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

m.id <- "PX-GSM" 
file.pre <- paste(wd, "res-pid-px-gsm/", "res", sep = "")
#----------------------------------------------------------------------
# step 5: run R to obtain results
#----------------------------------------------------------------------
beta.list <- as.list(1:S)
gama.list <- as.list(1:S)
r.list <- as.list(1:S)
bbeta.list <- as.list(1:S)
ggama.list <- as.list(1:S)

sigma.list <- as.list(1:S)
for(s in 1:S) {
	res.file <- paste(file.pre, s, "-", "CorrOrdBeta.dat", sep = "")
	beta.list[[s]] <- scan(res.file)
	bbeta.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdGama.dat", sep = "")
	gama.list[[s]] <- scan(res.file)
	ggama.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdR.dat", sep = "")
	r.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdSig.dat", sep = "")
	sigma.list[[s]] <- scan(res.file)
}
#-----------------------------------------------------------------------
# step 6: adjust estimates
#-----------------------------------------------------------------------
if(m.id == "PX-GS" || m.id == "PX-GSM") {
	for(s in 1:S) {
		cat("results adjustment", s,"of", S, "...\n")
		var.vec <- find.var(sigma.list[[s]], G = G, K = K)
		sd.vec <- sqrt(var.vec)
		beta.list[[s]] <- beta.adj(beta.list[[s]], sd.vec, G = G, K = K)
		gama.list[[s]] <- gama.adj(gama.list[[s]], sd.vec, G = G, K = K, NC = NC)
	}
}
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# step 7: obtain summary statisic for parameters
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# summary statistic for beta	
beta.stat <- sapply(beta.list, FUN = para.stat, G = G, burn.in = burn.in)
beta.stat <- t(beta.stat)
beta.out1 <- round(apply(beta.stat, 2, mean), 4)
beta.out <- matrix(beta.out1, 2, P, byrow=TRUE)
beta.out

# summary statistic for gama
gama.stat <- sapply(gama.list, FUN = para.stat, G = G, burn.in = burn.in)
gama.stat <- t(gama.stat)
gama.out1 <- round(apply(gama.stat, 2, mean), 2)
gama.out <- matrix(gama.out1, 2, cut, byrow=TRUE)
gama.out

# summary statistic for correlation matrix
r.stat <- sapply(r.list, FUN = para.stat, G = G, burn.in = burn.in)
r.stat <- t(r.stat)
r.out1 <- round(apply(r.stat, 2, mean), 2)
r.out <- matrix(r.out1, 2*K, K, byrow=TRUE)
r.out

# summary statistic for covariance matrix
sigma.stat <- sapply(sigma.list, FUN = para.stat, G = G)
sigma.stat <- t(sigma.stat)
w.out1 <- round(apply(sigma.stat, 2, mean), 2)
w.out <- matrix(w.out1, 2*K, K, byrow=TRUE)
w.out

# acf
beta.acf <- para.list.acf(beta.list, G = G, burn.in = burn.in)
gama.acf <- para.list.acf(gama.list, G = G, burn.in = burn.in)
r.acf <- para.list.acf(r.list, G = G, burn.in = burn.in)
sigma.acf <- para.list.acf(sigma.list, G = G, burn.in = burn.in)

# 95% CI
beta.ci = para.CI(beta.list[[1]], G = G, burn.in = burn.in, alpha = alpha)
beta.ci = round(matrix(unlist(beta.ci), ncol = 2), 3)
beta.ci

#----------------------------------------------------------------------
# step 8: Output and prepare acf plots for PX-GSM
#----------------------------------------------------------------------
beta.gsm <- beta.out
gama.gsm <- gama.out
r.gsm <- r.out
w.gsm <- w.out
beta.acfgsm <- beta.acf 
gama.acfgsm <- gama.acf
r.acfgsm <- r.acf
sigma.acfgsm <- sigma.acf
beta.cigsm <- beta.ci
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#----------Begin Diagnostics PX-GSM------------------------------------
#----------------------------------------------------------------------

beta.m <- matrix(beta.list[[1]], ncol=P, byrow=T)
head(beta.m)
apply(beta.m, MARGIN=2, FUN=mean)

beta.m = beta.m[(burn.in+1):G,]
apply(beta.m, MARGIN=2, FUN=mean)
Dbeta.gsm = beta.m

#----------------------------------------------------------------------
#----------Gama PX-GSM-------------------------------------------------
#----------------------------------------------------------------------
gama.m <- matrix(gama.list[[1]], ncol=K*CP, byrow=T)
head(gama.m)
apply(gama.m, MARGIN=2, FUN=mean)

gama.m = as.matrix(gama.m[(burn.in+1):G,Lgama])
apply(gama.m, MARGIN=2, FUN=mean)
Dgama.gsm = gama.m

#----------------------------------------------------------------------
#----------R PX-GSM----------------------------------------------------
#----------------------------------------------------------------------
r.m <- matrix(r.list[[1]], ncol=K*K, byrow=T)
head(r.m)
apply(r.m, MARGIN=2, FUN=mean)

r.m = r.m[(burn.in+1):G,Lr]
apply(r.m, MARGIN=2, FUN=mean)
Dr.gsm = r.m
#----------------------------------------------------------------------
#----------W PX-GSM-----------------------------------------------------
#----------------------------------------------------------------------
w.m <- matrix(sigma.list[[1]], ncol=K*K, byrow=T)
head(w.m)
apply(w.m, MARGIN=2, FUN=mean)

w.m = as.matrix(w.m[(burn.in+1):G, Lw])
apply(w.m, MARGIN=2, FUN=mean)
Dw.gsm = w.m
#----------------------------------------------------------------------
#----------End Diagnostics PX-GSM--------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#---PX-GS--------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

m.id <- "PX-GS" 
file.pre <- paste(wd, "res-pid-px-gs/", "res", sep = "")

beta.list <- as.list(1:S)
gama.list <- as.list(1:S)
r.list <- as.list(1:S)
bbeta.list <- as.list(1:S)
ggama.list <- as.list(1:S)

sigma.list <- as.list(1:S)
for(s in 1:S) {
	res.file <- paste(file.pre, s, "-", "CorrOrdBeta.dat", sep = "")
	beta.list[[s]] <- scan(res.file)
	bbeta.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdGama.dat", sep = "")
	gama.list[[s]] <- scan(res.file)
	ggama.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdR.dat", sep = "")
	r.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdSig.dat", sep = "")
	sigma.list[[s]] <- scan(res.file)
}
#-----------------------------------------------------------------------
# step 6: adjust estimates
#-----------------------------------------------------------------------
if(m.id == "PX-GS" || m.id == "PX-GSM") {
	for(s in 1:S) {
		cat("results adjustment", s,"of", S, "...\n")
		var.vec <- find.var(sigma.list[[s]], G = G, K = K)
		sd.vec <- sqrt(var.vec)
		beta.list[[s]] <- beta.adj(beta.list[[s]], sd.vec, G = G, K = K)
		gama.list[[s]] <- gama.adj(gama.list[[s]], sd.vec, G = G, K = K, NC = NC)
	}
}
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# step 7: obtain summary statisic for parameters
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# summary statistic for beta	
beta.stat <- sapply(beta.list, FUN = para.stat, G = G, burn.in = burn.in)
beta.stat <- t(beta.stat)
beta.out1 <- round(apply(beta.stat, 2, mean), 4)
beta.out <- matrix(beta.out1, 2, P*K, byrow=TRUE)
beta.out

# summary statistic for gama
gama.stat <- sapply(gama.list, FUN = para.stat, G = G, burn.in = burn.in)
gama.stat <- t(gama.stat)
gama.out1 <- round(apply(gama.stat, 2, mean), 2)
gama.out <- matrix(gama.out1, 2, cut, byrow=TRUE)
gama.out

# summary statistic for correlation matrix
r.stat <- sapply(r.list, FUN = para.stat, G = G, burn.in = burn.in)
r.stat <- t(r.stat)
r.out1 <- round(apply(r.stat, 2, mean), 2)
r.out <- matrix(r.out1, 2*K, K, byrow=TRUE)
r.out

# summary statistic for covariance matrix
sigma.stat <- sapply(sigma.list, FUN = para.stat, G = G)
sigma.stat <- t(sigma.stat)
w.out1 <- round(apply(sigma.stat, 2, mean), 2)
w.out <- matrix(w.out1, 2*K, K, byrow=TRUE)
w.out

# acf
beta.acf <- para.list.acf(beta.list, G = G, burn.in = burn.in)
gama.acf <- para.list.acf(gama.list, G = G, burn.in = burn.in)
r.acf <- para.list.acf(r.list, G = G, burn.in = burn.in)
sigma.acf <- para.list.acf(sigma.list, G = G, burn.in = burn.in)

# 95% CI
beta.ci = para.CI(beta.list[[1]], G = G, burn.in = burn.in, alpha = alpha)
beta.ci = round(matrix(unlist(beta.ci), ncol = 2), 3)
beta.ci

#----------------------------------------------------------------------
# step 8: Output and prepare acf plots for PX-GS
#----------------------------------------------------------------------
beta.gs <- beta.out
gama.gs <- gama.out
r.gs <- r.out
w.gs <- w.out
beta.acfgs <- beta.acf 
gama.acfgs <- gama.acf
r.acfgs <- r.acf
sigma.acfgs <- sigma.acf
beta.cigs <- beta.ci

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------



w.gsm <- matrix(sigma.list[[1]], ncol=K*K, byrow=T)

beta.gsm <- matrix(beta.list[[1]], ncol=P, byrow=T)
bbeta.gsm <- matrix(bbeta.list[[1]], ncol=P, byrow=T)

gama.gsm <- matrix(gama.list[[1]], ncol=K*CP, byrow=T)
ggama.gsm <- matrix(ggama.list[[1]], ncol=K*CP, byrow=T)

r.gsm <- matrix(r.list[[1]], ncol=K*K, byrow=T)



#----------------------------------------------------------------------
#----------Begin Diagnostics PX-GS-------------------------------------
#----------------------------------------------------------------------

beta.m <- matrix(beta.list[[1]], ncol=P, byrow=T)
head(beta.m)
apply(beta.m, MARGIN=2, FUN=mean)

beta.m = beta.m[(burn.in+1):G,]
apply(beta.m, MARGIN=2, FUN=mean)
Dbeta.gs = beta.m

#----------------------------------------------------------------------
#----------Gama PX-GS--------------------------------------------------
#----------------------------------------------------------------------
gama.m <- matrix(gama.list[[1]], ncol=CP*K, byrow=T)
head(gama.m)
apply(gama.m, MARGIN=2, FUN=mean)

gama.m = as.matrix(gama.m[(burn.in+1):G,Lgama])
apply(gama.m, MARGIN=2, FUN=mean)
Dgama.gs = gama.m
#----------------------------------------------------------------------
#----------R PX-GS-----------------------------------------------------
#----------------------------------------------------------------------
r.m <- matrix(r.list[[1]], ncol=K*K, byrow=T)
head(r.m)
apply(r.m, MARGIN=2, FUN=mean)

r.m = r.m[(burn.in+1):G,Lr]
apply(r.m, MARGIN=2, FUN=mean)
Dr.gs = r.m

#----------------------------------------------------------------------
#----------W PX-GS-----------------------------------------------------
#----------------------------------------------------------------------
w.m <- matrix(sigma.list[[1]], ncol=K*K, byrow=T)
head(w.m)
apply(w.m, MARGIN=2, FUN=mean)

w.m = as.matrix(w.m[(burn.in+1):G, Lw])
apply(w.m, MARGIN=2, FUN=mean)
Dw.gs = w.m
#----------------------------------------------------------------------
#----------End Diagnostics PX-GS---------------------------------------
#----------------------------------------------------------------------


w.gs <- matrix(sigma.list[[1]], ncol=K*K, byrow=T)

beta.gs <- matrix(beta.list[[1]], ncol=P, byrow=T)
bbeta.gs <- matrix(bbeta.list[[1]], ncol=P, byrow=T)

gama.gs <- matrix(gama.list[[1]], ncol=K*CP, byrow=T)
ggama.gs <- matrix(ggama.list[[1]], ncol=K*CP, byrow=T)

r.gs <- matrix(r.list[[1]], ncol=K*K, byrow=T)

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#--------------------Figure 5-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
tiff("Figure5.tif", width = 4800, height = 4500,res = 800, compression = "lzw")

par(mfrow=c(2,4), mar=c(4,3,1,1))

plot(Dw.gsm[,1], type="l", ylim=c(0,4),  col = "black",  xlab = bquote(sigma[~.(11)]), ylab=" ")
lines(Dw.gs[,1], type="l", col="grey50")

plot(Dw.gsm[,2], type="l", ylim=c(0,4),  col = "black",  xlab = bquote(sigma[~.(22)]), ylab=" ")
lines(Dw.gs[,2], type="l", col="grey50")

plot(Dw.gsm[,3], type="l", ylim=c(0,4),  col = "black",  xlab = bquote(sigma[~.(33)]), ylab=" ")
lines(Dw.gs[,3], type="l", col="grey50")

plot(Dw.gsm[,4], type="l", ylim=c(0,4),  col = "black",  xlab = bquote(sigma[~.(44)]), ylab=" ")
lines(Dw.gs[,4], type="l", col="grey50")

plot(Dw.gsm[,5], type="l", ylim=c(0,4),  col = "black",  xlab = bquote(sigma[~.(55)]), ylab=" ")
lines(Dw.gs[,5], type="l", col="grey50")

plot(Dw.gsm[,6], type="l", ylim=c(0,4),  col = "black",  xlab = bquote(sigma[~.(66)]), ylab=" ")
lines(Dw.gs[,6], type="l", col="grey50")

plot(Dw.gsm[,7], type="l", ylim=c(0,4),  col = "black",  xlab = bquote(sigma[~.(77)]), ylab=" ")
lines(Dw.gs[,7], type="l", col="grey50")


dev.off()
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------



