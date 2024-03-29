


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


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


burn.in = 5000
Output="Table3"


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
sigma.list <- as.list(1:S)
for(s in 1:S) {
	res.file <- paste(file.pre, s, "-", "CorrOrdBeta.dat", sep = "")
	beta.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdGama.dat", sep = "")
	gama.list[[s]] <- scan(res.file)
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
# beta.out and beat.ci
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# summary statistic for beta	
beta.stat <- sapply(beta.list, FUN = para.stat, G = G, burn.in = burn.in)
beta.stat <- t(beta.stat)
beta.out1 <- round(apply(beta.stat, 2, mean), 4)
beta.out <- matrix(beta.out1, 2, P, byrow=TRUE)
beta.out

# 95% CI
beta.ci = para.CI(beta.list[[1]], G = G, burn.in = burn.in, alpha = alpha)
beta.ci = round(matrix(unlist(beta.ci), ncol = 2), 3)
beta.ci


beta.gsm <- beta.out
beta.cigsm <- beta.ci

#----------------------------------------------------------------------
#----------------------------------------------------------------------
#---PX-GS--------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

m.id <- "PX-GS" 
file.pre <- paste(wd, "res-pid-px-gs/", "res", sep = "")

#----------------------------------------------------------------------
# step 5: run R to obtain results
#----------------------------------------------------------------------
beta.list <- as.list(1:S)
gama.list <- as.list(1:S)
r.list <- as.list(1:S)
sigma.list <- as.list(1:S)
for(s in 1:S) {
	res.file <- paste(file.pre, s, "-", "CorrOrdBeta.dat", sep = "")
	beta.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdGama.dat", sep = "")
	gama.list[[s]] <- scan(res.file)
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
# beta.out and beat.ci
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# summary statistic for beta	
beta.stat <- sapply(beta.list, FUN = para.stat, G = G, burn.in = burn.in)
beta.stat <- t(beta.stat)
beta.out1 <- round(apply(beta.stat, 2, mean), 4)
beta.out <- matrix(beta.out1, 2, P, byrow=TRUE)
beta.out

# 95% CI
beta.ci = para.CI(beta.list[[1]], G = G, burn.in = burn.in, alpha = alpha)
beta.ci = round(matrix(unlist(beta.ci), ncol = 2), 3)
beta.ci

beta.gs <- beta.out
beta.cigs <- beta.ci


#----------------------------------------------------------------------
#----------------------------------------------------------------------
#---PX-MH--------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

m.id <- "PX-MH" 
file.pre <- paste(wd, "res-pid-px-mh/", "res", sep = "")
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------
beta.list <- as.list(1:S)
gama.list <- as.list(1:S)
r.list <- as.list(1:S)
sigma.list <- as.list(1:S)
for(s in 1:S) {
	res.file <- paste(file.pre, s, "-", "CorrOrdBeta.dat", sep = "")
	beta.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdGama.dat", sep = "")
	gama.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdR.dat", sep = "")
	r.list[[s]] <- scan(res.file)
	res.file <- paste(file.pre, s, "-", "CorrOrdSig.dat", sep = "")
	sigma.list[[s]] <- scan(res.file)
}

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# beta.out and beat.ci
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# summary statistic for beta	
beta.stat <- sapply(beta.list, FUN = para.stat, G = G, burn.in = burn.in)
beta.stat <- t(beta.stat)
beta.out1 <- round(apply(beta.stat, 2, mean), 4)
beta.out <- matrix(beta.out1, 2, P, byrow=TRUE)
beta.out
# 95% CI
beta.ci = para.CI(beta.list[[1]], G = G, burn.in = burn.in, alpha = alpha)
beta.ci = round(matrix(unlist(beta.ci), ncol = 2), 3)
beta.ci


beta.mh <- beta.out
beta.cimh <- beta.ci





#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#---------------Table 3----------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------

setwd(""C:/ZXiao/MVP-Ordinal/Submission/JCGS/R-prog/RLMS")

texname<-paste0(Output, "Output", ".txt")
sink(file=texname)
x <- "Beta for PX-GS"
print(x)
beta.gs
x <- "Beta for PX-GSM"
print(x)
beta.gsm
x <- "Beta for PX-MH"
print(x)
beta.mh
x <- "Beta CI for PX-GS"
print(x)
beta.cigs
x <- "Beta CI for PX-GSM"
print(x)
beta.cigsm
x <- "Beta CI for PX-MH"
print(x)
beta.cimh
sink()






















