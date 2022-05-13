#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

setwd("C:/ZXiao/MVP-Ordinal/SimuStudy-X5")

source("simu-funs.R")

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#------------------Load files-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


setwd("C:/ZXiao/MVP-Ordinal/SimuStudy-X5/Simu-2000/Outputs/SavedResults-5000")

#------------------------------------------------------------------------------------
#------------------------PX-GS-------------------------------------------------------
#------------------------------------------------------------------------------------


load("beta.gs.RData")
load("gama.gs.RData")
load("r.gs.RData")
load("sigma.gs.RData")

load("beta.acfgs.RData")
load("gama.acfgs.RData")
load("r.acfgs.RData")
load("sigma.acfgs.RData")

load("beta.resgs.RData")
load("gama.resgs.RData")
load("r.resgs.RData")
load("sigma.resgs.RData")

load("Bgs.RData")
load("Ggs.RData")
load("Rgs.RData")
load("Sgs.RData")


#------------------------------------------------------------------------------------
#------------------------PX-GSM------------------------------------------------------
#------------------------------------------------------------------------------------


load("beta.gsm.RData")
load("gama.gsm.RData")
load("r.gsm.RData")
load("sigma.gsm.RData")

load("beta.acfgsm.RData")
load("gama.acfgsm.RData")
load("r.acfgsm.RData")
load("sigma.acfgsm.RData")

load("beta.resgsm.RData")
load("gama.resgsm.RData")
load("r.resgsm.RData")
load("sigma.resgsm.RData")

load("Bgsm.RData")
load("Ggsm.RData")
load("Rgsm.RData")
load("Sgsm.RData")


#------------------------------------------------------------------------------------
#------------------------PX-MH-------------------------------------------------------
#------------------------------------------------------------------------------------

load("beta.mh.RData")
load("gama.mh.RData")
load("r.mh.RData")
load("sigma.mh.RData")

load("beta.acfmh.RData")
load("gama.acfmh.RData")
load("r.acfmh.RData")
load("sigma.acfmh.RData")

load("beta.resmh.RData")
load("gama.resmh.RData")
load("r.resmh.RData")
load("sigma.resmh.RData")

load("Bmh.RData")
load("Gmh.RData")
load("Rmh.RData")
load("Smh.RData")


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
N <- 2000
burn.in = 5000
S <- 500 
K <- 5
P <- 5 
G <- 20000 
vt <- "1,1,1,1,1"
alpha = 0.05
NC = c(4, 4, 4, 4, 4) #cutpoints
cut = sum(NC) #cutpoint number

theta = c(0.5, 0.5, 1, 1, 1.5)
cuts <- list(c(-Inf, 0, 0.5, 1, Inf),
	c(-Inf, 0, 0.5, 1, Inf),
	c(-Inf, 0, 1, 2, Inf),
	c(-Inf, 0, 1, 2, Inf),
	c(-Inf, 0, 1.5, 2.5, Inf))

C <- rbind(c(1.0000, 0.500, 0.25, 0.125, 0.0625),
	c(0.5000, 1.000, 0.50, 0.250, 0.1250),
	c(0.2500, 0.500, 1.00, 0.500, 0.2500),
	c(0.1250, 0.250, 0.50, 1.000, 0.5000),
	c(0.0625, 0.125, 0.25, 0.500, 1.0000))

Tcuts <- c(0, 0.5, 1, 10000, 0, 0.5, 1, 10000, 0, 1, 2, 10000, 0, 1, 2, 10000, 0, 1.5, 2.5, 10000)

D = diag(c(1, 1, 1, 1, 1))

Sig = D%*%C%*%D

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
Nbeta = c("b1", "b2")
Ngama = c("g11", "g12", "g21", "g22", "g31", "g32", "g41", "g42", "g51", "g52")
Nr <- c("r12", "r13", "r14", "r15", 
"r23", "r24", "r25",
"r34", "r35", 
"r45")

Nw <- c("w11","w22", "w33", "w44", "w55")


Lgama = c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)
Lr = c(2:5, 8:10, 14, 15, 20)
Lw = diag(t(matrix(c(1:25), 5, 5)))
#----------------------------------------------------------------------
#----------------------------------------------------------------------
gg = length(Ngama)
rr = length(Nr)
ww = length(Nw)
indexG = c(11, 12, 21, 22, 31, 32, 41, 42, 51, 52)
indexR = c(12:15, 23:25, 34:35, 45)
indexW = c(11, 22, 33, 44, 55)

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------
#--------------------------------------PX-GS---------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

beta.biasgs = round(apply(X=beta.resgs$abs.bias, FUN=mean, MARGIN=2),3)
beta.biasgs
beta.rbiasgs = round(apply(X=beta.resgs$rel.bias, FUN=mean, MARGIN=2),3)
beta.rbiasgs
beta.rmsegs = round(apply(X=beta.resgs$rmse, FUN=mean, MARGIN=2),3)
beta.rmsegs

gama.bias = round(apply(X=gama.resgs$abs.bias, FUN=mean, MARGIN=2),3)
gama.bias = gama.bias[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.biasgs = matrix(gama.bias, 1, 10, byrow=T)
gama.biasgs
gama.rbias = round(apply(X=gama.resgs$rel.bias, FUN=mean, MARGIN=2),3)
gama.rbias = gama.rbias[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.rbiasgs = matrix(gama.rbias, 1, 10, byrow=T)
gama.rbiasgs
gama.rmse = round(apply(X=gama.resgs$rmse, FUN=mean, MARGIN=2),3)
gama.rmse = gama.rmse[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.rmsegs = matrix(gama.rmse, 1, 10, byrow=T)
gama.rmsegs

r.bias = round(apply(X=r.resgs$abs.bias, FUN=mean, MARGIN=2),3)
r.biasgs = matrix(r.bias, K, K, byrow=T)
r.biasgs
r.rbias = round(apply(X=r.resgs$rel.bias, FUN=mean, MARGIN=2),3)
r.rbiasgs = matrix(r.rbias, K, K, byrow=T)
r.rbiasgs
r.rmse = round(apply(X=r.resgs$rmse, FUN=mean, MARGIN=2),3)
r.rmsegs = matrix(r.rmse, K, K, byrow=T)
r.rmsegs

sigma.bias = round(apply(X=sigma.resgs$abs.bias, FUN=mean, MARGIN=2),3)
sigma.biasgs = matrix(sigma.bias, K, K, byrow=T)
sigma.biasgs
sigma.rbias = round(apply(X=sigma.resgs$rel.bias, FUN=mean, MARGIN=2),3)
sigma.rbiasgs = matrix(sigma.rbias, K, K, byrow=T)
sigma.rbiasgs
sigma.rmse = round(apply(X=sigma.resgs$rmse, FUN=mean, MARGIN=2),3)
sigma.rmsegs = matrix(sigma.rmse, K, K, byrow=T)
sigma.rmsegs

beta.cigs = apply(X=beta.resgs$ci.cov, FUN=mean, MARGIN=2)
beta.cigs
gama1.cigs = apply(X=gama.resgs$ci.cov, FUN=mean, MARGIN=2)
gama1.cigs = gama1.cigs[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.cigs = matrix(gama1.cigs, 1, 10, byrow=T)
gama.cigs
r1.cigs = apply(X=r.resgs$ci.cov, FUN=mean, MARGIN=2)
r.cigs  = matrix(r1.cigs, K, K, byrow=T)
r.cigs
sigma1.cigs = apply(X=sigma.resgs$ci.cov, FUN=mean, MARGIN=2)
sigma.cigs  = matrix(sigma1.cigs, K, K, byrow=T)
sigma.cigs


#----------------------------------------------------------------------------------------------------------------------------
#--------------------------------------PX-GSM--------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

beta.biasgsm = round(apply(X=beta.resgsm$abs.bias, FUN=mean, MARGIN=2),3)
beta.biasgsm
beta.rbiasgsm = round(apply(X=beta.resgsm$rel.bias, FUN=mean, MARGIN=2),3)
beta.rbiasgsm
beta.rmsegsm = round(apply(X=beta.resgsm$rmse, FUN=mean, MARGIN=2),3)
beta.rmsegsm

gama.bias = round(apply(X=gama.resgsm$abs.bias, FUN=mean, MARGIN=2),3)
gama.bias = gama.bias[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.biasgsm = matrix(gama.bias, 1, 10, byrow=T)
gama.biasgsm
gama.rbias = round(apply(X=gama.resgsm$rel.bias, FUN=mean, MARGIN=2),3)
gama.rbias = gama.rbias[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.rbiasgsm = matrix(gama.rbias, 1, 10, byrow=T)
gama.rbiasgsm
gama.rmse = round(apply(X=gama.resgsm$rmse, FUN=mean, MARGIN=2),3)
gama.rmse = gama.rmse[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.rmsegsm = matrix(gama.rmse, 1, 10, byrow=T)
gama.rmsegsm

r.bias = round(apply(X=r.resgsm$abs.bias, FUN=mean, MARGIN=2),3)
r.biasgsm = matrix(r.bias, K, K, byrow=T)
r.biasgsm
r.rbias = round(apply(X=r.resgsm$rel.bias, FUN=mean, MARGIN=2),3)
r.rbiasgsm = matrix(r.rbias, K, K, byrow=T)
r.rbiasgsm
r.rmse = round(apply(X=r.resgsm$rmse, FUN=mean, MARGIN=2),3)
r.rmsegsm = matrix(r.rmse, K, K, byrow=T)
r.rmsegsm

sigma.bias = round(apply(X=sigma.resgsm$abs.bias, FUN=mean, MARGIN=2),3)
sigma.biasgsm = matrix(sigma.bias, K, K, byrow=T)
sigma.biasgsm
sigma.rbias = round(apply(X=sigma.resgsm$rel.bias, FUN=mean, MARGIN=2),3)
sigma.rbiasgsm = matrix(sigma.rbias, K, K, byrow=T)
sigma.rbiasgsm
sigma.rmse = round(apply(X=sigma.resgsm$rmse, FUN=mean, MARGIN=2),3)
sigma.rmsegsm = matrix(sigma.rmse, K, K, byrow=T)
sigma.rmsegsm

beta.cigsm = apply(X=beta.resgsm$ci.cov, FUN=mean, MARGIN=2)
beta.cigsm
gama1.cigsm = apply(X=gama.resgsm$ci.cov, FUN=mean, MARGIN=2)
gama1.cigsm = gama1.cigsm[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.cigsm = matrix(gama1.cigsm, 1, 10, byrow=T)
gama.cigsm
r1.cigsm = apply(X=r.resgsm$ci.cov, FUN=mean, MARGIN=2)
r.cigsm  = matrix(r1.cigsm, K, K, byrow=T)
r.cigsm
sigma1.cigsm = apply(X=sigma.resgsm$ci.cov, FUN=mean, MARGIN=2)
sigma.cigsm  = matrix(sigma1.cigsm, K, K, byrow=T)
sigma.cigsm


#----------------------------------------------------------------------------------------------------------------------------
#--------------------------------------PX-MH---------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

beta.biasmh = round(apply(X=beta.resmh$abs.bias, FUN=mean, MARGIN=2),3)
beta.biasmh
beta.rbiasmh = round(apply(X=beta.resmh$rel.bias, FUN=mean, MARGIN=2),3)
beta.rbiasmh
beta.rmsemh = round(apply(X=beta.resmh$rmse, FUN=mean, MARGIN=2),3)
beta.rmsemh

gama.bias = round(apply(X=gama.resmh$abs.bias, FUN=mean, MARGIN=2),3)
gama.bias = gama.bias[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.biasmh = matrix(gama.bias, 1, 10, byrow=T)
gama.biasmh
gama.rbias = round(apply(X=gama.resmh$rel.bias, FUN=mean, MARGIN=2),3)
gama.rbias = gama.rbias[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.rbiasmh = matrix(gama.rbias, 1, 10, byrow=T)
gama.rbiasmh
gama.rmse = round(apply(X=gama.resmh$rmse, FUN=mean, MARGIN=2),3)
gama.rmse = gama.rmse[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.rmsemh = matrix(gama.rmse, 1, 10, byrow=T)
gama.rmsemh

r.bias = round(apply(X=r.resmh$abs.bias, FUN=mean, MARGIN=2),3)
r.biasmh = matrix(r.bias, K, K, byrow=T)
r.biasmh
r.rbias = round(apply(X=r.resmh$rel.bias, FUN=mean, MARGIN=2),3)
r.rbiasmh = matrix(r.rbias, K, K, byrow=T)
r.rbiasmh
r.rmse = round(apply(X=r.resmh$rmse, FUN=mean, MARGIN=2),3)
r.rmsemh = matrix(r.rmse, K, K, byrow=T)
r.rmsemh

sigma.bias = round(apply(X=sigma.resmh$abs.bias, FUN=mean, MARGIN=2),3)
sigma.biasmh = matrix(sigma.bias, K, K, byrow=T)
sigma.biasmh
sigma.rbias = round(apply(X=sigma.resmh$rel.bias, FUN=mean, MARGIN=2),3)
sigma.rbiasmh = matrix(sigma.rbias, K, K, byrow=T)
sigma.rbiasmh
sigma.rmse = round(apply(X=sigma.resmh$rmse, FUN=mean, MARGIN=2),3)
sigma.rmsemh = matrix(sigma.rmse, K, K, byrow=T)
sigma.rmsemh

beta.cimh = apply(X=beta.resmh$ci.cov, FUN=mean, MARGIN=2)
beta.cimh
gama1.cimh = apply(X=gama.resmh$ci.cov, FUN=mean, MARGIN=2)
gama1.cimh = gama1.cimh[c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
gama.cimh = matrix(gama1.cimh, 1, 10, byrow=T)
gama.cimh
r1.cimh = apply(X=r.resmh$ci.cov, FUN=mean, MARGIN=2)
r.cimh  = matrix(r1.cimh, K, K, byrow=T)
r.cimh
sigma1.cimh = apply(X=sigma.resmh$ci.cov, FUN=mean, MARGIN=2)
sigma.cimh  = matrix(sigma1.cimh, K, K, byrow=T)
sigma.cimh

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#----------------End of Analysis--------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#--------------------Start Output-------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------



setwd("C:/ZXiao/MVP-Ordinal/SimuStudy-X5/Simu-2000/Outputs/Burn-5000-1")

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------
sink("Table 2.txt")
#------------------------------------------------------------------------------------------------
x <- "Beta for PX-GS"
print(x)
beta.gs
x <- "Beta for PX-GSM"
print(x)
beta.gsm
x <- "Beta for PX-MH"
print(x)
beta.mh
x <- "Gamma for PX-GS"
print(x)
gama.gs
x <- "Gamma for PX-GSM"
print(x)
gama.gsm
x <- "Gamma for PX-MH"
print(x)
gama.mh
x <- "R for PX-GS"
print(x)
r.gs
x <- "R for PX-GSM"
print(x)
r.gsm
x <- "R for PX-MH"
print(x)
r.mh
x <- "Sigma for PX-GS"
print(x)
sigma.gs
x <- "Sigma for PX-GSM"
print(x)
sigma.gsm
x <- "Sigma for PX-MH"
print(x)
sigma.mh
sink()
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
sink("PID-CI.txt")
#------------------------------------------------------------------------------------------------
x <- "Beta for PX-GS"
print(x)
beta.cigs
x <- "Beta for PX-GSM"
print(x)
beta.cigsm
x <- "Beta for PX-MH"
print(x)
beta.cimh
x <- "Gamma for PX-GS"
print(x)
gama.cigs
x <- "Gamma for PX-GSM"
print(x)
gama.cigsm
x <- "Gamma for PX-MH"
print(x)
gama.cimh
x <- "R for PX-GS"
print(x)
r.cigs
x <- "R for PX-GSM"
print(x)
r.cigsm
x <- "R for PX-MH"
print(x)
r.cimh
x <- "Sigma for PX-GS"
print(x)
sigma.cigs
x <- "Sigma for PX-GSM"
print(x)
sigma.cigsm
x <- "Sigma for PX-MH"
print(x)
sigma.cimh
sink()
#---------------------------------------------------------------------------------------------------------------------------
