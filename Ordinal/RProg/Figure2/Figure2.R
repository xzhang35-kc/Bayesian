
N = 20000
B = 5000
B=0
n = N-B
###############################################################################################################
###############################################################################################################
###################PX-GS#######################################################################################
###############################################################################################################
###############################################################################################################

setwd("C:/ZXiao/MVP-Ordinal-P1/Figure2/PX-GS-SchStudy-PID")

BBeta <- matrix(scan("CorrOrdBeta.dat"),ncol=4,byrow=T)
Beta <- BBeta[(B+1):N,]
b.m <- apply(X = Beta, MARGIN = 2, FUN = mean)
b.m
b.sd <- apply(X = Beta, MARGIN = 2, FUN = sd)
b.sd

RR <- matrix(scan("CorrOrdR.dat"),ncol=16,byrow=T)
R <- RR[(B+1):N,]
R.m <- apply(X = R, MARGIN = 2, FUN = mean)
R.mean = matrix(R.m, 4, 4)
R.mean
R.sd <- apply(X = R, MARGIN = 2, FUN = sd)
R.sd = matrix(R.sd, 4, 4)
R.sd

WWR <- matrix(scan("CorrOrdSig.dat"),ncol=16,byrow=T)
WR <- WWR[(B+1):N,]
WR.m <- apply(X = WR, MARGIN = 2, FUN = mean)
WR.mean = matrix(WR.m, 4, 4)
WR.mean
WR.sd <- apply(X = WR, MARGIN = 2, FUN = sd)
WR.sd = matrix(WR.sd, 4, 4)
WR.sd

GGama <- matrix(scan("CorrOrdGama.dat"),ncol=16,byrow=T)
Gama <- GGama[(B+1):N,]
Gama.m <- apply(X = Gama, MARGIN = 2, FUN = mean)
Gama.sd <- apply(X = Gama, MARGIN = 2, FUN = sd)
G1 = round(Gama.m, 2)
G1
G2 = round(Gama.sd, 2)
G2

rS12 <- as.matrix(R[,2])
rS13 <- as.matrix(R[,3])
rS14 <- as.matrix(R[,4])
rS23 <- as.matrix(R[,7])
rS24 <- as.matrix(R[,8])
rS34 <- as.matrix(R[,12])

bS12 = rep(0,n)
bS13 = rep(0,n)
bS14 = rep(0,n)

bS22 = rep(0,n)
bS23 = rep(0,n)
bS24 = rep(0,n)

bS32 = rep(0,n)
bS33 = rep(0,n)
bS34 = rep(0,n)

bS42 = rep(0,n)
bS43 = rep(0,n)
bS44 = rep(0,n)

gS11 = rep(0,n)
gS12 = rep(0,n)
gS21 = rep(0,n)
gS22 = rep(0,n)
gS31 = rep(0,n)
gS32 = rep(0,n)
gS41 = rep(0,n)
gS42 = rep(0,n)

for(i in 1:n) 
{	
	bS12[i]=Beta[i,2]/sqrt(WR[i,1])
	bS13[i]=Beta[i,3]/sqrt(WR[i,1])
	bS14[i]=Beta[i,4]/sqrt(WR[i,1])
	bS22[i]=Beta[i,2]/sqrt(WR[i,6])
	bS23[i]=Beta[i,3]/sqrt(WR[i,6])
	bS24[i]=Beta[i,4]/sqrt(WR[i,6])
	bS32[i]=Beta[i,2]/sqrt(WR[i,11])
	bS33[i]=Beta[i,3]/sqrt(WR[i,11])
	bS34[i]=Beta[i,4]/sqrt(WR[i,11])
	bS42[i]=Beta[i,2]/sqrt(WR[i,16])
	bS43[i]=Beta[i,3]/sqrt(WR[i,16])
	bS44[i]=Beta[i,4]/sqrt(WR[i,16])
	
	gS11[i]=Gama[i,2]/sqrt(WR[i,1])
	gS12[i]=Gama[i,3]/sqrt(WR[i,1])
	gS21[i]=Gama[i,6]/sqrt(WR[i,6])
	gS22[i]=Gama[i,7]/sqrt(WR[i,6])
	gS31[i]=Gama[i,10]/sqrt(WR[i,11])
	gS32[i]=Gama[i,11]/sqrt(WR[i,11])
	gS41[i]=Gama[i,14]/sqrt(WR[i,16])
	gS42[i]=Gama[i,15]/sqrt(WR[i,16])
}

BetaM = c(mean(bS12), mean(bS13), mean(bS14))
BetaSD = c(sd(bS12), sd(bS13), sd(bS14))
BetaS1Out = round(rbind(BetaM, BetaSD), 2)
BetaS1Out

c1=round(quantile(bS12,probs=c(0.025,0.975)), 2)
c2=round(quantile(bS13,probs=c(0.025,0.975)), 2)
c3=round(quantile(bS14,probs=c(0.025,0.975)), 2)

CIS1Out=c(c1, c2, c3)

BetaM = c(mean(bS22), mean(bS23), mean(bS24))
BetaSD = c(sd(bS22), sd(bS23), sd(bS24))
BetaS2Out = round(rbind(BetaM, BetaSD), 2)
BetaS2Out

c1=round(quantile(bS22,probs=c(0.025,0.975)), 2)
c2=round(quantile(bS23,probs=c(0.025,0.975)), 2)
c3=round(quantile(bS24,probs=c(0.025,0.975)), 2)

CIS2Out=c(c1, c2, c3)


BetaM = c(mean(bS32), mean(bS33), mean(bS34))
BetaSD = c(sd(bS32), sd(bS33), sd(bS34))
BetaS3Out = round(rbind(BetaM, BetaSD), 2)
BetaS3Out

c1=round(quantile(bS32,probs=c(0.025,0.975)), 2)
c2=round(quantile(bS33,probs=c(0.025,0.975)), 2)
c3=round(quantile(bS34,probs=c(0.025,0.975)), 2)

CIS3Out=c(c1, c2, c3)


BetaM = c(mean(bS42), mean(bS43), mean(bS44))
BetaSD = c(sd(bS42), sd(bS43), sd(bS44))
BetaS4Out = round(rbind(BetaM, BetaSD), 2)
BetaS4Out

c1=round(quantile(bS42,probs=c(0.025,0.975)), 2)
c2=round(quantile(bS43,probs=c(0.025,0.975)), 2)
c3=round(quantile(bS44,probs=c(0.025,0.975)), 2)

CIS4Out=c(c1, c2, c3)

GamaM = c(mean(gS11), mean(gS12), mean(gS21), mean(gS22), mean(gS31), mean(gS32), mean(gS41), mean(gS42))
GamaSD = c(sd(gS11), sd(gS12), sd(gS21), sd(gS22), sd(gS31), sd(gS32), sd(gS41), sd(gS42))
GamaSOut = round(rbind(GamaM, GamaSD), 2)
GamaSOut

RSMOut=round(R.mean, 2)
RSMOut
RSSDOut=round(R.sd, 2)
RSSDOut

WRSMOut=round(WR.mean, 2)
WRSMOut
WRSSDOut=round(WR.sd, 2)
WRSSDOut
###############################################################################################################
###############################################################################################################
###################PX-GSM######################################################################################
###############################################################################################################
###############################################################################################################

setwd("C:/ZXiao/MVP-Ordinal-P1/Figure2/PX-GSM-SchStudy-PID")

BBeta <- matrix(scan("CorrOrdBeta.dat"),ncol=4,byrow=T)
Beta <- BBeta[(B+1):N,]
b.m <- apply(X = Beta, MARGIN = 2, FUN = mean)
b.m
b.sd <- apply(X = Beta, MARGIN = 2, FUN = sd)
b.sd

RR <- matrix(scan("CorrOrdR.dat"),ncol=16,byrow=T)
R <- RR[(B+1):N,]
R.m <- apply(X = R, MARGIN = 2, FUN = mean)
R.mean = matrix(R.m, 4, 4)
R.mean
R.sd <- apply(X = R, MARGIN = 2, FUN = sd)
R.sd = matrix(R.sd, 4, 4)
R.sd

WWR <- matrix(scan("CorrOrdSig.dat"),ncol=16,byrow=T)
WR <- WWR[(B+1):N,]
WR.m <- apply(X = WR, MARGIN = 2, FUN = mean)
WR.mean = matrix(WR.m, 4, 4)
WR.mean
WR.sd <- apply(X = WR, MARGIN = 2, FUN = sd)
WR.sd = matrix(WR.sd, 4, 4)
WR.sd

GGama <- matrix(scan("CorrOrdGama.dat"),ncol=16,byrow=T)
Gama <- GGama[(B+1):N,]
Gama.m <- apply(X = Gama, MARGIN = 2, FUN = mean)
Gama.sd <- apply(X = Gama, MARGIN = 2, FUN = sd)
G1 = round(Gama.m, 2)
G1
G2 = round(Gama.sd, 2)
G2

rP12 <- as.matrix(R[,2])
rP13 <- as.matrix(R[,3])
rP14 <- as.matrix(R[,4])
rP23 <- as.matrix(R[,7])
rP24 <- as.matrix(R[,8])
rP34 <- as.matrix(R[,12])


bP12 = rep(0,n)
bP13 = rep(0,n)
bP14 = rep(0,n)

bP22 = rep(0,n)
bP23 = rep(0,n)
bP24 = rep(0,n)

bP32 = rep(0,n)
bP33 = rep(0,n)
bP34 = rep(0,n)

bP42 = rep(0,n)
bP43 = rep(0,n)
bP44 = rep(0,n)

gP11 = rep(0,n)
gP12 = rep(0,n)
gP21 = rep(0,n)
gP22 = rep(0,n)
gP31 = rep(0,n)
gP32 = rep(0,n)
gP41 = rep(0,n)
gP42 = rep(0,n)

for(i in 1:n) 
{
	bP12[i]=Beta[i,2]/sqrt(WR[i,1])
	bP13[i]=Beta[i,3]/sqrt(WR[i,1])
	bP14[i]=Beta[i,4]/sqrt(WR[i,1])
	bP22[i]=Beta[i,2]/sqrt(WR[i,6])
	bP23[i]=Beta[i,3]/sqrt(WR[i,6])
	bP24[i]=Beta[i,4]/sqrt(WR[i,6])
	bP32[i]=Beta[i,2]/sqrt(WR[i,11])
	bP33[i]=Beta[i,3]/sqrt(WR[i,11])
	bP34[i]=Beta[i,4]/sqrt(WR[i,11])
	bP42[i]=Beta[i,2]/sqrt(WR[i,16])
	bP43[i]=Beta[i,3]/sqrt(WR[i,16])
	bP44[i]=Beta[i,4]/sqrt(WR[i,16])

	gP11[i]=Gama[i,2]/sqrt(WR[i,1])
	gP12[i]=Gama[i,3]/sqrt(WR[i,1])
	gP21[i]=Gama[i,6]/sqrt(WR[i,6])
	gP22[i]=Gama[i,7]/sqrt(WR[i,6])
	gP31[i]=Gama[i,10]/sqrt(WR[i,11])
	gP32[i]=Gama[i,11]/sqrt(WR[i,11])
	gP41[i]=Gama[i,14]/sqrt(WR[i,16])
	gP42[i]=Gama[i,15]/sqrt(WR[i,16])
}

BetaM = c(mean(bP12), mean(bP13), mean(bP14))
BetaSD = c(sd(bP12), sd(bP13), sd(bP14))
BetaP1Out = round(rbind(BetaM, BetaSD), 2)
BetaP1Out

c1=round(quantile(bP12,probs=c(0.025,0.975)), 2)
c2=round(quantile(bP13,probs=c(0.025,0.975)), 2)
c3=round(quantile(bP14,probs=c(0.025,0.975)), 2)

CIP1Out=c(c1, c2, c3)

BetaM = c(mean(bP22), mean(bP23), mean(bP24))
BetaSD = c(sd(bP22), sd(bP23), sd(bP24))
BetaP2Out = round(rbind(BetaM, BetaSD), 2)
BetaP2Out

c1=round(quantile(bP22,probs=c(0.025,0.975)), 2)
c2=round(quantile(bP23,probs=c(0.025,0.975)), 2)
c3=round(quantile(bP24,probs=c(0.025,0.975)), 2)

CIP2Out=c(c1, c2, c3)


BetaM = c(mean(bP32), mean(bP33), mean(bP34))
BetaSD = c(sd(bP32), sd(bP33), sd(bP34))
BetaP3Out = round(rbind(BetaM, BetaSD), 2)
BetaP3Out

c1=round(quantile(bP32,probs=c(0.025,0.975)), 2)
c2=round(quantile(bP33,probs=c(0.025,0.975)), 2)
c3=round(quantile(bP34,probs=c(0.025,0.975)), 2)

CIP3Out=c(c1, c2, c3)


BetaM = c(mean(bP42), mean(bP43), mean(bP44))
BetaSD = c(sd(bP42), sd(bP43), sd(bP44))
BetaP4Out = round(rbind(BetaM, BetaSD), 2)
BetaP4Out

c1=round(quantile(bP42,probs=c(0.025,0.975)), 2)
c2=round(quantile(bP43,probs=c(0.025,0.975)), 2)
c3=round(quantile(bP44,probs=c(0.025,0.975)), 2)

CIP4Out=c(c1, c2, c3)

GamaM = c(mean(gP11), mean(gP12), mean(gP21), mean(gP22), mean(gP31), mean(gP32), mean(gP41), mean(gP42))
GamaSD = c(sd(gP11), sd(gP12), sd(gP21), sd(gP22), sd(gP31), sd(gP32), sd(gP41), sd(gP42))
GamaPOut = round(rbind(GamaM, GamaSD), 2)
GamaPOut

RPMOut=round(R.mean, 2)
RPMOut
RPSDOut=round(R.sd, 2)
RPSDOut

WRPMOut=round(WR.mean, 2)
WRPMOut
WRPSDOut=round(WR.sd, 2)
WRPSDOut

###############################################################################################################
###############################################################################################################
###################PX-MH#######################################################################################
###############################################################################################################
###############################################################################################################

setwd("C:/ZXiao/MVP-Ordinal-P1/Figure2/PX-MH-SchStudy-PID")

BBeta <- matrix(scan("CorrOrdBeta.dat"),ncol=4,byrow=T)
Beta <- BBeta[(B+1):N,]
b.m <- apply(X = Beta, MARGIN = 2, FUN = mean)
b.m
b.sd <- apply(X = Beta, MARGIN = 2, FUN = sd)
b.sd

round(quantile(Beta[,2],probs=c(0.025,0.975)), 2)
round(quantile(Beta[,3],probs=c(0.025,0.975)), 2)
round(quantile(Beta[,4],probs=c(0.025,0.975)), 2)


RR <- matrix(scan("CorrOrdR.dat"),ncol=16,byrow=T)
R <- RR[(B+1):N,]
R.m <- apply(X = R, MARGIN = 2, FUN = mean)
R.mean = matrix(R.m, 4, 4)
round(R.mean, 2)
R.sd <- apply(X = R, MARGIN = 2, FUN = sd)
R.sd = matrix(R.sd, 4, 4)
round(R.sd, 2)

WWR <- matrix(scan("CorrOrdWR.dat"),ncol=16,byrow=T)
WR <- WWR[(B+1):N,]
WR.m <- apply(X = WR, MARGIN = 2, FUN = mean)
WR.mean = matrix(WR.m, 4, 4)
WR.mean
WR.sd <- apply(X = WR, MARGIN = 2, FUN = sd)
WR.sd = matrix(WR.sd, 4, 4)
WR.sd

bW1 <- as.matrix(Beta[,1])
bW2 <- as.matrix(Beta[,2])
bW3 <- as.matrix(Beta[,3])
bW4 <- as.matrix(Beta[,4])

rW12 <- as.matrix(R[,2])
rW13 <- as.matrix(R[,3])
rW14 <- as.matrix(R[,4])
rW23 <- as.matrix(R[,7])
rW24 <- as.matrix(R[,8])
rW34 <- as.matrix(R[,12])


GGama <- matrix(scan("CorrOrdGama.dat"),ncol=16,byrow=T)
Gama <- GGama[(B+1):N,]
Gama.m <- apply(X = Gama, MARGIN = 2, FUN = mean)
Gama.sd <- apply(X = Gama, MARGIN = 2, FUN = sd)

G1 = round(Gama.m, 2)
G1
G2 = round(Gama.sd, 2)
G2
	
gW11 = as.matrix(Gama[,2])
gW12 = as.matrix(Gama[,3])
gW21 = as.matrix(Gama[,6])
gW22 = as.matrix(Gama[,7])
gW31 = as.matrix(Gama[,10])
gW32 = as.matrix(Gama[,11])
gW41 = as.matrix(Gama[,14])
gW42 = as.matrix(Gama[,15])

BetaM = c(mean(bW2), mean(bW3), mean(bW4))
BetaSD = c(sd(bW2), sd(bW3), sd(bW4))
BetaWOut = round(rbind(BetaM, BetaSD), 2)
BetaWOut

c1=round(quantile(bW2,probs=c(0.025,0.975)), 2)
c2=round(quantile(bW3,probs=c(0.025,0.975)), 2)
c3=round(quantile(bW4,probs=c(0.025,0.975)), 2)

CIWOut=c(c1, c2, c3)


GamaM = c(mean(gW11), mean(gW12), mean(gW21), mean(gW22),mean(gW31), mean(gW32), mean(gW41), mean(gW42))
GamaSD = c(sd(gW11), sd(gW12), sd(gW21), sd(gW22), sd(gW31), sd(gW32), sd(gW41), sd(gW42))
GamaWOut = round(rbind(GamaM, GamaSD), 2)
GamaWOut

RWMOut=round(R.mean, 2)
RWMOut
RWSDOut=round(R.sd, 2)
RWSDOut

WRWMOut=round(WR.mean, 2)
WRWMOut
WRWSDOut=round(WR.sd, 2)
WRWSDOut



############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
#tiff("Figure_2.tif", width = 4800, height = 4800,res = 800, compression = "lzw")
par(mfrow=c(3,5), mar=c(4,4,1,0.1))


##par(mfrow=c(3,5))
plot(bW3, type = "l", xlab = bquote(beta[3]),ylab = "PX-MH")
plot(rW14, type = "l", xlab = bquote(r[14]), ylab = " ")
plot(rW34, type = "l", xlab = bquote(r[34]),ylab = " ")
plot(gW11, type = "l", xlab = bquote(zeta[11]),ylab = " ")
plot(gW42, type = "l", xlab = bquote(zeta[42]),ylab = " ")


plot(bS13, type = "l", xlab = bquote(beta[3]),ylab = "PX-GS")
plot(rS14, type = "l", xlab = bquote(r[14]), ylab = " ")
plot(rS34, type = "l", xlab = bquote(r[34]),ylab = " ")
plot(gS11, type = "l", xlab = bquote(zeta[11]),ylab = " ")
plot(gS42, type = "l", xlab = bquote(zeta[42]),ylab = " ")

plot(bP13, type = "l", xlab = bquote(beta[3]),ylab = "PX-GSM")
plot(rP14, type = "l", xlab = bquote(r[14]), ylab = " ")
plot(rP34, type = "l", xlab = bquote(r[34]),ylab = " ")
plot(gP11, type = "l", xlab = bquote(zeta[11]),ylab = " ")
plot(gP42, type = "l", xlab = bquote(zeta[42]),ylab = " ")



###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
#dev.off() 
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################


############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################





