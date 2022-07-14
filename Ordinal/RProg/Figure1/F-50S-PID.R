
NumIte=10000
NumSamp=50
nn=41

#####################################################################################################################
################SGDA#################################################################################################
#####################################################################################################################
NumSamp=50

setwd("C:/ZXiao/MVP-Ordinal-P1/Figure1/PX-GS-Simu-PID/Simu50S5D-UX")
Beta <- matrix(scan("CorrOrdBeta.dat"),ncol=2,byrow=T)
R <- matrix(scan("CorrOrdR.dat"),ncol=25,byrow=T)
WR <- matrix(scan("CorrOrdSig.dat"),ncol=25,byrow=T)
Gama <- matrix(scan("CorrOrdGama.dat"), ncol=20, byrow=T)

sa1 = rep(0, nn)
sa2 = rep(0, nn)
sr12 = rep(0, nn)
sr13 = rep(0, nn)
sr14 = rep(0, nn)
sr15 = rep(0, nn)
sr23 = rep(0, nn)
sr24 = rep(0, nn)
sr25 = rep(0, nn)
sr34 = rep(0, nn)
sr35 = rep(0, nn)
sr45 = rep(0, nn)

sg11 = rep(0, nn)
sg12 = rep(0, nn)
sg21 = rep(0, nn)
sg22 = rep(0, nn)
sg31 = rep(0, nn)
sg32 = rep(0, nn)
sg41 = rep(0, nn)
sg42 = rep(0, nn)
sg51 = rep(0, nn)
sg52 = rep(0, nn)

b0=rep(0, NumIte)
b1=rep(0, NumIte)
g11 = rep(0, NumIte)
g12 = rep(0, NumIte)
g21 = rep(0, NumIte)
g22 = rep(0, NumIte)
g31 = rep(0, NumIte)
g32 = rep(0, NumIte)
g41 = rep(0, NumIte)
g42 = rep(0, NumIte)
g51 = rep(0, NumIte)
g52 = rep(0, NumIte)


for(i in 1:NumSamp)
{
	n1=(i-1)*NumIte+1
	n2=i*NumIte
	BBeta = Beta[n1:n2,]
	GGama = Gama[n1:n2,]
	WWR = WR[n1:n2,]
	for(j in 1:NumIte)
	{
		b0[j]=BBeta[j,1]/sqrt(WWR[j,1])
		b1[j]=BBeta[j,2]/sqrt(WWR[j,1])
		g11[j]=GGama[j,2]/sqrt(WWR[j,1])
		g12[j]=GGama[j,3]/sqrt(WWR[j,1])
		g21[j]=GGama[j,6]/sqrt(WWR[j,7])
		g22[j]=GGama[j,7]/sqrt(WWR[j,7])
		g31[j]=GGama[j,10]/sqrt(WWR[j,13])
		g32[j]=GGama[j,11]/sqrt(WWR[j,13])
		g41[j]=GGama[j,14]/sqrt(WWR[j,19])
		g42[j]=GGama[j,15]/sqrt(WWR[j,19])
		g51[j]=GGama[j,18]/sqrt(WWR[j,25])
		g52[j]=GGama[j,19]/sqrt(WWR[j,25])	
	}	
	a1=acf(b0, plot = FALSE)
	a2=acf(b1, plot = FALSE)
	sa1=sa1+a1$acf
	sa2=sa2+a2$acf
	g1=acf(g11, plot = FALSE)
	g2=acf(g12, plot = FALSE)
	sg11=sg11+g1$acf
	sg12=sg12+g2$acf
	g3=acf(g21, plot = FALSE)
	g4=acf(g22, plot = FALSE)
	sg21=sg21+g3$acf
	sg22=sg22+g4$acf
	g5=acf(g31, plot = FALSE)
	g6=acf(g32, plot = FALSE)
	sg31=sg31+g5$acf
	sg32=sg32+g6$acf
	g7=acf(g41, plot = FALSE)
	g8=acf(g42, plot = FALSE)
	sg41=sg41+g7$acf
	sg42=sg42+g8$acf
	g9=acf(g51, plot = FALSE)
	g10=acf(g52, plot = FALSE)
	sg51=sg51+g9$acf
	sg52=sg52+g10$acf	
	RR = R[n1:n2,]
	r12=acf(RR[,2], plot = FALSE)
	sr12=sr12+r12$acf
	r13=acf(RR[,3], plot = FALSE)
	sr13=sr13+r13$acf
	r14=acf(RR[,4], plot = FALSE)
	sr14=sr14+r14$acf
	r15=acf(RR[,5], plot = FALSE)
	sr15=sr15+r15$acf
	r23=acf(RR[,8], plot = FALSE)
	sr23=sr23+r23$acf
	r24=acf(RR[,9], plot = FALSE)
	sr24=sr24+r24$acf
	r25=acf(RR[,10], plot = FALSE)
	sr25=sr25+r25$acf
	r34=acf(RR[,14], plot = FALSE)
	sr34=sr14+r14$acf
	r35=acf(RR[,15], plot = FALSE)
	sr35=sr35+r35$acf
	r45=acf(RR[,20], plot = FALSE)
	sr45=sr45+r45$acf	
}
acfsb0 = sa1/NumSamp
acfsb1 = sa2/NumSamp
acfsr12 = sr12/NumSamp
acfsr13 = sr13/NumSamp
acfsr14 = sr14/NumSamp
acfsr15 = sr15/NumSamp
acfsr23 = sr23/NumSamp
acfsr24 = sr24/NumSamp
acfsr25 = sr25/NumSamp
acfsr34 = sr34/NumSamp
acfsr35 = sr35/NumSamp
acfsr45 = sr45/NumSamp
acfsg11 = sg11/NumSamp
acfsg12 = sg12/NumSamp
acfsg21 = sg21/NumSamp
acfsg22 = sg22/NumSamp
acfsg31 = sg31/NumSamp
acfsg32 = sg32/NumSamp
acfsg41 = sg41/NumSamp
acfsg42 = sg42/NumSamp
acfsg51 = sg51/NumSamp
acfsg52 = sg52/NumSamp

#####################################################################################################################
################PMDA#################################################################################################
#####################################################################################################################
setwd("C:/ZXiao/MVP-Ordinal-P1/Figure1/PX-GSM-Simu-PID/Simu50S5D-UX")
Beta <- matrix(scan("CorrOrdBeta.dat"),ncol=2,byrow=T)
R <- matrix(scan("CorrOrdR.dat"),ncol=25,byrow=T)
WR <- matrix(scan("CorrOrdSig.dat"),ncol=25,byrow=T)
Gama <- matrix(scan("CorrOrdGama.dat"), ncol=20, byrow=T)

sa1 = rep(0, nn)
sa2 = rep(0, nn)
sr12 = rep(0, nn)
sr13 = rep(0, nn)
sr14 = rep(0, nn)
sr15 = rep(0, nn)
sr23 = rep(0, nn)
sr24 = rep(0, nn)
sr25 = rep(0, nn)
sr34 = rep(0, nn)
sr35 = rep(0, nn)
sr45 = rep(0, nn)

sg11 = rep(0, nn)
sg12 = rep(0, nn)
sg21 = rep(0, nn)
sg22 = rep(0, nn)
sg31 = rep(0, nn)
sg32 = rep(0, nn)
sg41 = rep(0, nn)
sg42 = rep(0, nn)
sg51 = rep(0, nn)
sg52 = rep(0, nn)

b0=rep(0, NumIte)
b1=rep(0, NumIte)
g11 = rep(0, NumIte)
g12 = rep(0, NumIte)
g21 = rep(0, NumIte)
g22 = rep(0, NumIte)
g31 = rep(0, NumIte)
g32 = rep(0, NumIte)
g41 = rep(0, NumIte)
g42 = rep(0, NumIte)
g51 = rep(0, NumIte)
g52 = rep(0, NumIte)

for(i in 1:NumSamp)
{
	n1=(i-1)*NumIte+1
	n2=i*NumIte
	BBeta = Beta[n1:n2,]
	GGama = Gama[n1:n2,]
	WWR = WR[n1:n2,]
	for(j in 1:NumIte)
	{
		b0[j]=BBeta[j,1]/sqrt(WWR[j,1])
		b1[j]=BBeta[j,2]/sqrt(WWR[j,1])
		g11[j]=GGama[j,2]/sqrt(WWR[j,1])
		g12[j]=GGama[j,3]/sqrt(WWR[j,1])
		g21[j]=GGama[j,6]/sqrt(WWR[j,7])
		g22[j]=GGama[j,7]/sqrt(WWR[j,7])
		g31[j]=GGama[j,10]/sqrt(WWR[j,13])
		g32[j]=GGama[j,11]/sqrt(WWR[j,13])
		g41[j]=GGama[j,14]/sqrt(WWR[j,19])
		g42[j]=GGama[j,15]/sqrt(WWR[j,19])
		g51[j]=GGama[j,18]/sqrt(WWR[j,25])
		g52[j]=GGama[j,19]/sqrt(WWR[j,25])	
	}	
	a1=acf(b0, plot = FALSE)
	a2=acf(b1, plot = FALSE)
	sa1=sa1+a1$acf
	sa2=sa2+a2$acf	
	g1=acf(g11, plot = FALSE)
	g2=acf(g12, plot = FALSE)
	sg11=sg11+g1$acf
	sg12=sg12+g2$acf
	g3=acf(g21, plot = FALSE)
	g4=acf(g22, plot = FALSE)
	sg21=sg21+g3$acf
	sg22=sg22+g4$acf
	g5=acf(g31, plot = FALSE)
	g6=acf(g32, plot = FALSE)
	sg31=sg31+g5$acf
	sg32=sg32+g6$acf
	g7=acf(g41, plot = FALSE)
	g8=acf(g42, plot = FALSE)
	sg41=sg41+g7$acf
	sg42=sg42+g8$acf
	g9=acf(g51, plot = FALSE)
	g10=acf(g52, plot = FALSE)
	sg51=sg51+g9$acf
	sg52=sg52+g10$acf	
	RR = R[n1:n2,]
	r12=acf(RR[,2], plot = FALSE)
	sr12=sr12+r12$acf
	r13=acf(RR[,3], plot = FALSE)
	sr13=sr13+r13$acf
	r14=acf(RR[,4], plot = FALSE)
	sr14=sr14+r14$acf
	r15=acf(RR[,5], plot = FALSE)
	sr15=sr15+r15$acf
	r23=acf(RR[,8], plot = FALSE)
	sr23=sr23+r23$acf
	r24=acf(RR[,9], plot = FALSE)
	sr24=sr24+r24$acf
	r25=acf(RR[,10], plot = FALSE)
	sr25=sr25+r25$acf
	r34=acf(RR[,14], plot = FALSE)
	sr34=sr14+r14$acf
	r35=acf(RR[,15], plot = FALSE)
	sr35=sr35+r35$acf
	r45=acf(RR[,20], plot = FALSE)
	sr45=sr45+r45$acf	
}

acfpb0 = sa1/NumSamp
acfpb1 = sa2/NumSamp
acfpr12 = sr12/NumSamp
acfpr13 = sr13/NumSamp
acfpr14 = sr14/NumSamp
acfpr15 = sr15/NumSamp
acfpr23 = sr23/NumSamp
acfpr24 = sr24/NumSamp
acfpr25 = sr25/NumSamp
acfpr34 = sr34/NumSamp
acfpr35 = sr35/NumSamp
acfpr45 = sr45/NumSamp
acfpg11 = sg11/NumSamp
acfpg12 = sg12/NumSamp
acfpg21 = sg21/NumSamp
acfpg22 = sg22/NumSamp
acfpg31 = sg31/NumSamp
acfpg32 = sg32/NumSamp
acfpg41 = sg41/NumSamp
acfpg42 = sg42/NumSamp
acfpg51 = sg51/NumSamp
acfpg52 = sg52/NumSamp

#####################################################################################################################
################W-MH#################################################################################################
#####################################################################################################################
NumSamp=50

setwd("C:/ZXiao/MVP-Ordinal-P1/Figure1/PX-MH-Simu-PID/Simu50S5D-UX")
Beta <- matrix(scan("CorrOrdBeta.dat"),ncol=2,byrow=T)
R <- matrix(scan("CorrOrdR.dat"),ncol=25,byrow=T)
WR <- matrix(scan("CorrOrdWR.dat"),ncol=25,byrow=T)
Gama <- matrix(scan("CorrOrdGama.dat"), ncol=20, byrow=T)

sa1 = rep(0, nn)
sa2 = rep(0, nn)
sr12 = rep(0, nn)
sr13 = rep(0, nn)
sr14 = rep(0, nn)
sr15 = rep(0, nn)
sr23 = rep(0, nn)
sr24 = rep(0, nn)
sr25 = rep(0, nn)
sr34 = rep(0, nn)
sr35 = rep(0, nn)
sr45 = rep(0, nn)

sg11 = rep(0, nn)
sg12 = rep(0, nn)
sg21 = rep(0, nn)
sg22 = rep(0, nn)
sg31 = rep(0, nn)
sg32 = rep(0, nn)
sg41 = rep(0, nn)
sg42 = rep(0, nn)
sg51 = rep(0, nn)
sg52 = rep(0, nn)

for(i in 1:NumSamp)
{
	n1=(i-1)*NumIte+1
	n2=i*NumIte
	BBeta = Beta[n1:n2,]
	a1=acf(BBeta[,1], plot = FALSE)
	a2=acf(BBeta[,2], plot = FALSE)
	sa1=sa1+a1$acf
	sa2=sa2+a2$acf
	RR = R[n1:n2,]
	r12=acf(RR[,2], plot = FALSE)
	sr12=sr12+r12$acf
	r13=acf(RR[,3], plot = FALSE)
	sr13=sr13+r13$acf
	r14=acf(RR[,4], plot = FALSE)
	sr14=sr14+r14$acf
	r15=acf(RR[,5], plot = FALSE)
	sr15=sr15+r15$acf
	r23=acf(RR[,8], plot = FALSE)
	sr23=sr23+r23$acf
	r24=acf(RR[,9], plot = FALSE)
	sr24=sr24+r24$acf
	r25=acf(RR[,10], plot = FALSE)
	sr25=sr25+r25$acf
	r34=acf(RR[,14], plot = FALSE)
	sr34=sr14+r14$acf
	r35=acf(RR[,15], plot = FALSE)
	sr35=sr35+r35$acf
	r45=acf(RR[,20], plot = FALSE)
	sr45=sr45+r45$acf	
	GGama = Gama[n1:n2,]
	g11=acf(GGama[,2], plot = FALSE)
	sg11=sg11+g11$acf
	g12=acf(GGama[,3], plot = FALSE)
	sg12=sg12+g12$acf
	g21=acf(GGama[,6], plot = FALSE)
	sg21=sg11+g21$acf
	g22=acf(GGama[,7], plot = FALSE)
	sg22=sg22+g22$acf
	g31=acf(GGama[,10], plot = FALSE)
	sg31=sg11+g31$acf
	g32=acf(GGama[,11], plot = FALSE)
	sg32=sg32+g32$acf
	g41=acf(GGama[,14], plot = FALSE)
	sg41=sg41+g41$acf
	g42=acf(GGama[,15], plot = FALSE)
	sg42=sg42+g42$acf
	g51=acf(GGama[,18], plot = FALSE)
	sg51=sg51+g51$acf
	g52=acf(GGama[,19], plot = FALSE)
	sg52=sg52+g52$acf
}

acfwb0 = sa1/NumSamp
acfwb1 = sa2/NumSamp
acfwr12 = sr12/NumSamp
acfwr13 = sr13/NumSamp
acfwr14 = sr14/NumSamp
acfwr15 = sr15/NumSamp
acfwr23 = sr23/NumSamp
acfwr24 = sr24/NumSamp
acfwr25 = sr25/NumSamp
acfwr34 = sr34/NumSamp
acfwr35 = sr35/NumSamp
acfwr45 = sr45/NumSamp
acfwg11 = sg11/NumSamp
acfwg12 = sg12/NumSamp
acfwg21 = sg21/NumSamp
acfwg22 = sg22/NumSamp
acfwg31 = sg31/NumSamp
acfwg32 = sg32/NumSamp
acfwg41 = sg41/NumSamp
acfwg42 = sg42/NumSamp
acfwg51 = sg51/NumSamp
acfwg52 = sg52/NumSamp


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################


pb150ID=acfpb1
sb150ID=acfsb1
wb150ID=acfwb1

pg1150ID=acfpg11
sg1150ID=acfsg11
wg1150ID=acfwg11

pr2350ID=acfpr23
sr2350ID=acfsr23
wr2350ID=acfwr23

pr1550ID=acfpr15
sr1550ID=acfsr15
wr1550ID=acfwr15

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

