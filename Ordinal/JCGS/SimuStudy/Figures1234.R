#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

setwd("C:/ZXiao/MVP-Ordinal/Submission/JCGS/R-prog")

source("simu-funs.R")


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
N <- 500
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


# function for boxplot: only three methods
# number of parameters can be 2, 3, etc.
# for each parameter, boxplot for two methods at the same 
my.box5 <- function(mh, gs, gsm, para.type = c("beta", "gamma", "r", "sigma")) {
   # number of simulations/data sets
   nd <- nrow(gs) 
   
   # create data matrix for plot
   plot.y <- rbind(mh, gs, gsm)
   plot.y <- matrix(plot.y, nrow = nd) # data matrix 

   # number of parameters
   nk <- ncol(plot.y) / 3 
   
   # positions for boxplot
   plot.pos <- c(1:(3 * nk))
   plot.pos <- plot.pos + rep(0:(nk - 1), rep(3, nk))
   
   # colors for boxplot
   plot.col <- rep(c("grey40", "grey70", "grey90"), nk)
   #plot.col <- rep(c("orchid", "slateblue", "gray48"), nk)
   
   # plot without x-axis labels
   boxplot(plot.y, at = plot.pos, col = plot.col, xaxt = "n",
      main = "", sub = " ")
   abline(h=0, col="black")
   
   # add legend
   # legend(x = 1, y = 0.3, fill = c("grey80", "grey90"), legend = c("gs", "gsm"))
   
   # add x-axis labels
   plot.pos <- matrix(plot.pos, nrow = 3)
   plot.pos <- apply(plot.pos, MARGIN = 2, FUN = mean)
   for (k in 1:nk) {
      if (para.type == "beta") {
         plot.names <- bquote(beta[.(k)])
      } else if (para.type == "gamma") {
		a = indexG[k]
         plot.names <- bquote(gamma[.(a)])      
      } else if (para.type == "r") {
		b = indexR[k]
         plot.names <- bquote(r[.(b)])
      }
	  else {
		plot.names <- bquote(sigma[.(k)])
	  }
      axis(side = 1, at = plot.pos[k], labels = plot.names, tick = TRUE)
   }
}

# function for boxplot: only three methods
# number of parameters can be 2, 3, etc.
# for each parameter, boxplot for two methods at the same 
my.box51 <- function(mh, gs, gsm, para.type = c("beta", "gamma", "r", "sigma")) {
   # number of simulations/data sets
   nd <- nrow(gs) 
   
   # create data matrix for plot
   plot.y <- rbind(mh, gs, gsm)
   plot.y <- matrix(plot.y, nrow = nd) # data matrix 

   # number of parameters
   nk <- ncol(plot.y) / 3 
   
   # positions for boxplot
   plot.pos <- c(1:(3 * nk))
   plot.pos <- plot.pos + rep(0:(nk - 1), rep(3, nk))
   
   # colors for boxplot
   plot.col <- rep(c("grey40", "grey70", "grey90"), nk)
   #plot.col <- rep(c("orchid", "slateblue", "gray48"), nk)
   
   # plot without x-axis labels
   boxplot(plot.y, at = plot.pos, col = plot.col, xaxt = "n",
      main = "", sub = " ", ylab = "N=500")
   abline(h=0, col="black")
   
   # add legend
   # legend(x = 1, y = 0.3, fill = c("grey80", "grey90"), legend = c("gs", "gsm"))
   
   # add x-axis labels
   plot.pos <- matrix(plot.pos, nrow = 3)
   plot.pos <- apply(plot.pos, MARGIN = 2, FUN = mean)
   for (k in 1:nk) {
      if (para.type == "beta") {
         plot.names <- bquote(beta[.(k)])
      } else if (para.type == "gamma") {
		a = indexG[k]
         plot.names <- bquote(gamma[.(a)])      
      } else if (para.type == "r") {
		b = indexR[k]
         plot.names <- bquote(r[.(b)])
      }
	  else {
		plot.names <- bquote(sigma[.(k)])
	  }
      axis(side = 1, at = plot.pos[k], labels = plot.names, tick = TRUE)
   }
}
my.box52 <- function(mh, gs, gsm, para.type = c("beta", "gamma", "r", "sigma")) {
   # number of simulations/data sets
   nd <- nrow(gs) 
   
   # create data matrix for plot
   plot.y <- rbind(mh, gs, gsm)
   plot.y <- matrix(plot.y, nrow = nd) # data matrix 

   # number of parameters
   nk <- ncol(plot.y) / 3 
   
   # positions for boxplot
   plot.pos <- c(1:(3 * nk))
   plot.pos <- plot.pos + rep(0:(nk - 1), rep(3, nk))
   
   # colors for boxplot
   plot.col <- rep(c("grey40", "grey70", "grey90"), nk)
   #plot.col <- rep(c("orchid", "slateblue", "gray48"), nk)
   
   # plot without x-axis labels
   boxplot(plot.y, at = plot.pos, col = plot.col, xaxt = "n",
      main = "", sub = " ", ylab = "N=2000")
   abline(h=0, col="black")
   
   # add legend
   # legend(x = 1, y = 0.3, fill = c("grey80", "grey90"), legend = c("gs", "gsm"))
   
   # add x-axis labels
   plot.pos <- matrix(plot.pos, nrow = 3)
   plot.pos <- apply(plot.pos, MARGIN = 2, FUN = mean)
   for (k in 1:nk) {
      if (para.type == "beta") {
         plot.names <- bquote(beta[.(k)])
      } else if (para.type == "gamma") {
		a = indexG[k]
         plot.names <- bquote(gamma[.(a)])      
      } else if (para.type == "r") {
		b = indexR[k]
         plot.names <- bquote(r[.(b)])
      }
	  else {
		plot.names <- bquote(sigma[.(k)])
	  }
      axis(side = 1, at = plot.pos[k], labels = plot.names, tick = TRUE)
   }
}



#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#------------------Load files-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


setwd("C:/ZXiao/MVP-Ordinal/Submission/JCGS/R-prog/SavedResults-500")

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



beta.resmh1 = beta.resmh
beta.resgs1 = beta.resgs
beta.resgsm1 = beta.resgsm

gama.resmh1 = gama.resmh
gama.resgs1 = gama.resgs
gama.resgsm1 = gama.resgsm

r.resmh1 = r.resmh
r.resgs1 = r.resgs
r.resgsm1 = r.resgsm


beta.acfmh1 = beta.acfmh
beta.acfgs1 = beta.acfgs
beta.acfgsm1 = beta.acfgsm

gama.acfmh1 = gama.acfmh
gama.acfgs1 = gama.acfgs
gama.acfgsm1 = gama.acfgsm


r.acfmh1 = r.acfmh
r.acfgs1 = r.acfgs
r.acfgsm1 = r.acfgsm


sigma.acfmh1 = sigma.acfmh
sigma.acfgs1 = sigma.acfgs
sigma.acfgsm1 = sigma.acfgsm


Bmh1 = Bmh
Gmh1 = Gmh
Rmh1 = Rmh
Smh1 = Smh

Bgs1 = Bgs
Ggs1 = Ggs
Rgs1 = Rgs
Sgs1 = Sgs

Bgsm1 = Bgsm
Ggsm1 = Ggsm
Rgsm1 = Rgsm
Sgsm1 = Sgsm



#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#------------------Load files-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------



setwd("C:/ZXiao/MVP-Ordinal/Submission/JCGS/R-prog/SavedResults-2000")


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




beta.resmh2 = beta.resmh
beta.resgs2 = beta.resgs
beta.resgsm2 = beta.resgsm

gama.resmh2 = gama.resmh
gama.resgs2 = gama.resgs
gama.resgsm2 = gama.resgsm

r.resmh2 = r.resmh
r.resgs2 = r.resgs
r.resgsm2 = r.resgsm



beta.acfmh2 = beta.acfmh
beta.acfgs2 = beta.acfgs
beta.acfgsm2 = beta.acfgsm

gama.acfmh2 = gama.acfmh
gama.acfgs2 = gama.acfgs
gama.acfgsm2 = gama.acfgsm


r.acfmh2 = r.acfmh
r.acfgs2 = r.acfgs
r.acfgsm2 = r.acfgsm

sigma.acfmh2 = sigma.acfmh
sigma.acfgs2 = sigma.acfgs
sigma.acfgsm2 = sigma.acfgsm


Bmh2 = Bmh
Gmh2 = Gmh
Rmh2 = Rmh
Smh2 = Smh

Bgs2 = Bgs
Ggs2 = Ggs
Rgs2 = Rgs
Sgs2 = Sgs

Bgsm2 = Bgsm
Ggsm2 = Ggsm
Rgsm2 = Rgsm
Sgsm2 = Sgsm


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------




setwd("C:/ZXiao/MVP-Ordinal/Submission/JCGS/R-prog")



#---------------------------------------------------------------------------------------------------------------------------
#--------------------Figure 1-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

tiff("Figure1.tif", width = 4800, height = 4700,res = 800, compression = "lzw")
par(mfrow=c(3,2), mar=c(4,4,1,1))


#-----------relative bias-----------------------------
a = beta.resmh1$rel.bias
b = beta.resgs1$rel.bias
c = beta.resgsm1$rel.bias
my.box51(mh=a, gs = b, gsm = c, para.type = "beta")


#-----------relative bias-----------------------------
a = beta.resmh2$rel.bias
b = beta.resgs2$rel.bias
c = beta.resgsm2$rel.bias
my.box52(mh=a, gs = b, gsm = c, para.type = "beta")

#-----------relative bias-----------------------------
a = gama.resmh1$rel.bias[,c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
b = gama.resgs1$rel.bias[,c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
c = gama.resgsm1$rel.bias[,c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
my.box51(mh = a, gs = b, gsm = c, para.type = "gamma")


#-----------relative bias-----------------------------
a = gama.resmh2$rel.bias[,c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
b = gama.resgs2$rel.bias[,c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
c = gama.resgsm2$rel.bias[,c(2, 3, 6, 7, 10, 11, 14, 15, 18, 19)]
my.box52(mh = a, gs = b, gsm = c, para.type = "gamma")

#-----------relative bias-----------------------------
a = r.resmh1$rel.bias[,c(2:5, 8:10, 14, 15, 20)]
b = r.resgs1$rel.bias[,c(2:5, 8:10, 14, 15, 20)]
c = r.resgsm1$rel.bias[,c(2:5, 8:10, 14, 15, 20)]
my.box51(mh = a, gs = b, gsm = c, para.type = "r")


#-----------relative bias-----------------------------
a = r.resmh2$rel.bias[,c(2:5, 8:10, 14, 15, 20)]
b = r.resgs2$rel.bias[,c(2:5, 8:10, 14, 15, 20)]
c = r.resgsm2$rel.bias[,c(2:5, 8:10, 14, 15, 20)]
my.box52(mh = a, gs = b, gsm = c, para.type = "r")

dev.off()
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------






#---------------------------------------------------------------------------------------------------------------------------
#---------Figure 2: N=500---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


tiff("Figure2.tif", width = 4500, height = 4800,res = 800, compression = "lzw")
par(mfrow=c(6,5), mar=c(4,3,1.1,0.1))


for (i in 1:P)
{
	plot(beta.acfgsm1[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(beta[~.(i)]), ylab = " ", lty = 1, lwd = 2)
	lines(beta.acfgs1[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	lines(beta.acfmh1[,i], pch = 18, type = "l", lty = 3, lwd = 2)		
}


Gama.acfgsm = gama.acfgsm1[, Lgama]
Gama.acfgs = gama.acfgs1[, Lgama]
Gama.acfmh = gama.acfmh1[, Lgama]
for (i in 1:10)
{
	j=indexG[i]
	plot(Gama.acfgsm[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(gamma[~.(j)]), ylab = " ", lty = 1, lwd = 2)
	lines(Gama.acfgs[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	lines(Gama.acfmh[,i], pch = 18, type = "l", lty = 3, lwd = 2)	
}


R.acfgsm = r.acfgsm1[, Lr]
R.acfgs = r.acfgs1[, Lr]
R.acfmh = r.acfmh1[, Lr]

for (i in 1:10)
{
	j=indexR[i]
	plot(R.acfgsm[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(r[~.(j)]), ylab = " ", lty = 1, lwd = 2)
	lines(R.acfgs[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	lines(R.acfmh[,i], pch = 18, type = "l", lty = 3, lwd = 2)
}


Sig.acfgsm = sigma.acfgsm1[, Lw]
Sig.acfgs = sigma.acfgs1[, Lw]
Sig.acfmh = sigma.acfmh1[, Lw]

for (i in 1:5)
{
	j=indexW[i]
	plot(Sig.acfgsm[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(sigma[~.(j)]), ylab = " ", lty = 1, lwd = 2)
	lines(Sig.acfgs[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	#lines(Sig.acfmh[,i], pch = 18, type = "l", lty = 3, lwd = 2)
}

dev.off()

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------
#---------Figure 3: N=2000--------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------




tiff("Figure3.tif", width = 4500, height = 4800,res = 800, compression = "lzw")
par(mfrow=c(6,5), mar=c(4,3,1.1,0.1))


for (i in 1:P)
{
	plot(beta.acfgsm2[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(beta[~.(i)]), ylab = " ", lty = 1, lwd = 2)
	lines(beta.acfgs2[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	lines(beta.acfmh2[,i], pch = 18, type = "l", lty = 3, lwd = 2)		
}


Gama.acfgsm = gama.acfgsm2[, Lgama]
Gama.acfgs = gama.acfgs2[, Lgama]
Gama.acfmh = gama.acfmh2[, Lgama]
for (i in 1:10)
{
	j=indexG[i]
	plot(Gama.acfgsm[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(gamma[~.(j)]), ylab = " ", lty = 1, lwd = 2)
	lines(Gama.acfgs[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	lines(Gama.acfmh[,i], pch = 18, type = "l", lty = 3, lwd = 2)	
}


R.acfgsm = r.acfgsm2[, Lr]
R.acfgs = r.acfgs2[, Lr]
R.acfmh = r.acfmh2[, Lr]

for (i in 1:10)
{
	j=indexR[i]
	plot(R.acfgsm[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(r[~.(j)]), ylab = " ", lty = 1, lwd = 2)
	lines(R.acfgs[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	lines(R.acfmh[,i], pch = 18, type = "l", lty = 3, lwd = 2)
}


Sig.acfgsm = sigma.acfgsm2[, Lw]
Sig.acfgs = sigma.acfgs2[, Lw]
Sig.acfmh = sigma.acfmh2[, Lw]

for (i in 1:5)
{
	j=indexW[i]
	plot(Sig.acfgsm[,i], type = "l", frame = FALSE, pch = 19, xlab = bquote(sigma[~.(j)]), ylab = " ", lty = 1, lwd = 2)
	lines(Sig.acfgs[,i], pch = 18, type = "l", lty = 5, lwd = 2)
	#lines(Sig.acfmh[,i], pch = 18, type = "l", lty = 3, lwd = 2)
}

dev.off()

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------
#--------------------Figure 4-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

tiff("Figure4.tif", width = 4800, height = 4500,res = 800, compression = "lzw")
par(mfrow=c(4,5), mar=c(4,4,1,1))


plot(Sgs1[,1], type="l", ylim=c(0,3.5),  col = "black",  xlab = bquote(sigma[~.(11)]), ylab="N=500")
plot(Sgs1[,7], type="l", ylim=c(0,3.5),  col = "black",  xlab = bquote(sigma[~.(22)]), ylab=" ")
plot(Sgs1[,13], type="l", ylim=c(0,3.5),  col = "black",  xlab = bquote(sigma[~.(33)]), ylab=" ")
plot(Sgs1[,19], type="l", ylim=c(0,3.5),  col = "black",  xlab = bquote(sigma[~.(44)]), ylab=" ")
plot(Sgs1[,25], type="l", ylim=c(0,2),  col = "black",  xlab = bquote(sigma[~.(55)]), ylab=" ")


plot(Sgs2[,1], type="l", ylim=c(0,1),  col = "black",  xlab = bquote(sigma[~.(11)]), ylab=" N=2000")
plot(Sgs2[,7], type="l", ylim=c(0,1),  col = "black",  xlab = bquote(sigma[~.(22)]), ylab=" ")
plot(Sgs2[,13], type="l", ylim=c(0,0.5),  col = "black",  xlab = bquote(sigma[~.(33)]), ylab=" ")
plot(Sgs2[,19], type="l", ylim=c(0,0.4),  col = "black",  xlab = bquote(sigma[~.(44)]), ylab=" ")
plot(Sgs2[,25], type="l", ylim=c(0,0.3),  col = "black",  xlab = bquote(sigma[~.(55)]), ylab=" ")



plot(Sgsm1[,1], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(11)]), ylab="N=500")
plot(Sgsm1[,7], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(22)]), ylab=" ")
plot(Sgsm1[,13], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(33)]), ylab=" ")
plot(Sgsm1[,19], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(44)]), ylab=" ")
plot(Sgsm1[,25], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(55)]), ylab=" ")


plot(Sgsm2[,1], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(11)]), ylab=" N=2000")
plot(Sgsm2[,7], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(22)]), ylab=" ")
plot(Sgsm2[,13], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(33)]), ylab=" ")
plot(Sgsm2[,19], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(44)]), ylab=" ")
plot(Sgsm2[,25], type="l", ylim=c(0,6),  col = "black",  xlab = bquote(sigma[~.(55)]), ylab=" ")


dev.off()

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
