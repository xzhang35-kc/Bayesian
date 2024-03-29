



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
Output="PID-B5000-"
Plot="PID-B5000-"


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







#----------------------------------------------------------------------
#----------------------------------------------------------------------
#---PX-MH--------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------

m.id <- "PX-MH" 
file.pre <- paste(wd, "res-pid-px-mh/", "res", sep = "")
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
#----------------------------------------------------------------------
#----------------------------------------------------------------------
#----------------------------------------------------------------------




#----------------------------------------------------------------------
#----------Begin Diagnostics PX-MH-------------------------------------
#----------------------------------------------------------------------

beta.m <- matrix(beta.list[[1]], ncol=P, byrow=T)
head(beta.m)
apply(beta.m, MARGIN=2, FUN=mean)

beta.m = beta.m[(burn.in+1):G,]
apply(beta.m, MARGIN=2, FUN=mean)
Dbeta.mh = beta.m
#----------------------------------------------------------------------
#----------Gama PX-MH--------------------------------------------------
#----------------------------------------------------------------------
gama.m <- matrix(gama.list[[1]], ncol=K*CP, byrow=T)
head(gama.m)
apply(gama.m, MARGIN=2, FUN=mean)

gama.m = as.matrix(gama.m[(burn.in+1):G,Lgama])
apply(gama.m, MARGIN=2, FUN=mean)
Dgama.mh = gama.m

#----------------------------------------------------------------------
#----------R PX-MH-----------------------------------------------------
#----------------------------------------------------------------------
r.m <- matrix(r.list[[1]], ncol=K*K, byrow=T)
head(r.m)
apply(r.m, MARGIN=2, FUN=mean)

r.m = r.m[(burn.in+1):G,Lr]
apply(r.m, MARGIN=2, FUN=mean)
Dr.mh = r.m
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
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
      main = "", sub = " ")
  # abline(h=0, col="black")
   
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



#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------------------------
#--------------------Figure 6-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
tiff("Figure6.tif", width = 4800, height = 4500,res = 800, compression = "lzw")
my.box51(mh = Dgama.mh, gs = Dgama.gs, gsm = Dgama.gsm, para.type = "gamma")


dev.off()
#---------------------------------------------------------------------------------------------------------------------------
#--------------------Figure 7-----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

tiff("Figure7.tif", width = 4800, height = 4500,res = 800, compression = "lzw")
my.box51(mh = Dr.mh, gs = Dr.gs, gsm = Dr.gsm, para.type = "r")


dev.off()
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------




