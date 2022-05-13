#----------------------------------------------------------------------
# functions used in simulation studies
#----------------------------------------------------------------------

# load packages
library(mvtnorm)
library(boa)
library(coda)
library(mcmcplots)


# function to simulate the data without missing
#   X - data matrix, supplied by users
#   theta - regression coefficient, supplied by users
#      it is a length K list
#   C - covariance matrix, supplied by users
#   K - number of repeats, supplied by users   
#   N - sample size, default is 500
#   cuts - cutoff points for ordinal variables, a list
# return: a list
simu.ordinal <- function(X, theta, C, K, N = 500, cuts)
{
   # check data consistency
   if (ncol(X) != length(theta) || K != nrow(C) || K != ncol(C) || K  != length(cuts)
      || nrow(X) < K * N)
      stop("data input error in simu.ordinal")
   
   # randomly select n samples from original data
   n.s <- nrow(X) / K # number of sample
   P <- ncol(X) # number of covariates
   if (n.s < N)
      stop("the orginal sample size is too small")
   bn <- sample(1:n.s, N)
   bX <- matrix(t(X), ncol = n.s)
   bX <- bX[, bn]
   bX <- matrix(bX, ncol = P, byrow = TRUE)
   
   # create latent variables
   Z <- matrix(NA, N, K)
   for (n in 1:N) 
   {
      t1 <- (n - 1) * K + 1
      t2 <- (n - 1) * K + K
      mu <- bX[t1:t2, ] %*% theta
      Z[n,] <- rmvnorm(1, mu, C)
   }
   
   # create responses
   Y <- as.data.frame(matrix(NA, N, K))
   for (k in 1:K) {
      Y[, k] <- cut(Z[, k], breaks = cuts[[k]], ordered_result = TRUE)
      Y[, k] <- as.numeric(Y[, k])
   }
   Y <- Y - 1 # ordered data is from 0
   
   # for output
   Z <- as.data.frame(Z)
   Y <- as.data.frame(Y)
   names(Z) <- paste("r", 1:K, sep = "")
   names(Y) <- paste("r", 1:K, sep = "")
   # vectorized results
   Z2 <- as.vector(t(Z))
   Y2 <- as.vector(t(Y))
   
   # type of responses
   n <- sapply(cuts, length)
   rType <- rep("Ord", K)
   rType[n <= 2] <- "Con"
   
   compl.dat <- list(Y = Y, X = bX, Z = Z, Y2 = Y2, Z2 = Z2, OX = X, rType = rType, 
      bn = bn, C = C, cuts = cuts, theta = theta, N = N, K = K, P = P)
   
   return(compl.dat)
}

# function for simulating dropouts and missing
#   compl.dat - list from function simu.ordinal
#      list contains the complete data    
#   reg.drop - logistic regression coefficients for dropouts
#   prop.drop - dropping rates for each year
#   reg.miss - logistic regression coefficients for missing data 
#   prop.miss - missing rates for each year after removing dropouts
# return: a list
add.missing <- function(compl.dat, reg.drop, prop.drop,
   reg.miss, prop.miss) {
   N <- compl.dat$N # number of samples
   K <- compl.dat$K # number of repeats 
   P <- compl.dat$P # number of covaiates
   
   # check data consistency
   if (length(reg.drop) != P || length(reg.miss) != P || 
      length(prop.drop) != K || length(prop.miss) != K)
      stop("data input errors in add.drop.missing")
   
   # first simulate dropouts
   Y <- compl.dat$Y
   X <- compl.dat$X
   year.dropped <- rep(NA, N)
   years.dropped <- rep(0, N)
   for (n in 1:N) {
      p.drop <- X[n * K - K + 1, ] * reg.drop
      p.drop <- sum(p.drop)
      p.drop <- exp(p.drop) / (1 + exp(p.drop))
      if.drop <- sample(c(FALSE, TRUE), size = 1, prob = c(1 - p.drop, p.drop))
      if (if.drop) {
         year.dropped[n] <- sample(1:K, size = 1, prob = prop.drop)
         years.dropped[n] <-  K -  year.dropped[n] + 1
         Y[n, year.dropped[n]:K] <- 999 
      }
   }
   
   # second simulate missing data
   years.missed <- rep(0, N)
   for (k in 3:(K - 1)) {
      for (n in 1:n) {
         if (years.dropped[n] >= (K - k))
            next
         p.miss <- X[n * K - K + 1, ] * reg.miss
         p.miss <- sum(p.miss)
         p.miss <- exp(p.miss) / (1 + exp(p.miss))
         if.miss <- sample(c(FALSE, TRUE), size = 1, prob = c(1 - p.miss, p.miss))
         if (if.miss) {
            years.missed[n] = years.missed[n] + 1
            Y[n, k] <- 999
         }
      }
   }
   
   # vectorize Y
   Y2 <- as.vector(t(Y))
   
   drop.dat <- list(Y = Y, Y2 = Y2, X = X, year.dropped = year.dropped, years.dropped = years.dropped,
      years.missed = years.missed, reg.drop = reg.drop, prop.drop = prop.drop,   
      reg.miss = reg.miss, prop.miss = prop.miss)
   
   return(list(compl.dat = compl.dat, drop.dat = drop.dat))
}

# function for finding variances
#   sigma.vec: covariance matrix
#   G: number of Gibbs
#   K: number of repeats
find.var <- function(sigma.vec, G, K) {
   var.index <-  (1:K) * K - K + 1:K
   var.index <- expand.grid((1:G) * K * K - K * K, var.index)
   var.index <- apply(var.index, FUN = sum, MARGIN = 1)
   var.index <- sort(var.index)
   var.vec <- sigma.vec[var.index]   
   return(var.vec) 
}

# function for adjusting beta
#   beta.vec: regression coefficients
#   sd.vec: standard deviations
#   G: number of Gibbs
#   K: number of repeats
beta.adj <- function(beta.vec, sd.vec, G, K) {
   beta.m <- matrix(beta.vec, nrow = G, byrow = T)
   sd.m <- sd.vec[(1:G)* K - K + 1]
   beta.m <- beta.m / sd.m
   return(as.vector(t(beta.m)))   
}
beta.adj1 <- function(beta.vec, sd.vec, G, K, P) {
   if(length(beta.vec) != G * P | length(sd.vec) != G * K) {
      stop("check length of beta or sd, and G, K, or P!\n"); 
   }
   # create matrix for results
   res.m <- matrix(NA, nrow = K, ncol = G * P)
   # adjustment
   for(k in 1:K) {
      sd.k <- sd.vec[(1:G) * K - K + k]
      sd.k <- rep(sd.k, times = rep(P, times = G))     
      res.m[k, ] <- beta.vec / sd.k
   }
   return(as.vector(res.m))   
}

beta.adj2 <- function(beta.vec, sd.vec, G, K) {
   beta.m <- beta.vec / sd.vec
   return(beta.m)   
}

# function for adjusting gama
#   gama.vec: cut-off points
#   sd.vec: standard deviations
#   G: number of Gibbs
#   K: number of repeats
gama.adj <- function(gama.vec, sd.vec, G, K) {
   gama.m <- matrix(gama.vec, nrow = G, byrow = T)
   NC <- ncol(gama.m) / K
   sd.m <- rep(sd.vec, rep(NC, length(sd.vec)))
   sd.m <- matrix(sd.m, nrow = G, byrow = T)
   gama.m <- gama.m / sd.m
   return(as.vector(t(gama.m)))
}

# function for adjusting gama
#   gama.vec: cut-off points
#   sd.vec: standard deviations
#   G: number of Gibbs
#   K: number of repeats
#   NC: number of cut points
gama.adj <- function(gama.vec, sd.vec, G, K, NC) {
   gama.m <- matrix(gama.vec, nrow = G, byrow = T)
   sd.m <- rep(sd.vec, rep(NC, G))
   sd.m <- matrix(sd.m, nrow = G, byrow = T)
   gama.m <- gama.m / sd.m
   return(as.vector(t(gama.m)))
}

# obtain summary statisic for parameters
#   para.vec: parameter vector
#   G: number of repeats
para.stat <- function(para.vec, G, burn.in = 0) {
   para.m <- matrix(para.vec, nrow = G, byrow = T)
   if (burn.in >= G)
      stop("burn in is too large")
   if (burn.in >= 1) {      
      para.m <- para.m[-(1:burn.in), ]
   }
   para.mean <- apply(para.m, FUN = mean, MARGIN = 2)
   para.sd <- apply(para.m, FUN = sd, MARGIN = 2)
   c(para.mean, para.sd)
} 

# obtain acf for parameters
#   para.vec: parameter vector
#   G: number of repeats
#   burn.in: number of burn in
my.acf <- function(x) {
   res <- acf(x, plot = FALSE)
   res$acf
}
para.vec.acf <- function(para.vec, G, burn.in = 0) {
   para.m <- matrix(para.vec, nrow = G, byrow = T)
   if (burn.in >= G)
      stop("burn in is too large")
   if (burn.in >= 1) {      
      para.m <- para.m[-(1:burn.in), ]
   }
   res.acf <- apply(para.m, FUN = my.acf, MARGIN = 2)
   return(res.acf)
}
para.list.acf <- function(para.list, G, burn.in = 0) {
   S <- length(para.list)
   if (S <= 0) {
      stop("list is empty in para.list.acf")
   }
   res.acf <- para.vec.acf(para.list[[1]], G = G, burn.in = burn.in)
   if (S >= 2) {
      for (s in 2:S) {
         res.acf <- res.acf + para.vec.acf(para.list[[s]], G = G, burn.in = burn.in)
      }
   }   
   res.acf <- res.acf / S
   return(res.acf)
}

# obtain CI for parameters
#   para.vec: parameter vector
#   G: number of repeats
para.CI <- function(para.vec, G, burn.in = 0, alpha = 0.05) {
   para.m <- matrix(para.vec, nrow = G, byrow = T)
   if (burn.in >= G)
      stop("burn in is too large")
   if (burn.in >= 1) {      
      para.m <- para.m[-(1:burn.in), ]
   }
   ci.lower <- apply(para.m, FUN = quantile, MARGIN = 2, probs = alpha / 2);
   ci.upper <- apply(para.m, FUN = quantile, MARGIN = 2, probs = 1 - alpha / 2);
   res <- list(ci.lower = ci.lower, ci.upper = ci.upper)
   return(res)
} 

# obtain coverage for parameters
#   para.vec: parameter vector
#   true.vec: vector for true parameters 
#   G: number of repeats
para.CI.cov <- function(para.vec, true.vec, G, burn.in = 0, alpha = 0.05) {
   ci.res <- para.CI(para.vec, G = G, burn.in = burn.in, alpha = alpha)
   ci.lower <- ci.res$ci.lower
   ci.upper <- ci.res$ci.upper
   ci.cov <- (ci.lower < true.vec) & (ci.upper > true.vec)
   return(ci.cov)
} 

# obtain relative bias, coverage of 95% CI, and RMSE
#   para.list: a list of parameter vector
#   G: number of repeats
#   true.para: a vector for true parameter
rmse.simu <- function(para.list, true.para, G, burn.in = 0) {
   if (burn.in >= G)
      stop("burn in is too large")
   P <- length(true.para)
   # find point estimate
   para.est <- sapply(para.list, FUN = para.stat, G = G, burn.in = burn.in)
   para.est <- para.est[1:P, ]
   # absolute and relative bias
   S <- length(para.list)
   para.diff <- para.est - true.para
   abs.bias <- apply(para.diff, FUN = mean, MARGIN = 1)
   rel.bias <- abs.bias / true.para * 100
   # rmse
   rmse <- apply(para.diff * para.diff, FUN = sum, MARGIN = 1)
   rmse <- rmse / (S - 1)
   # coverage proability
   ci.cov <- sapply(para.list, true.para, FUN = para.CI.cov, G = G, burn.in = burn.in)
   ci.cov <- apply(ci.cov, FUN = mean, MARGIN = 1)
   # for average
   res <- list(abs.bias = abs.bias, rel.bias = rel.bias, rmse = rmse, ci.cov = ci.cov)
   # for each data
   abs.bias <- t(para.diff)
   rel.bias <- para.diff / true.para * 100
   rel.bias <- t(rel.bias)
   rmse <- t(para.diff * para.diff)
   ci.cov <- sapply(para.list, true.para, FUN = para.CI.cov, G = G, burn.in = burn.in)
   ci.cov <- t(ci.cov)
   res <- list(abs.bias = abs.bias, rel.bias = rel.bias, rmse = rmse, ci.cov = ci.cov)
   return(res)
}

# box plot  
my.box <- function(in.dat, sum.stat = c("abs.bias", "rel.bias", "rmse"), para.index = 1, 
   para.type = "beta", para.sub = para.index, main.title = "") {
   if (length(in.dat) != 09 && length(in.dat) != 10 && length(in.dat) != 12){
      stop("in my.box, length of in.dat must be 9, 10, or 12\n")
   }
   num.g <- length(in.dat) # number of data sets
   names.g <- names(in.dat) # names for data set
   n <- rep(NA, num.g) # number of rows for each data set
   for (i in 1:num.g) {
      n[i] <- nrow(in.dat[[i]][[sum.stat]])
   }
   if (any(is.na(n))) {
      stop("in my.box, some data sets are incomplete \n");
   }
   plot.g <- rep(names.g, n) # data for groups
   plot.g <- factor(plot.g, levels = names.g) 
   plot.y <- rep(NA, sum(n)) # data for stats
   upper.index <- cumsum(n)
   lower.index <- upper.index - n + 1
   for (i in 1:num.g) {
      plot.y[lower.index[i]:upper.index[i]] <- in.dat[[i]][[sum.stat]][ , para.index] 
   }
   
   # plot boxplot
   # for labels of y-axis
   if (sum.stat == "abs.bias")
      y.txt <- "Absolute Bias"
   if (sum.stat == "rel.bias")
      y.txt <- "Relative Bias"
   if (sum.stat == "rmse")
      y.txt <- "Root Mean Square Error"
   
   if (length(in.dat) == 10) {
      # for color
      plot.col <- c(2, rep(3:5, 3))
      # for positions
      plot.pos <- c(1, 3:5, 7:9, 11:13) 
      # for tick names of x-axis
      plot.names <- names.g
      plot.names <- c("Complete", "", "Xiao", "", "", "Mice", "", "", "Norm", "")
      plot.names <- rep("", 10)
      t.pos <- c(1, 3, 6, 9) + 0.5
      t.names <- c("Complete", "Xiao", "Mice", "Norm")
   }
   if (length(in.dat) == 12) {
      # for color
      plot.col <- rep(2:5, 3)
      # for positions
      plot.pos <- c(1:4, 6:9, 11:14) 
      # for tick names of x-axis
      plot.names <- c("", "PX-GSM", "", "", "", "PX-GS", "", "", "", "PX-MH", "", "")
      plot.names <- rep("", 12)
      t.pos <- c(2, 7, 12) + 0.5
      t.names <- c("PX-GSM", "PX-GS", "PX-MH")
   }
   if (length(in.dat) == 9) {
      # for color
      plot.col <- rep(3:5, 3)
      # for positions
      plot.pos <- c(1:3, 5:7, 9:11) 
      # for tick names of x-axis
      plot.names <- c("", "PX-GSM", "", "", "PX-GS", "", "", "PX-MH", "")
      plot.names <- rep("", 9)
      t.pos <- c(2, 6, 10) + 0.0
      t.names <- c("PX-GSM", "PX-GS", "PX-MH")
   }   
   
   k <- eval(para.sub)
   nt <- eval(main.title)
   if (para.type == "r")
      new.main <- substitute(expression(paste(main.title, " ", r[k])))
   if (para.type == "beta")
      new.main <- substitute(expression(paste(main.title, " ", beta[k])))
   if (para.type == "gamma")
      new.main <- substitute(expression(paste(main.title, " ", gamma[k])))

   plot(plot.y ~ plot.g, ylab = y.txt, xlab = "", main = new.main, 
      col = plot.col, at = plot.pos, names  = plot.names, xaxt = "n")
   mtext(side = 1, line = 0.5, at = t.pos, text = t.names, cex = 0.75)
}

# box plot
my.box.3 <- function(in.dat, sum.stat = c("abs.bias", "rel.bias", "rmse"), para.index = 1, 
   main.title = "Boxplots for Coefficients") {
   if (length(in.dat) != 12) {
      stop("in my.box.3, length of in.dat is not 12\n")
   }   
   num.g <- length(in.dat) # number of data sets
   names.g <- names(in.dat) # names for data sets
   n <- rep(NA, num.g) # number of rows for each data set
   for (i in 1:num.g) {
      n[i] <- nrow(in.dat[[i]][[sum.stat]])
   }
   if (any(is.na(n))) {
      stop("in my.box, some data sets are incomplete \n")
   }
   plot.g <- rep(names.g, n) # data for groups
   plot.g <- factor(plot.g, levels = names.g) 
   plot.y <- rep(NA, sum(n)) # data for stats
   upper.index <- cumsum(n)
   lower.index <- upper.index - n + 1
   for (i in 1:num.g) {
      plot.y[lower.index[i]:upper.index[i]] <- in.dat[[i]][[sum.stat]][ , para.index] 
   }
   
   # plot boxplot
   # for labels of y-axis
   if (sum.stat == "abs.bias")
      y.txt <- "Absolute Bias"
   if (sum.stat == "rel.bias")
      y.txt <- "Relative Bias"
   if (sum.stat == "rmse")
      y.txt <- "Root Mean Square Error"
   # for color
   plot.col <- rep(2:5, 3)
   # for positions
   plot.pos <- c(1:4, 6:9, 11:14) 
   # for tick names of x-axis
   plot.names <- names.g
   plot.names <- c("", "PX-1", "", "", "", "PX-2", "", "", "", "PX-3", "", "")
   plot(plot.y ~ plot.g, ylab = y.txt, xlab = "Methods", main = main.title,
      at = plot.pos, col = plot.col, names = plot.names)
}

# function for boxplot: only two methods
# number of parameters can be 2, 3, etc.
# for each parameter, boxplot for two methods at the same 
my.box4 <- function(gs, gsm, para.type = c("beta", "gamma", "r", "sigma")) {
   # number of simulations/data sets
   nd <- nrow(gs) 
   
   # create data matrix for plot
   plot.y <- rbind(gs, gsm)
   plot.y <- matrix(plot.y, nrow = nd) # data matrix 

   # number of parameters
   nk <- ncol(plot.y) / 2 
   
   # positions for boxplot
   plot.pos <- c(1:(2 * nk))
   plot.pos <- plot.pos + rep(0:(nk - 1), rep(2, nk))
   
   # colors for boxplot
   plot.col <- rep(c("grey80", "grey90"), nk)
   
   # plot without x-axis labels
   boxplot(plot.y, at = plot.pos, col = plot.col, xaxt = "n",
      main = "", sub = "Parameters")
   abline(h=0, col="black")
	    
   # add legend
   # legend(x = 1, y = 0.3, fill = c("grey80", "grey90"), legend = c("gs", "gsm"))
   
   # add x-axis labels
   plot.pos <- matrix(plot.pos, nrow = 2)
   plot.pos <- apply(plot.pos, MARGIN = 2, FUN = mean)
   for (k in 1:nk) {
      if (para.type == "beta") {
         plot.names <- bquote(beta[.(k)])
      } else if (para.type == "gamma") {
         plot.names <- bquote(gamma[.(k)])      
      } else if (para.type == "r") {
         plot.names <- bquote(r[.(k)])
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
   #plot.col <- rep(c("grey70", "grey80", "grey90"), nk)
   plot.col <- rep(c("orchid", "slateblue", "gray48"), nk)
   
   # plot without x-axis labels
   boxplot(plot.y, at = plot.pos, col = plot.col, xaxt = "n",
      main = "", sub = "Parameters")
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
         plot.names <- bquote(gamma[.(k)])      
      } else if (para.type == "r") {
         plot.names <- bquote(r[.(k)])
      }
	  else {
		plot.names <- bquote(sigma[.(k)])
	  }
      axis(side = 1, at = plot.pos[k], labels = plot.names, tick = TRUE)
   }
}


simu.Ord <- function(n, k, theta, rho, C = NULL)
{
	# simulate data for X
	#p <- length(theta)
	#if(p != 2) stop("length of coefficient must be at 2 at this stage!/n")
	X <- matrix(NA, n * k, 2)
	#X[, 1] <- rep(1, n * k) # for intercept
	#a = runif(k,min=-0.1,max=0.8)
      #X[, 2] = rep(a, n)
	X[,1] = runif(n*k,min=-0.5,max=0.5)
      X[,2] = runif(n*k,min=-0.5,max=0.5)

	# simulate data for Y and Z
	# idealy the covariance matrix should be supplied by users
	corr <- diag(1, k)	
	for(i in 1:k)
		for(j in 1:k)
		{
			corr[i,j] <- rho ^(abs(i-j))					
		}
	#mu0 <- as.vector(X[1:k,] %*% theta)

	Z <- matrix(0, n, k)
	#Z <- matrix(0, n, k)
	for(i in 1:n) 
	{
		t1 = (i-1)*5+1
		t2 = (i-1)*5+5
		mu = X[t1:t2,]%*% theta
		Z[i,] <- rmvnorm(1, mu, corr)
		#Z[i,] <- D%*%Z[i,]
	}
	Y <- matrix(1, n, k)
	for(i in 1:n)
	{
		for (j in 1:k)
		{
			if (Z[i,j]< 0) Y[i,j]=0
			else if (Z[i,j]<1) Y[i,j]=1
			else if (Z[i,j]<2) Y[i,j]=2
			else Y[i,j]=3
		}	
	}		
	# for output
	Z <- as.data.frame(Z)
	Y <- as.data.frame(Y)
	names(Z) <- paste("r", 1:k, sep = "")
	names(Y) <- paste("r", 1:k, sep = "")
	
	Z2 <- as.vector(t(Z))
	Y2 <- as.vector(t(Y))
	return(list(Y = Y, X = X, Z = Z, Y2 = Y2, Z = Z2))
}



simu.Ord1 <- function(n, k, theta, rho, C = NULL)
{
	# simulate data for X
	X <- matrix(NA, n * k, 3)
		
	X[,1] = rep(1,n*k)
      X[,2] = runif(n*k,min=-0.5,max=0.5)
      X[,3] = runif(n*k,min=-0.5,max=0.5)

	# simulate data for Y and Z
	# idealy the covariance matrix should be supplied by users
	corr <- diag(1, k)	
	for(i in 1:k)
		for(j in 1:k)
		{
			corr[i,j] <- rho ^(abs(i-j))					
		}
	#mu0 <- as.vector(X[1:k,] %*% theta)

	Z <- matrix(0, n, k)
	#Z <- matrix(0, n, k)
	for(i in 1:n) 
	{
		t1 = (i-1)*5+1
		t2 = (i-1)*5+5
		mu = X[t1:t2,]%*% theta
		Z[i,] <- rmvnorm(1, mu, corr)
		#Z[i,] <- D%*%Z[i,]
	}
	Y <- matrix(1, n, k)
	for(i in 1:n)
	{
		for (j in 1:k)
		{
			if (Z[i,j]< 0) Y[i,j]=0
			else if (Z[i,j]<1) Y[i,j]=1
			else if (Z[i,j]<2) Y[i,j]=2
			else Y[i,j]=3
		}	
	}		
	# for output
	Z <- as.data.frame(Z)
	Y <- as.data.frame(Y)
	names(Z) <- paste("r", 1:k, sep = "")
	names(Y) <- paste("r", 1:k, sep = "")
	
	Z2 <- as.vector(t(Z))
	Y2 <- as.vector(t(Y))
	return(list(Y = Y, X = X, Z = Z, Y2 = Y2, Z = Z2))
}


simu.Ord2 <- function(n, k, theta, rho, C = NULL)
{
	# simulate data for X
		
	X1 = rep(c(1,0,0,0,0),n)
      X2 = rep(c(0,1,0,0,0),n)
	X3 = rep(c(0,0,1,0,0),n)
	X4 = rep(c(0,0,0,1,0),n)
	X5 = rep(c(0,0,0,0,1),n)

	X = cbind(X1, X2, X3, X4, X5)


	# simulate data for Y and Z
	# idealy the covariance matrix should be supplied by users
	corr <- diag(1, k)	
	for(i in 1:k)
		for(j in 1:k)
		{
			corr[i,j] <- rho ^(abs(i-j))					
		}
	#mu0 <- as.vector(X[1:k,] %*% theta)

	Z <- matrix(0, n, k)
	#Z <- matrix(0, n, k)
	for(i in 1:n) 
	{
		t1 = (i-1)*5+1
		t2 = (i-1)*5+5
		mu = X[t1:t2,]%*% theta
		Z[i,] <- rmvnorm(1, mu, corr)
		#Z[i,] <- D%*%Z[i,]
	}
	Y <- matrix(1, n, k)
	for(i in 1:n)
	{
		for (j in 1:k)
		{
			if (Z[i,j]< 0) Y[i,j]=0
			else if (Z[i,j]<1) Y[i,j]=1
			else if (Z[i,j]<2) Y[i,j]=2
			else Y[i,j]=3
		}	
	}		
	# for output
	Z <- as.data.frame(Z)
	Y <- as.data.frame(Y)
	names(Z) <- paste("r", 1:k, sep = "")
	names(Y) <- paste("r", 1:k, sep = "")
	
	Z2 <- as.vector(t(Z))
	Y2 <- as.vector(t(Y))
	return(list(Y = Y, X = X, Z = Z, Y2 = Y2, Z = Z2))
}


simu.Ord3 <- function(n = 100, k = 5, theta = c(0.5, 0.5, 1, 1, 1.5), rho = 0.5, C = NULL)
{
	# simulate data for X
		
	X1 = rep(c(1,0,0,0,0),n)
      X2 = rep(c(0,1,0,0,0),n)
	X3 = rep(c(0,0,1,0,0),n)
	X4 = rep(c(0,0,0,1,0),n)
	X5 = rep(c(0,0,0,0,1),n)

	X = cbind(X1, X2, X3, X4, X5)


	# simulate data for Y and Z
	# idealy the covariance matrix should be supplied by users
	corr <- diag(1, k)	
	for(i in 1:k)
		for(j in 1:k)
		{
			corr[i,j] <- rho ^(abs(i-j))					
		}
	#mu0 <- as.vector(X[1:k,] %*% theta)

	Z <- matrix(0, n, k)
	#Z <- matrix(0, n, k)
	for(i in 1:n) 
	{
		t1 = (i-1)*5+1
		t2 = (i-1)*5+5
		mu = X[t1:t2,]%*% theta
		Z[i,] <- rmvnorm(1, mu, corr)
		#Z[i,] <- D%*%Z[i,]
	}
	Y <- matrix(1, n, k)
	for(i in 1:n)
	{
		if (Z[i,1]< 0) Y[i,1]=0
		else if (Z[i,1]<0.5) Y[i,1]=1
		else if (Z[i,1]<1) Y[i,1]=2
		else Y[i,1]=3
	}	
	#table(Y[,1])

	for(i in 1:n)
	{
		if (Z[i,2]< 0) Y[i,2]=0
		else if (Z[i,2]<0.5) Y[i,2]=1
		else if (Z[i,2]<1) Y[i,2]=2
		else Y[i,2]=3
	}
	#table(Y[,2])
	
	for(i in 1:n)
	{
		if (Z[i,3]< 0) Y[i,3]=0
		else if (Z[i,3]<1) Y[i,3]=1
		else if (Z[i,3]<2) Y[i,3]=2
		else Y[i,3]=3
	}
	#table(Y[,3])

	for(i in 1:n)
	{
		if (Z[i,4]< 0) Y[i,4]=0
		else if (Z[i,4]<1) Y[i,4]=1
		else if (Z[i,4]<2) Y[i,4]=2
		else Y[i,4]=3
	}
	#table(Y[,4])

	for(i in 1:n)
	{
		if (Z[i,5]< 0) Y[i,5]=0
		else if (Z[i,5]<1.5) Y[i,5]=1
		else if (Z[i,5]<2.5) Y[i,5]=2
		else Y[i,5]=3
	}
	#table(Y[,5])
	

	
	# for output
	Z <- as.data.frame(Z)
	Y <- as.data.frame(Y)
	names(Z) <- paste("r", 1:k, sep = "")
	names(Y) <- paste("r", 1:k, sep = "")
	
	Z2 <- as.vector(t(Z))
	Y2 <- as.vector(t(Y))
	return(list(Y = Y, X = X, Z = Z, Y2 = Y2, Z = Z2))
}


# obtain average of each sampling step for parameters
#   para.list: lists of sampled parameter vector
para.list.step <- function(para.list) {
   S <- length(para.list)
   L <- length(para.list[[1]])
   res.m <- matrix(NA, nrow = S, ncol = L)
   for(s in 1:S) {
      res.m[s, ] <- para.list[[s]]
   }
   
   res <- apply(res.m, FUN = mean, MARGIN = 2, na.rm = TRUE)
   
   return(res)   
}
