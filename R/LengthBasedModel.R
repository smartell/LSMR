# --------------------------------------------------------------------------- #
#                             LENGTH BASED MODEL                              #
# Author     : Steven Martell                                                 #
# Date       : July18, 2012                                                   #
# Description: R-code for a length-based model initialized at equilibrium.    #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# SIMULATION MODEL PARAMETERS                                                 #
# --------------------------------------------------------------------------- #
rbar     <- 100.		# Steady-state average recruitment.
minf     <- 0.30		# Asymptotic instantaneous natural mortality rate.
linf     <- 40.0		# Asymptotic length
vbk      <- 0.20		# von Bertalanffy growth coefficient.
beta     <- 0.95		# Variance parameter in length-distribution.
mu.r     <- 12.0		# Mean length of new recruits.
cv.r     <- 0.15		# CV of mean length of new recruits.
xmid     <- seq(5.5, 45.5, by=1)	# mid point of length intervals

theta    <- list(rbar=rbar, minf=minf, linf=linf, vbk=vbk, beta=beta
                 , mu.r=mu.r, cv.r=cv.r)

# --------------------------------------------------------------------------- #
# SIZE TRANSITION MATRIX                                                      #
# args: x is the mid point of the length intervals                            #
# args: pars is a list object with linf, vbk, and beta params                 #
# --------------------------------------------------------------------------- #
.calcSizeTransitionMatrix <- function(x, pars)
{
	with(as.list(pars), {
		lambda  <- exp(log((linf-x)*(1.0-exp(-vbk))+1))
		alpha   <- lambda/beta
		fn <- function(i)
		{
			z2 <- pgamma(xmid-xmid[i]+0.5,alpha[i])
			z1 <- pgamma(xmid-xmid[i]-0.5,alpha[i])
			return(z2-z1)
		}
		P = t(sapply(1:n,fn))
		P = P/rowSums(P)
		return(P)
	})
}

# --------------------------------------------------------------------------- #
# LENGTH-BASED MORTALITY                                                      #
# --------------------------------------------------------------------------- #
.calcLengthBasedMortality <- function(x, pars)
{
	with(as.list(pars), {
		mx   <- minf*linf/x
		return(mx)
	})
}

# --------------------------------------------------------------------------- #
# LENGTH DISTRIBUTION OF NEW RECRUITS                                         #
# --------------------------------------------------------------------------- #
.calcRecruitsAtLength <- function(x, pars)
{
	with(as.list(pars), {
		a    <- 1/(cv.r)^2
		b    <- mu.r/a
		rx   <- dgamma(x, shape=a, scale=b)
		rx   <- rx/sum(rx)
		
		return(rx)
	})
}


# --------------------------------------------------------------------------- #
# LENGTH-BASED POPULATION MODEL                                               #
# --------------------------------------------------------------------------- #
.calcNumbersAtLength <- function(x, pars)
{
	with(as.list(pars), {
		# indexes
		n     <- length(x)
		nyr   <- 100
		
		mx    <- .calcLengthBasedMortality(x, pars)
		P     <- .calcSizeTransitionMatrix(x, pars)
		rx    <- .calcRecruitsAtLength(x, pars)
		
		# Calculate per-recruit survivorship matrix
		phi      <- matrix(0, nrow=n, ncol=n)
		phi[1, ] <- rx
		for(i in 2:n)
		{
			phi[i, ] <- as.vector((phi[i-1, ]*exp(-mx)) %*% P)
			if( i==n )
			{
				phi[i, ] <- phi[i, ]/( 1.0 - exp(-mx) )
			}
		}
		
		# Calculate initial numbers-at-length
		ri     <- rep(rbar, length=n)
		rt     <- rep(rbar, length=nyr)
		N      <- matrix(0, nrow=nyr, ncol=n)
		N[1, ] <- as.vector(ri %*% phi)
		
		# Calculate number-at-length over time
		for(i in 2:nyr)
		{
			N[i, ]  <- as.vector((N[i-1, ]*exp(-mx)) %*% P) + rt[i]*rx
		}
		
		# Return N matrix
		return(N)
	})
}






