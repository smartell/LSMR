# A Length Transition Model based on the von Bertalanffy growth curve.
# Author: Steven Martell
# When: Aboard AC flight 033 YVR to SYD
require(hacks)
# Model Dimensions
tau <- 1							# time step
n   <- 50							# Number of length bins
xl  <- 5							# Smallest bin
bw  <- 1							# width of bin
xi  <- seq(xl, by=bw, length=(n+1))	# bin intervals
xm  <- xi[1:n]+0.5*bw				# bin midpoints

# Growth parameters
linf  <- 45							# asymptotic length
k     <- 0.15 * tau					# Brody growth coefficient
beta  <- 0.75						# Variance scaler in growth increment

# Growth increment & parameters for the Incomplete gamma function
ginc  <- log(exp((linf-xm)*(1-exp(-k)))+1)
alfa  <- ginc/beta

P     <- matrix(0, nrow=n, ncol=n)
for(i in 1:n)
{
	dx <- xi[i:(n+1)] - xi[i]
	z  <- Igamma(alfa[i], dx)
	P[i, i:n] <- diff(z)/sum(diff(z))
}


# Graphic for the length transition matrix
matplot(P, xlim=range(xm), ylim=range(xm), type="n", xlab="Length (cm)", ylab="", yaxt="n")
for(i in 1:n)
{
	p = P[i,]/(0.8*max(P[i,]))
	#lines(xm,((n-i)+min(xi))+p)
	xx = c(xm, rev(xm))
	yy = c((n-i)+min(xi) + p,(n-i)+min(xi)+0*p )
	polygon(xx, yy, border=NA, col=colr(4, 0.5))
}