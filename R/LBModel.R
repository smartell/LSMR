## Testing some Ideas for a length Transistion Matrix for
## a Humpback Chub assessment based on statistical catch-at-length.

## Note that the pgamma function is closely related to the 
## incomplete gamma function (which is used in ADMB cumd_gamma() ).
## the incomplete gamma function is available in the package
## zipfR as Igamma(a, x)

## pgamma(x, a, b=1) is equivalent to Igamma(a, x)/gamma(a)

require(zipfR)

## Growth parameters
linf	<- 40.0
k		<- 0.25
beta	<- 0.75
n		<- 35	#number of intervals


## Bin intervals (mm)
x	<- seq(5.0, 1.0*linf, length=n+1)
d	<- 0.5*(x[2]-x[1])
xp	<- seq(x[1]+d, by=2*d, length=n)

dl	<- (linf-xp)*(1-exp(-k))
alpha <- dl/beta*6/12
P	<- matrix(0, nrow=n, ncol=n)
for(i in 1:n)
{
	dx = x-x[i]+d
	print(dx)
	z <- Igamma(alpha[i], dx)#-Igamma(alpha[i], dx-d)
	print(z)
	P[i-1,] <- diff(z)/sum(diff(z)) 

}
matplot(xp, t(P), type="l")

f = rep(0, length=n)
f[7] = 10
#par(mfcol=c(1, 3))
plot(xp, f, type="S", ylim=c(0, 12), lwd=2, xlab="Size interval (cm)", ylab="Frequency")
for(i in 2:4)
{
	f=f%*%P
	lines(xp+(i-2)*3/12, f, type="S", col=i, lwd=2)
}
legend("top",c("Age 0", "Age 0.5", "Age 1.0", "Age 1.5"), lty=1, col=1:4, ncol=4, lwd=2, bty="n")
dev.copy2pdf(file="../FIGS/fig:lengthTransition.pdf")