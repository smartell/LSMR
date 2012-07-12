# --------------------------------------------------------------------------- #
# R-script for plotting output from LSMR.rep
# Author: Steven Martell
# DATE: July 9, 2012
# LOCATION: Wailea Hawaii
# DEPENDENCIES: read.admb.R
#





# --------------------------------------------------------------------------- #
source("read.admb.R", echo=FALSE)

.SAVEFIGS		<- FALSE
.PRESENTATION	<- FALSE
.FILENAME		<- "../ADMB/srcLSMR/lsmr"
.FIGDIR			<- "../FIGS/LSMR/"

obj				<- read.admb(.FILENAME)
class(obj)		<- c(class(obj), "lsmr")





# --------------------------------------------------------------------------- #
# S3 Method for class 'lsmr'                                                  #
# --------------------------------------------------------------------------- #
print.lsmr <- function(obj, ...)
{
	print(attributes(obj))
}

plot.lsmr <- function(obj, ..1000.)
{
	opar <- par(no.readonly=TRUE)
	
	with(obj, {
		par(mfcol=c(2, 2), las=1)
		plot(yr, Nt, type="l", ylim=c(0, max(Nt))
		, xlab="Year", ylab="Abundance (numbers > 30 mm)")
		if(exists("true_Nt")) lines(yr, true_Nt, col=2, lwd=2)
		gletter(1)
		
		plot(yr, Rt, type="h", ylim=c(0, max(Rt)) 
		, xlab="Year", ylab="Annual recruits ")
		if(exists("true_Rt"))points(yr, true_Rt, pch=20, col=2)
		gletter(2)
		
		
		plot(yr, fi, type="o", ylim=c(0, max(fi))
		, xlab="Year", ylab="Capture probability")
		if(exists("true_fi"))lines(yr, true_fi, type="o", pch=2, col=2)
		gletter(3)
		
		plot(xmid, mx, type="o", ylim=c(0, max(mx))
		, xlab="Size class (cm)", ylab="Natural mortality")
		abline(h=m_infty, col=2)
		gletter(4)
		par(opar)
		
		.plotLF(xmid, i_C, Chat)
		
		.plotLF(xmid, i_M, Mhat)
		
		.plotLF(xmid, i_R, Rhat)
	})
	par(opar)
}

.plotLF <- function(x, O, P, ...)
{
	# This funciton plots the observed (O) and predicted (P) 
	# length frequency distributions.
	opar <- par(no.readonly=TRUE)
	n  <- dim(O)[1]
	nr <- round(sqrt(n))
	nc <- round(n/nr)
	par(mfcol=c(nr, nc), mar=c(0, 0, 0, 0))
	par(oma=c(5, 5, 3, 1), font.main=1)
	ymax = max(O[, -1:-2], P)
	for(i in 1:n)
	{
		plot(x, O[i, -1:-2], type="h", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0, ymax))
		lines(x, P[i, ], col=4)
		title(main=O[i, 1], line=-1)
		mfg <- par(no.readonly=T)$mfg
		if(mfg[2]==1) axis(2)
		if(mfg[1]==nr) axis(1)
	}
	mtext(c("Size class (cm)", "Frequency"), side=c(1, 2), outer=TRUE, line=2.5, las=0)
	par(opar)
}

