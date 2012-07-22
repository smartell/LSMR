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
		par(mfcol=c(2, 2), las=1, mar=c(5.1, 4.1, 4.1, 2.1))
		plot(yr, Nt, type="l", ylim=c(0, max(Nt))
		, xlab="Year", ylab="Abundance (numbers > 30 mm)")
		if(exists("true_Nt")) lines(yr, true_Nt, col=2, lwd=2)
		gletter(1)
		
		plot(yr, Rt, type="h", ylim=c(0, max(Rt)) 
		, xlab="Year", ylab="Annual recruits ")
		if(exists("true_Rt"))points(yr, true_Rt, pch=20, col=2)
		gletter(2)
		
		f_yr = seq(min(yr), max(yr), length=dim(fi)[2])
		matplot(f_yr, t(fi), type="l", ylim=c(0, max(fi))
		, xlab="Year", ylab="Capture probability", col=1)
		if(exists("true_fi"))
			matlines(f_yr, t(true_fi), type="l", col=2)
		gletter(3)
		
		plot(xmid, mx, type="o", ylim=c(0, max(mx))
		, xlab="Size class (cm)", ylab="Natural mortality")
		abline(h=m_infty, col=2)
		gletter(4)
		par(opar)
		
		
		ir = 0
		for(i in 1:ngear)
		{
			ir = 1:irow[i] + max(ir)
			.plotLF(xmid, i_C[ir, ], Chat[ir, ])
		
			.plotLF(xmid, i_M[ir, ], Mhat[ir, ])
		
			.plotLF(xmid, i_R[ir, ], Rhat[ir, ])
		}
	})
	par(opar)
}

.staircase <- function(x, y, ...)
{
	# Add a shaded staircase polygon to a plot.
	dx <- 0.5*(x[2]-x[1])
	xp <- as.vector(rbind(x-dx, x+dx))
	yp <- as.vector(rbind(y, y))
	xp <- c(min(xp), xp, max(xp))
	yp <- c(0, yp, 0)
	polygon(xp, yp, col=colr(1, 0.25), ...)
}


.plotLF <- function(x, O, P, ...)
{
	# This funciton plots the observed (O) and predicted (P) 
	# length frequency distributions.
	opar <- par(no.readonly=TRUE)
	ir   <- which(rowSums(P)>0,arr.ind=T)
	O    <- O[ir, ]
	P    <- P[ir, ]
	n    <- dim(O)[1]
	nr   <- round(sqrt(n))
	nc   <- round(n/nr)
	
	par(mfcol=c(nr, nc), mar=c(0, 0, 0, 0))
	par(oma=c(5, 5, 3, 4), font.main=1)
	ymax = max(O[, -1:-2], P)
	for(i in 1:n)
	{
		plot(x, O[i, -1:-2], type="n", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0, ymax))
		.staircase(x, O[i, -1:-2], border=NA)
		lines(x, P[i, ], col=1)
		title(main=O[i, 1], line=-1)
		mfg <- par(no.readonly=T)$mfg
		if(mfg[2]==1 && mfg[1]%%2) axis(2)
		if(mfg[2]==nc && !mfg[1]%%2) axis(4)
		if(mfg[1]==nr && mfg[2]%%2) axis(1)
		if(mfg[1]==1 && !mfg[2]%%2) axis(3)
	}
	mtext(c("Size class (cm)", "Frequency"), side=c(1, 2), outer=TRUE, line=2.5, las=0)
	par(opar)
}

