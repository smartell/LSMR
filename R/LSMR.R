# --------------------------------------------------------------------------- #
# R-script for plotting output from LSMR.rep
# Author: Steven Martell
# DATE: July 9, 2012
# LOCATION: Wailea Hawaii
# DEPENDENCIES: read.admb.R
#





# --------------------------------------------------------------------------- #
source("read.admb.R", echo=FALSE)

.SAVEFIGS		<- TRUE
.PRESENTATION	<- TRUE
.CEXLAB			<- 1.5
.FILENAME		<- "../ADMB/srcLSMR/lsmr"
.FIGDIR			<- "../FIGS/LSMR/SIMb/"



getObj <- function(fn=.FILENAME)
{
	obj				<- read.admb(.FILENAME)
	class(obj)		<- c(class(obj), "lsmr")
	return(obj)
}

obj <- getObj()

main <- function()
{
	obj <- getObj()
	plot(obj);
}

# --------------------------------------------------------------------------- #
# S3 Method for class 'lsmr'                                                  #
# --------------------------------------------------------------------------- #
print.lsmr <- function(obj, ...)
{
	print(attributes(obj))
}

plot.lsmr <- function(obj, ...)
{
	opar <- par(no.readonly=TRUE)
	
	
	with(obj, {		
		par(mfcol=c(2, 2), las=1, mar=c(5.1, 4.1, 4.1, 2.1), oma=c(1, 1, 1, 1))
		
		if(.SAVEFIGS) 
		{
			gfn <- paste(.FIGDIR, "fig:LSRM%d.png", sep="")
			png(gfn)
		}
		par(cex.lab = .CEXLAB)
		
		jmin  <- min(which(xmid>=15))
		nx    <- dim(N)[2]
		
		.plotNt(obj)
		gletter(1)
		
		.plotRt(obj)
		gletter(2)
		
		.plotFt(obj)
		gletter(3)
		
		.plotMx(obj)
		gletter(4)
		
		
		par(mfcol=c(1, 1))
		.plotSelex(xmid, t(sx))
		
		
		
		ir = 0
		for(i in 1:ngear)
		{
			ir = 1:irow[i] + max(ir)
			
			.plotLF(xmid, i_C[ir, ], Chat[ir, ], "Total Catch")
		
			.plotLF(xmid, i_M[ir, ], Mhat[ir, ], "New Marks")
		
			.plotLF(xmid, i_R[ir, ], Rhat[ir, ], "Recaptures")
		}
		
		dev.off()
	})
	par(opar)
}

.plotMx <- function(obj, ...)
{
	with(obj, { 
		yl <- c(0, 1.2*max(mx))
		if(exists("true_mx"))yl <- c(0, 1.2*max(mx, true_mx))
		
		plot(xmid, mx, type="l", ylim=yl
		, xlab="Size class (cm)", ylab="Natural mortality")
		abline(h=m_infty, col=2)
		points(l_infty, m_infty, pch=19)
		text(l_infty, m_infty, "M at Linf", pos=3)
		
		if(exists("true_mx"))
		{
			lines(xmid, true_mx, lwd=2, col=colr(2, 0.5))
			legend("top", c("estimated","true"),lty=1, lwd=c(1, 2)
			, col=c(1, colr(2, 0.5)), bty="n")
		}
		
	})
}

.plotFt <- function(obj, ...)
{
	with(obj, {
		f_yr = seq(min(yr), max(yr), length=dim(fi)[2])
		matplot(f_yr, t(fi), type="l", ylim=c(0, max(fi))
		, xlab="Year", ylab="Capture probability", col=1, lty=1:2)
		
		if(exists("true_fi"))
		{
			matlines(f_yr, t(true_fi), type="l",lwd=2, col=colr(2, 0.5))	
		}
		legend("top", c("Tramel", "Hoop"), lty=1:2, bty="n")
	})
}

.plotRt <- function(obj, scale=1000, ...)
{
	with(obj,{
		yl <- c(0, max(Rt)/scale)
		if(exists("true_Rt"))yl <- c(0, max(Rt, true_Rt)/scale)
		
		plot(yr, Rt/scale, type="h", ylim = yl 
		, xlab="Year", ylab="Annual recruits (1000s)")
		if(exists("true_Rt"))points(yr, true_Rt/scale, pch=20, col=colr(2, 0.5))
	})	
}

.plotNt <- function(obj, ...)
{
	# plot Numbers greater than 50 mm 
	with(obj, {
		yl <- c(0, max(Nt))
		if(exists("true_Nt"))
			yl <- c(0, max(Nt, true_Nt))
		
		plot(yr, Nt, type="l", ylim=yl
		, xlab="Year", ylab="Abundance (> 50 mm)")
	
		if(exists("true_Nt"))
		{
			lines(yr, true_Nt, lwd=2, col=colr(2, 0.5))
			legend("top", c("estimated","true"),lty=1, lwd=c(1, 2)
			, col=c(1, colr(2, 0.5)), bty="n")
		}
		 
	})	
}

.plotSelex <- function(x, y, ...)
{
	#Plot seletivities
	matplot(x, y, xlab="Length (cm)", ylab="Selectivity", type="l", col=1)
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


.plotLF <- function(x, O, P, main="", ...)
{
	# This funciton plots the observed (O) and predicted (P) 
	# length frequency distributions.
	opar <- par(no.readonly=TRUE)
	ir   <- which(rowSums(P)!=0,arr.ind=TRUE)
	
	O    <- O[ir, ]
	P    <- P[ir, ]
	n    <- dim(O)[1]
	nr   <- ceiling(sqrt(n))
	nc   <- ceiling(n/nr)
	
	par(mfcol=c(nr, nc), mar=c(0, 0, 0, 0))
	par(oma=c(5, 5, 5, 4), font.main=1)
	ymax = max(O[, -1:-2], P)
	for(i in 1:n)
	{
		plot(x, O[i, -1], type="n", xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0, ymax))
		.staircase(x, O[i, -1], border=NA)
		lines(x, P[i, ], col=1)
		abline(v=c(10, 15), lty=3)
		title(main=O[i, 1], line=-1)
		mfg <- par(no.readonly=TRUE)$mfg
		if(mfg[2]==1 && mfg[1]%%2) axis(2)
		if(mfg[2]==nc && !mfg[1]%%2) axis(4)
		if(mfg[1]==nr && mfg[2]%%2) axis(1)
		if(mfg[1]==1 && !mfg[2]%%2) axis(3)
	}
	mtext(c("Size class (cm)", "Frequency"), side=c(1, 2), outer=TRUE, line=2.5, las=0)
	mtext(main, side=3, outer=TRUE, line=2.5)
	par(opar)
}

