# --------------------------------------------------------------------------- #
# R-script for plotting output from ASMR.rep
# Author: Steven Martell
# DATE: March 27, 2012
# LOCATION: Seattle Washington
# DEPENDENCIES: read.admb.R






# --------------------------------------------------------------------------- #
source("read.admb.R", echo=FALSE)

.SAVEFIGS		<- TRUE
.PRESENTATION	<- TRUE
.FILENAME		<- "../ADMB/srcASMR/asmr"
.FIGDIR			<- "../FIGS/ASMR/"

obj				<- read.admb(.FILENAME)
obj$mc_nt4		<- read.table("../ADMB/srcASMR/nt4.mcmc", header=FALSE)
obj$mc_rt		<- read.table("../ADMB/srcASMR/rt.mcmc", header=FALSE)
obj$mc_pt		<- read.table("../ADMB/srcASMR/pt.mcmc", header=FALSE)
class(obj)		<- c(class(obj), "asmr")
class(obj$mta)		<- c("matrix", "marks")
class(obj$rta)		<- c("matrix", "marks")
class(obj$epsilon)	<- c("matrix", "marks")
class(obj$delta)	<- c("matrix", "marks")

obj$cmr <- list(mta=obj$mta, rcta=obj$rcta)
class(obj$cmr)		<- c("cmr", "list")

load("EstM_obj.RData")
load("FixM_obj.RData")

# --------------------------------------------------------------------------- #
# S3 method for class 'asmr'                                                  #
# --------------------------------------------------------------------------- #
plot.asmr <- function(obj, ask=dev.interactive(), ...)
{
	devAskNewPage(ask)
	with(obj,{
		plot(yr, nt4/1000, type="l", xlab="Year", ylab="Abundance (1000 age-4+)"
			,ylim=c(0, 1.2*max(nt4/1000)),  ...)
		grid(col=1)
		if(.SAVEFIGS)dev.copy2pdf(file=paste(.FIGDIR, "fig:nt4.pdf", sep=""))
		
			
		plot(byr, rt/1000, type="l", xlab="Brood year", ylab="Age-2 recruits", 
			,ylim=c(0, 1.2*max(rt/1000)), ...)
		grid(col=1)
		if(.SAVEFIGS)dev.copy2pdf(file=paste(.FIGDIR, "fig:rt.pdf", sep=""))
		
		plot(yr, pt, type="l", xlab="Year", ylab="Capture probability", ...)
		grid(col=1)
		if(.SAVEFIGS)dev.copy2pdf(file=paste(.FIGDIR, "fig:pt.pdf", sep=""))
	})
}

print.asmr <- function(obj, ...)
{
	print(attributes(obj))
}

plot.marks <- function(obj,x=NULL, y=NULL,  ...)
{
	# use PBS modelling to plot bubbles
	require(PBSmodelling)
	print(attributes(obj))
	
	z <- t(obj)
	plotBubbles(z, xval=x, yval=y, prettyaxis=TRUE, ...)
}

plot.cmr <- function(obj, x=NULL,  ...)
{
	# use PBS modelling to plot Bubbles.
	require(PBSmodelling)
	print(attributes(obj))
	
	# Construrct z-matrix with nyr rows and nyr+1 cols.
	# The first column is number of markes released by year.
	with(obj, {
		nyr 	<- dim(mta)[1]
		nage	<- dim(mta)[2]
		z   	<- matrix(0, ncol=nyr+1, nrow=nyr)
		z[, 1]	<- -rev((rowSums(mta)))
		
		# Recapture by cohort array
		R 		<- array(rcta,dim=c(nyr,nyr,nage))
		for (i in nyr:1)
		{
			z[i, 2:(nyr+1)] <- (rowSums(R[, nyr-i+1, ]))
		}
		
		# plotBubbles
		xv <- c(x[1], x)
		print(length(xv))
		plotBubbles(z,xval=xv,yval=paste(rev(x)),hide0=TRUE,prettyaxis=TRUE, ...)
	})

}

retrospective <- function(n=2)
{
	curwd <- getwd()
	setwd("../ADMB/srcASMR/")
	A<-list()
	i=0
	for(ryr in 0:n)
	{
		arg = paste("./ASMR"," -retro ", ryr, sep="")
		system(arg)
		
		i=i+1
		tmp	<- read.admb("asmr")
		A[[i]]	<- tmp
	}
	
	class(A) <- c("list", "retro")
	setwd(curwd)
	return(A)
}
plot.retro <- function(obj, ...)
{
	nretro=length(obj)

	# plot nt4 estimates
	yr <- obj[[1]]$yr
	nt4<- obj[[1]]$nt4
	plot(yr, nt4/1000, type="l", xlab="Year", ylab="Abundance (1000 age-4+)"
		,ylim=c(0, 1.2*max(nt4/1000)),  ...)
	
	for(i in 2:nretro)
	{
		yr <- obj[[i]]$yr
		nt4<- obj[[i]]$nt4
		ii <- length(yr)
		lines(yr, nt4/1000,col=i)
		points(yr[ii], nt4[ii]/1000, col=i, pch=20)
	}
	grid()
	if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:retro_nt4.pdf", sep=""))
	
	# plotr rt estimates
	yr <- obj[[1]]$byr
	rt<-  obj[[1]]$rt
	plot(yr, rt/1000, type="l", xlab="Brood year", ylab="Age-2 recruits"
		,ylim=c(0, 1.2*max(rt/1000)),  ...)
	
	for(i in 2:nretro)
	{
		yr <- obj[[i]]$byr
		rt <- obj[[i]]$rt
		ii <- length(yr)
		lines(yr, rt/1000,col=i)
		points(yr[ii], rt[ii]/1000, col=i, pch=20)
		#points(yr[ii], obj[[i]]$M, pch=19, col=i)
	}
	grid()
	if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:retro_rt.pdf", sep=""))
	
	# plotr rt estimates
	yr <- obj[[1]]$byr
	rt <- obj[[1]]$rt
	plot(yr, rt/1000, type="n", xlab="Brood year", ylab="Natural mortality"
		,ylim=c(0, 0.3,  ...))
	
	for(i in 2:nretro)
	{
		yr <- obj[[i]]$byr
		rt <- obj[[i]]$rt
		ii <- length(yr)
		#lines(yr, rt/1000,col=i)
		#points(yr[ii], rt[ii]/1000, col=i, pch=20)
		points(yr[ii], obj[[i]]$M, pch=19, col=i)
		
	}
	grid()
	
	
}

# --------------------------------------------------------------------------- #
# Saved Figures from ASMR                                                     #
# --------------------------------------------------------------------------- #
#
if(.PRESENTATION)
{
	require(hacks)
	par(las=1, bg=colr("white", 0.85), cex.lab=1.5, cex.axis=1.2,  
	    mar=c(5.1, 4.6, 4.1, 2.1))
}
# Number of marks released by age-year
plot(obj$mta,obj$yr,obj$age,hide0=TRUE,frange=0,xlab="Year",ylab="Age",size=0.25)
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:mta.pdf", sep=""))
#
# Number of recaptures by age-year
plot(obj$rta,obj$yr,obj$age,hide0=TRUE,frange=0,xlab="Year",ylab="Age",size=0.25,clrs=4)
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:rta.pdf", sep=""))
#
# Number of recaptures by cohort
#plot.cmr(obj$cmr,x=1989:2012,frange=0.1,xlab="Recapture year",ylab="Release year",size=0.25,clrs=4)
plot.cmr(obj$cmr,x=obj$yr, frange=0.02,xlab="Recapture year",size=0.4)
grid()
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:cmr.pdf", sep=""))
#
# Residuals in marks released by age-year
plot(obj$epsilon,obj$yr,obj$age,hide0=TRUE,frange=0,xlab="Year",ylab="Age")
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:epsilon.pdf", sep=""))
#
# Residuals in marks recaptured by age-year
plot(obj$delta,obj$yr,obj$age,hide0=TRUE,frange=0,xlab="Year",ylab="Age")
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:delta.pdf", sep=""))
#
# Plot retrospective estimates
if(!exists("RETRO"))RETRO	<- retrospective(8)
plot(RETRO)
#
# Uncertainty in age-4 abundance & recruits
boxplot(obj$mc_nt4/1000, range=0, names=obj$yr, ylim=c(0, max(obj$mc_nt4/1000)), xlab="Year", ylab="Abundance (1000)", lwd=0.25, col="blue")
boxplot(obj$mc_rt/1000, add=TRUE, range=0, names=NA, col="orange", lwd=0.25)
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:bxplt_Nt4_rt.pdf", sep=""))
#
# Uncertainty in age-4+ in 2011
hist(obj$mc_nt4[,23], prob=TRUE, yaxt="n", ylab="Posterior density", xlab="Age-4+ abundance in 2011", main="", col="tan")
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:marg_Nt4.pdf", sep=""))
#
# Sensitivity to M
# Uncertainty in age-4 abundance & recruits
boxplot(objM$mc_nt4/1000, range=0, names=obj$yr, ylim=c(0, max(objM$mc_nt4/1000)), xlab="Year", ylab="Abundance (1000)", lwd=0.25, col="lightblue")
boxplot(obj$mc_nt4/1000, add=TRUE, range=0, names=NA, col="blue", lwd=0.25)
boxplot(obj$mc_rt/1000, add=TRUE, range=0, names=NA, col="orange", lwd=0.25)
boxplot(objM$mc_rt/1000, add=TRUE, range=0, names=NA, col="yellow", lwd=0.25)
legend("top", c("Nt4 M=0.130","Rt  M=0.130","Nt4 M=0.094", "Rt  M=0.094")
, bty="n", pch=15, col=c("blue", "orange", "lightblue", "yellow"))
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:bxplt_Nt4_rt_estM.pdf", sep=""))
#
hist(obj$mc_nt4[,23], prob=TRUE, yaxt="n", ylab="Posterior density", xlab="Age-4+ abundance in 2011", main="", col="tan", xlim=c(8500, 13500))

hist(objM$mc_nt4[,23], prob=TRUE, yaxt="n", ylab="Posterior density", xlab="Age-4+ abundance in 2011", main="", col="light blue", add=TRUE)
legend("top", c("M=0.130", "M=0.094"), pch=15, col=c("tan", "lightblue"), bty="n")
text("Truth lies somewhere in\nbetween these distributions", x=10500, y=0.002, cex=1.65)
if(.SAVEFIGS) dev.copy2pdf(file=paste(.FIGDIR, "fig:marg_Nt4_estM.pdf", sep=""))







