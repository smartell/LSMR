# R-Script for LASLETT model output.
#require(hacks)
source("../../R/read.admb.R")
A <- read.admb("LASLETT")
class(A) <- "laslett"

# S3 class
plot.laslett <- function(A, ...)
{
	# Creat summer plot for model fit.
	opar <- par(no.readonly=TRUE)
	par(mfcol=c(2, 2))
	with(A, {
		# PDF of age at tagging
		h  <- hist(age,breaks=1:max(age+1),
			xlab="Age at tagging (years)", ylab="pdf", freq=FALSE)
		p1 <- log(muA); p2 <- sigA
		xl <- range(age)
		curve(unlist(lapply(x,dlnorm,p1,p2)),xl[1],xl[2], add=T, col=2, lwd=2)
		
		# Distribution of linf
		h <- hist(linf, xlab="Asymptotic length (mm)", ylab="pdf", freq=FALSE)
		p1 <- mu_linf; p2 <- sd_linf
		xl <- range(linf)
		curve(unlist(lapply(x,dnorm,p1,p2)),xl[1],xl[2], add=T, col=2, lwd=2)
		
		# length-at-age (tagging)
		iage <- seq(0, max(age+dt))
		ilen <- mu_linf*(1-exp(-k*(iage)))
		plot(iage, ilen,type="l",ylim=c(0,1.3*max(ilen)), 
			xlab="Age (years)", ylab="Length (mm)", lwd=5)
		
		segments(age, l1, age+dt, l2, col="grey")
		points(age,l1,pch=20, col="cyan")
		
		# Residuals
		xx = cbind(1:length(l1_res),1:length(l1_res)+0.25 )
		matplot(xx, cbind(l1_res, l2_res),type="h",  pch=".",col=c(1, 2), lty=1
			, xlab="Individual", ylab="Residual")
		
	})
	par(opar)
}