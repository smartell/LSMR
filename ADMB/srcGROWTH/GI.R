require(PBSmodelling)
require(ggplot2)
source("../../R/Read.admb.R")
A <- read.admb("GI")
class(A) <- c("list", "GI")

plot.GI <- function(A, ...)
{
	with(A, {
		par(mfcol=c(6, 4), mar=c(2, 2, 1, 1), oma=c(3, 3, 1, 1))
		for(i in 1:nobs)
		{
			y <- dl[i, ]/dt[i, ]
			plot(l1[i, ], y, xlab="", ylab="", col=loc[i, ]
			, xlim=c(min(l1, na.rm=T), max(400, na.rm=T))
			, ylim=c(min(dl/dt, na.rm=T), max(dl/dt, na.rm=T)))
			lines(sort(l1[i, ]), rev(sort(dl_hat[i, ]/dt[i, ])), pch=19, col=4)
			abline(h=0)
			title(main=iyr[i, 1])
		}
		
		par(mfcol=c(6, 4), mar=c(2, 2, 1, 1), oma=c(3, 3, 1, 1))
		fn <- function (x)
		{
			qqnorm(x, main="", xlab="", ylab="")
			qqline(x)
		}
		apply(epsilon, 1, fn)
		
	})
}



plot.marg <- function(A, ...)
{
	par(mfcol=c(2, 1))
	with(A, {
		linf= as.data.frame(exp(post.samp[, 1]+post.samp[, 7:27]))
		colnames(linf)<-paste(1990:2010)
		mdf <- melt(linf)
		q<-qplot(factor(variable),value,data=mdf,geom="violin", xlab="Year", ylab="Asymptotic length (mm)")
		print(q)
		dev.copy2pdf(file="../../FIGS/LSMR/fig:LinfPosteriors.pdf", height=5, width=10)
	})
	
	
	with(A, {
		k = as.data.frame(exp(post.samp[, 2]+post.samp[, 28:48]))
		colnames(k)<-paste(1990:2010)
		mdf <- melt(k)
		q<-qplot(factor(variable),value,data=mdf,geom="violin", xlab="Year", ylab="Metabolic rate paraemter (k)")
		print(q)
		dev.copy2pdf(file="../../FIGS/LSMR/fig:kPosteriors.pdf", height=5, width=10 )	
	})
	
}


plot.TransitionMatrix <- function(A, ... )
{
	L       = read.rep("TransitionMatrix.txt")
	names(L)=c("lbins",paste("T",1990:2010,sep=""))
	
	xx = yy = mean(L$lbins[1:2])+cumsum(diff(L$lbins))
	zz = t(L$T2005)
	iyr  <- 1990:2010
	linf <- A$linf
	jj  <<-1
	fn<-function(zz)
	{
		zz = t(zz)
		colnames(zz) = xx
		rownames(zz) = yy
		plotBubbles(zz,hide0=TRUE,cpro=FALSE, rpro=FALSE,prettyaxis=FALSE,size=0.05,frange=0.01,dnam=TRUE,
			 clrs=colr(4, 0.5), xlab="", ylab="", yaxt="n", xaxt="n", powr=1.0)
		title(main=iyr[jj], line=-1)
		abline(a=0, b=1, col=colr(1, 0.5), lty=1)
		points(linf[jj], linf[jj], pch=19, col=colr(2, 0.5))
		grid()
		jj <<- jj+1
		mfg <- par(no.readonly=TRUE)$mfg
		if(mfg[2]==1 && mfg[1]%%2) axis(2)
		if(mfg[2]==nc && !mfg[1]%%2) axis(4)
		if(mfg[1]==nr && mfg[2]%%2) axis(1)
		if(mfg[1]==1 && !mfg[2]%%2) axis(3)
		
	}
	nc<-4
	nr<-6
	par(bg="white")
	par(mfcol=c(nr,nc),mar=c(0,0,0,0),oma=c(4,4,4,4), cex.lab=0.5, mfg=c(1, 1, 6, 4), las=1)
	lapply(L[-1], fn)
	mtext("Total length (mm) in year t", 1, outer=TRUE, line=2.5)
	mtext("Total length (mm) in year t+1", 2, outer=TRUE, line=2.5, las=0)
	dev.copy2pdf(file="../../FIGS/LSMR/fig:TransitionMatrix.pdf")	
}

