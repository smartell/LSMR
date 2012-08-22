.plotMarginalPosteriors	<- function( admbObj )
{
	print("	.plotMarginalPosteriors")
	#Marginal distributions & priors for theta
	op=par(no.readonly=T)
	par(mfcol=c(3, 2), oma=c(1, 1, 1, 1)*4, mar=c(1, 4, 1, 1), las=0)

	## Read control file to get bounds and priors for theta
	## ctrl=read.table(A$control.file, header=F, skip=13, nrow=6)

	with(admbObj, {
		#mcmc = post.samp
		#colnames(mcmc) = fit$names[1:dim(mcmc)[2]]
		std=apply(mcmc[,1:8],2,sd)
		nr=length(std[std!=0])
		
		for(i in 1:8){
			if(std[i]!=0){
				ps = mcmc[, i]  #posterior samples
				xl=c(ctrl[i, 2], ctrl[i, 3])#range(ps)
				xl = range(ps)
				
				
				#browser()
				hist(ps,xlab=colnames(mcmc[i]),prob=T, 
					main="", ylab="", col="lightgrey", 
					xlim=xl)#, ...)
				
				
				
				## Add priors
				nfn=c("dunif","dnorm","dlnorm","dbeta","dgamma")
				pt = ctrl[i, 5]+1
				fn=match.fun(nfn[pt])
				p1=ctrl[i, 6]; p2=ctrl[i, 7]
				#browser()
				if(pt!=4){
					curve(unlist(lapply(x,fn,p1,p2)),
						xl[1],xl[2],add=TRUE, col=colr(4, 0.7), lwd=2)
					if(pt!=1){
						sr <- round(sd(ps)/p2, 2)
						if(pt==5)sr <- round(sd(ps)/sqrt(p1*(1/p2)^2), 2)
						legend("topright", paste(sr), bty="n")
					}
				}
				else{
					curve(unlist(lapply((x-ctrl[i,2])/
						 (ctrl[i,3]-ctrl[i,2])
						,fn,p1,p2)),xl[1],xl[2],add=TRUE, col=colr(4, 0.7), lwd=2)
					sr <- round(sd(ps)/sqrt( p1*p2/((p1+p2)^2 *(p1+p2+1)) ), 2)
					legend("topright", paste(sr), bty="n")
				}
						
				
			}
		}

		mtext(c("Parameter", "Probability density"), c(1, 2), outer=TRUE, line=c(2,2), las=0)
	})
	par(op)
}
