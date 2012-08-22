.mcmcTrace	<- function( admbObj, label=NULL )
{
	## this function examines the trace plots for the
	## estimated leading parameters
	print("	.mcmcTrace")
	par(mfcol=c(4, 3), oma=c(1, 1, 1, 1)*4, mar=c(1, 4, 1, 1), las=0)
	
	plotTrace <- function( obj )
	{
	  # Input "obj" is a VECTOR of MCMC samples.
	  # Produces one panel trace plot.

	  nSample <- length( obj )
	  plot( c(1:nSample), obj, type="n", axes=FALSE, xlab="", ylab="" )
	  points( c(1:nSample),obj, cex=0.5, pch=20, col="darkgray" )

	  lines( lowess( c(1:nSample),obj,f=1/4), lty=1, lwd=1 )
	  abline( h=mean(obj), lty=2 )

	  # Plot MPD point (1st element).
	  points( 1,obj[1], cex=1.0, pch=16, col="green" )
	  points( 1,obj[1], cex=1.0, pch=1 )    

	  axis( side=1 )
	  axis( side=2 )
	  box()
	}
  
	with(admbObj, {
	  # Find the active parameters.  If the chain is all equal, then the parameter
	  # was fixed in the model configuration.  This gets a Boolean vector that
	  # indicates which columns have fixed values.
	  mcmcObj=post.samp[, 1:12]
	  colnames(mcmcObj) <- admbObj$fit$names[1:dim(mcmcObj)[2]]
	  iPars <- apply( mcmcObj,2,function(x) { sum(diff(x))!=0.0 } )
	  nPars <- sum( iPars )     # Number of active parameters in mcmc output.

	  tmp <- mcmcObj[ ,iPars ]
	  tmpNames <- colnames( tmp )
	
	  for ( i in 1:ncol(tmp) )
	  {
	    plotTrace( tmp[,i] )
		title(ylab=tmpNames[i], line=2.25)
		print(tmpNames[i])
	  }

		mtext(c("Sample", "Parameter"), c(1, 2), outer=TRUE, line=c(2,2), las=0)
	  
	})
}
