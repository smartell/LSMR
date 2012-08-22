# Generic plots for LSMR
# In general, for each function an lsmr model object,  or a list of objects
# is pased as an argument,  where each object represents the entire contents
# of a unique model.  If a list of objects exists, the each element is either
# plotted on a unique plot,  or overlayed on a single plot (.OVERLAY==TRUE).

# Pseudcode for a typical plot:
# 1) Determine the number of objects
# 2) Construct a data frame with Model being the id variable.
# 3) Use ggplot2 to facet or overlay the plots.


.lsmrPlotNt <- function(M, what="Nt", ylbl=NULL, ...)
{
	n  <- length(M);
	# construct data.frame
	mdf <- NULL
	for(i in 1:n)
	{
		idx <- match(what, names(M[[i]]))
		print(idx)
		print(M[[i]][idx])
		df	<- data.frame(Year=M[[i]]$yr, Model=names(M)[i], Nt=M[[i]][idx])
		mdf <- rbind(mdf, melt(df, id.vars=c("Model", "Year")))
	}
	names(mdf) <- c("Model", "Year", "Variable", "Abundance")
	print(mdf)
	
	if(is.null(ylbl))
	{
		ylbl="Value"
	}
	
	
	p <- ggplot(mdf, aes(x=Year, y=Abundance, col=Model)) 
	p <- p + geom_line()
	p <- p + scale_y_continuous(ylbl, limits=c(0,max(mdf$Abundance)))
	if( ! .OVERLAY )
	{
		p <- p + facet_wrap(~ Model, scales="free")
	}
	print(p)
}



# plot Numbers greater than 50 mm 
.DEP.lsmrPlotNt <- function(M,  ...)
{
	
	n <- length(M)
	for(im in 1:n)
	{
		with(M[[im]], {
		
			yl <- c(0, 1.2*max(Nt))
			if(exists("true_Nt"))
				yl <- c(0, 1.2*max(Nt, true_Nt))
			
			plot(yr, Nt, type="l", ylim=yl
			, xlab="Year", ylab="Abundance (> 50 mm)", ...)
			
			grid()
			
			if(exists("true_Nt"))
			{
				lines(yr, true_Nt, lwd=2, col=colr(2, 0.5))
				legend("top", c("estimated","true"),lty=1, lwd=c(1, 2)
				, col=c(1, colr(2, 0.5)), bty="n")
			}
		
			# Add credible interval
			if(exists("Nt.ps"))
			{
				ci <- apply(Nt.ps, 2, quantile, probs=c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
				.polygon(yr, ci[c(1, 7), ])
				.polygon(yr, ci[c(2, 6), ])
				.polygon(yr, ci[c(3, 5), ])
			
			}
		 
		})
	}	
}