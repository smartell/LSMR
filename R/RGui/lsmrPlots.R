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
	print(".lsmrPlotNt")
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

.lsmrPlotViolin <- function(M, what="Nt.ps", ylbl=NULL, ...)
{
	# Violin plots of the marginal posterior distributions.
	# construct data.frame
	mdf <- NULL
	n   <- length(M)
	print(n)
	for(i in 1:n)
	{
		idx <- match(what, names(M[[i]]))
		N   <- t(data.frame(M[[i]][idx]))
		print(N)
		df  <- data.frame(Year=M[[i]]$yr, Model=names(M)[i], Nt=N )
		mdf <- rbind(mdf, melt(df, id.vars=c("Model", "Year")))
		print(mdf)
	}
	names(mdf) <- c("Model", "Year", "Variable", "Abundance")
	
	if(is.null(ylbl))
	{
		ylbl="Value"
	}
	
	p <- ggplot(mdf, aes(factor(Year), Abundance, col=Model))
	p <- p + geom_violin()
	if( ! .OVERLAY )
	{
		p <- p + facet_wrap(~ Model, scales="free")
	}
	print(p)
	
}