# ------------------------------------------------------------------------------- #
# R-script for creating summary data for LSMR 
# based on the data-based extraction from GCMRC
# DESCRIPTION OF ROUTINES:
# 	get.DF: reads csv file with raw data,  appends YEAR MONTH TAGNO TL columns
# 	
# 	getGrowthIncrement: extracts initial tag length and most recent recapture 
# 						length information from DF, at-least 1-year at large
# 
# 	getCapturesByLength: extacts capture history by length interval at each time
# 						 step.
# ------------------------------------------------------------------------------- #
require(ggplot2)
require(reshape)
require(Hmisc)
require(PBSmodelling)
dfile    <- "raw2012data.csv"
xbin     <- seq(50, 600, by=10)

# ------------------------------------------------------------------------------- #
# Read data file and create data.frame object for use. Returns (DF)
# ------------------------------------------------------------------------------- #
get.DF   <- function(dfile)
{
	DF       <- read.csv(dfile, header=T)
	t1       <- as.POSIXct(strptime(as.character(DF$START_DATETIME), "%m/%d/%Y"))
	DF$DATE  <- t1
	DF$YEAR  <- as.numeric(format(DF$DATE,"%Y"))
	DF$MONTH <- as.numeric(format(DF$DATE,"%m"))
	DF$TAGNO <- DF$TH_ENCOUNTER_RANKING
	DF$TL    <- DF$TOTAL_LENGTH
	DF$XI    <- xbin[findInterval(DF$TL, xbin)]
	
	# Read in Gear codes and convert to Gear Group
	gearType  <- read.table("gearCode.txt", header=T, sep="\t")
	gid       <- pmatch(DF$GEAR_CODE,gearType[,1], duplicates.ok=TRUE)
	DF$GROUP  <- gearType[gid, 2]
	
	# Sort data frame by date the tagno,  add recapture field
	DF        <- DF[order(DF$DATE, DF$TAGNO), ]
	DF$RECAP  <- duplicated(DF$TAGNO)
	
	
	return (DF)
}

# ------------------------------------------------------------------------------- #
# Construct growth increment data.frame for estimating growth using Laslett model.
# ------------------------------------------------------------------------------- #
getGrowthIncrement <- function(DF)
{
	fnGI     <- function(idNum)
	{
		sDF  <- subset(DF,DF$TAGNO==idNum,select=c(DATE,TL,YEAR,RIVER_CODE))
		if(dim(sDF)[1]==1) return(NULL)
		sDF  <- sDF[order(sDF[,1]),]
		ni   <- dim(sDF)[1]
		dt   <- as.integer(sDF[ni, 1] - sDF[1, 1]) # days at large
		if(dt<=365) return(NULL)                   # at least 1-year at liberty
		l1   <- sDF[1 , 2]                         # length at tagging
		l2   <- sDF[ni, 2]                         # length at recapture
		yr   <- sDF[1 , 3]                         # release year
		tl   <- sDF[1 , 4]                         # release location
		dl   <- (l2-l1)/(dt/365.25)                # annual growth increment
		return (c(yr, tl, l1, l2, dt, round(dl, 2)))
	}
	idNum    <- unique(DF$TAGNO)
	TRdata   <- NULL
	for(i in idNum) TRdata <- rbind(TRdata, fnGI(i))
	colnames(TRdata) <- c("#YEAR","RIVER","L1","L2","DT","DL")
	TRdata <- na.omit(TRdata)     #Remove missing recapture measurements
	ii <- which(TRdata$DL <= -5)  #Remove records with growth increments < -5 mm
	TR$RIVER[TR$RIVER==1] <- "COR"
	TR$RIVER[TR$RIVER==2] <- "LCR"
	return(as.data.frame(TRdata[-ii, ]))
}

# ------------------------------------------------------------------------------- #
# Get captures & recaptures by length interval
# ------------------------------------------------------------------------------- #
getCaptureByLength <- function(DF)
{
	xbin	<- seq(50, 600, by=5)
	nbin	<- length(xbin)
	nyr		<- length(unique(DF$YEAR))
	ctx		<- matrix(0, nrow=nyr,ncol=nbin-1)
	fnCH	<- function(idNum)
	{
		sDF	<- subset(DF, DF$TAGNO==idNum,select=c(YEAR,MONTH,RIVER_CODE,TL))
		sDF <- cbind(sDF, RECAP=FALSE)
		sDF <- sDF[order(sDF$YEAR,sDF$MONTH), ]
		sDF <- sDF[!duplicated(sDF,MARGIN=c(1,2)),]
		
		sDF$RECAP[-1] <- TRUE
		
		sDF$TL	<- factor(sDF$TL, levels=xbin)
		t_TL	<- table(sDF$TL)
		ctx		<- ctx + t_TL
		return(sDF)
	}
}

# ------------------------------------------------------------------------------- #
# Table of monthly captures by year
# ------------------------------------------------------------------------------- #
tableCaptures <- function(DF)
{
	
	# Table of total HBC captures by MONTH~YEAR
	names(DF) <- toupper(names(DF))
	DFm       <- melt(DF, id=c("YEAR","MONTH"), na.rm=F)
	tx        <- cast(DFm,YEAR~MONTH,length,subset=variable=="TL"
	                  ,margins=c("grand_row","grand_col"))
	fn        <- "../../TABLES/LSMR/tableCaptureNumbers.tex"
	cap       <- "Number of fish measured by year and month, sampled by all gears
	              in all reaches of both the LRC and COR."
	cgrp      <- c("", "MONTH")
	ncgrp     <- c(1, 13)
	d1        <- latex(tx,file=fn,rowname=NULL, caption=cap
		         , size="footnotesize",cgroup=cgrp, n.cgroup=ncgrp)
	
}

# ------------------------------------------------------------------------------- #
# Table of annual captures by gear type
# ------------------------------------------------------------------------------- #
tableGear <- function(DF)
{
	names(DF) <- toupper(names(DF))
	DFm       <- melt(DF, id=c("YEAR","GROUP"), na.rm=FALSE)
	tx        <- cast(DFm,YEAR~GROUP,length, subset=variable=="TL"
	                  ,margins=c("grand_row","grand_col"))
	fn        <- "../../TABLES/LSMR/tableCaptureGear.tex"
	cap       <- "Number of fish captured by gear type listed in the GCMRC database
	              for each year."
	cgrp      <- c("", "YEAR")
	ncgrp     <- c(1, 10)
	d1        <- latex(tx,file=fn,rowname=NULL, caption=cap
		         , size="footnotesize",cgroup=cgrp, n.cgroup=ncgrp)
	
}

# ------------------------------------------------------------------------------- #
# Table of annual captures by gear type
# ------------------------------------------------------------------------------- #
tableLf <- function(DF)
{
	names(DF) <- toupper(names(DF))
	DFm       <- melt(DF, id=c("YEAR","XI"), na.rm=TRUE)
	tx        <- cast(DFm,YEAR~XI,length, subset=variable=="TL")
	
	plotBubbles(t(tx[,c(-1,-48)]),xval=1989:2012,yval=sort(unique(DF$XI,na.rm=T))
	            , size=0.15, hide0=TRUE, frange=0.01
	            , xlab="Year", ylab="Length (mm)")
	rect(1988,100,2013,150,col=colr(4,0.2),border=NA)
	fig       <- "../../FIGS/LSMR/fig:CaptureLFbubbles.pdf"
	dev.copy2pdf(file=fig)
	
	fn        <- "../../TABLES/LSMR/tableCaptureLF.tex"
	cap       <- "Number of fish captured by length by all gear types 
	              listed in the GCMRC database for each year."
	cgrp      <- c("", "YEAR")
	ncgrp     <- c(1, 22)
	
	d1        <- latex(t(tx),file=fn, caption=cap
		         , size="tiny",cgroup=cgrp, n.cgroup=ncgrp)
	# use sed to make sideways table
	sed.exp <- paste("sed -i~ 's/table/sidewaystable/g' ", fn)
	system(sed.exp)

}

# ------------------------------------------------------------------------------- #
# Table of new marks released each year (RECAP=FALSE)
# ------------------------------------------------------------------------------- #
tableNewMarks <- function(DF)
{
	names(DF) <- toupper(names(DF))
	DFm       <- melt(DF, id=c("YEAR","MONTH","XI","RECAP","GROUP"), na.rm=TRUE)
	tx        <- cast(DFm,YEAR~XI~RECAP|GROUP, length, subset=c(variable=="TL"))
	gear      <- c("GILL", "HOOP") # Add gear group here to plot results
	
	# Barplots of marks released and recaptured.
	fn.plot <- function(x)
	{
		par(mfcol=c(6, 4), mar=c(0.5, 1, 0.5, 1), oma=c(4, 4, 1, 1))
		n   <- dim(x)[1]
		for(i in 1:n)
		{
			if(!i%%6 || i==n)par(xaxt="s") else par(xaxt="n")
			if(i==n)lgdtxt=c("Mark", "Recap") else lgdtxt=NULL
			barplot(t(x[i, , ]), xlim=c(0, length(xbin))
			        , legend.text=lgdtxt
			        , args.legend=list(bty="n", cex=0.75))
			title( main=rownames(x)[i],   font.main=2, cex.main=0.75, line=-1)
		}
		mtext(c("Length (mm)","Frequency"),side=c(1, 2), outer=T, line=2)
		
		fig=paste("../../FIGS/LSMR/fig:MarksAtLength",gear[jj],".pdf", sep="")
		print(fig)
		jj<<-jj+1
		dev.copy2pdf(file=fig)
	}
	
	event    <- c("Mark", "Recapture")
	cap      <- c("Number of new marks released by year and size interval."
	            , "Number of recaptured marks by year and size interval.")
	fdir     <- "../../TABLES/LSMR/"
	cgrp      <- c("YEAR")
	

	fn.write <- function(x)
	{
		for(i in 1:2)
		{
			ncgrp <- c(dim(x)[1])
			fi    <- paste(fdir,"table:",event[i],":", gear[jj],".tex", sep="")
			LI    <- t(x[, , i])
			d1    <- latex(LI, file=fi, caption=cap[i]
			               , size="tiny", cgroup=cgrp, n.cgroup=ncgrp)
			
			# use sed to replace table with sidewaystable
			sed.exp <- paste('sed -i~ "s/table/sidewaystable/g" ', fi)
			system(sed.exp)
		}
		jj <<- jj+1
		print(jj)
	}
	
	ii = names(tx) %in% gear
	jj        <<-1
	lapply(tx[ii], fn.plot)
	jj       <<- 1
	lapply(tx[ii], fn.write)
}

# ------------------------------------------------------------------------------- #
# Read in data frames and obtain tag-recapture data
# ------------------------------------------------------------------------------- #
if(!exists("DF")) DF <- get.DF(dfile)
if(!exists("TR"))
{
	TR <- getGrowthIncrement(DF)
	TR <- TR[order(TR$YEAR,TR$RIVER),]
	write(dim(TR)[1],file="HBC_Tag_Recapture.dat")
	write(c("#YEAR","RIVER","L1","L2","DT","DL"), ncol=6, 
			file="HBC_Tag_Recapture.dat", append=TRUE)
	write.table(TR,file="HBC_Tag_Recapture.dat",row.names=F, 
			col.names=F, append=TRUE)
	
	fig="../../FIGS/LSMR/fig:GrowthIncrements.pdf"
	plot.gi(TR, file=fig)
}

# ------------------------------------------------------------------------------- #
# Plotting routines
# ------------------------------------------------------------------------------- #
plot.gi <- function(TR, file=NULL,  ...)
{
	p<-ggplot(TR,aes(L1,DL))
	p<-p + geom_point(aes(colour=factor(-RIVER), levels=2),alpha=0.5)
	p<-p + stat_quantile(alpha=0.7, col="black")+facet_wrap(~YEAR) 
	p<-p + labs(x="Release Length (mm)", y="Annual growth increment (mm)")
	p<-p + labs(colour="River")+theme_grey(base_size =  12, base_family = "")

	
	if(!is.null(file))
		dev.copy2pdf(file=file)
	
	return(p)
}