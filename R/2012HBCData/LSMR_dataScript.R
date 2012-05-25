# R-script for creating summary data for LSMR based on the data-based extraction from GCMRC
require(ggplot2)
dfile    <- "raw2012data.csv"
# ------------------------------------------------------------------------------- #
# Read data file and create data.frame object for use. Returns (DF)
# ------------------------------------------------------------------------------- #
get.DF   <- function(dfile)
{
	DF       <- read.csv(dfile, header=T)
	t1       <- as.POSIXct(strptime(as.character(DF$START_DATETIME), "%m/%d/%Y"))
	DF$DATE  <- t1
	DF$YEAR  <- as.numeric(format(DF$DATE,"%Y"))
	DF$MONTH <- as.numeric(format(DF$DATE,"%M"))
	DF$TAGNO <- DF$TH_ENCOUNTER_RANKING
	DF$TL    <- DF$TOTAL_LENGTH
	
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
		return (c(yr, tl, l1, l2, dt, dl))
	}
	idNum    <- unique(DF$TAGNO)
	TRdata   <- NULL
	for(i in idNum) TRdata <- rbind(TRdata, fnGI(i))
	colnames(TRdata) <- c("#YEAR","RIVER","L1","L2","DT","DL")
	TR <- na.omit(TR)           #Remove missing recapture measurements
	ii <- which(TR$DL <= -5)  #Remove records with growth increments < -5 mm
	return(as.data.frame(TRdata[-ii, ]))
}

if(!exists("DF")) DF <- get.DF(dfile)
if(!exists("TR"))
{
	TR <- getGrowthIncrement(DF)
	write(dim(TR)[1],file="HBC_Tag_Recapture.dat")
	write.table(TR,file="HBC_Tag_Recapture.dat",row.names=F, append=TRUE)
}

p<-ggplot(TR,aes(L1,DL))
p+geom_point(aes(colour=factor(-RIVER)),alpha=0.5)+stat_quantile(alpha=0.7, col="black")+facet_wrap(~YEAR) + labs(x="Release Length (mm)", y="Annual growth increment (mm)")+labs(colour="River")+theme_grey(base_size =  12, base_family = "")

dev.copy2pdf(file="../../FIGS/LSMR/fig:GrowthIncrements.pdf")
