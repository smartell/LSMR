## Estimation of growth parameters for HBC using Hillary,  2011 paper.
source("../R/read.admb.R")
gdata <- read.table("HBC.GI.data.txt", header=F)
A <- read.admb("HillaryGrowth")

plot(A$ltag,A$dl/A$dt, pch=20, cex=0.5, xlab="Length at release", ylab="Relative growth increment", xlim=c(0, 500))
abline(h=0)

x=seq(0, 500, by=10)
#linf = 390.878; k=0.183129
#lines(x, (linf-x)*(1-exp(-k)), col=2)
linf = A$linf; k=A$k
lines(x, (linf-x)*(1-exp(-k)), col=4)

points(A$ltag,(A$dl_hat+A$epsilon)/A$dt, pch=1, cex=0.25, col=3)

plot(A$ltag,A$epsilon, type="h")

## Extract the HBC data from the large CSV file.
get.HBC.data<-function()
{
HBC <- read.csv(file = "GCMRC_ASMR_TAG_HISTORY_ALL.csv", header = TRUE, stringsAsFactors = FALSE)

## Reformat the date column
dates <- as.POSIXct(strptime(as.character(HBC$DATES), "%m/%d/%Y"))
HBC$DATES <- dates

id <- unique(HBC$TAGNO)

growth.increment <- function(tag.no)
{
	df <- subset(HBC, TAGNO==tag.no)
	if(dim(df)[1]==1) return(NULL)
	df <- df[order(df$DATE), ]

	n <- dim(df)[1]
	dt <- as.numeric(diff(c(df$DATE[1],df$DATE[n])))/365.25
	if(dt<=0.5) return(NULL)
	dl <- df$TL[n] - df$TL[1]
	return(c(df$TL[1], dl, dt))
}

dlist <- sapply(id, growth.increment)
gdata <- matrix(unlist(dlist),ncol=3,byrow=T)
plot(gdata[,1],gdata[,2]/gdata[,3], pch=".")
abline(h=0)
write.table(gdata,file="HBC.GI.data.txt",row.names=F,col.names=c("RL","dl","dt"))
}