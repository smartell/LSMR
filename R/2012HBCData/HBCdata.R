#R Script for processing ASMR data files based on data queries 
#provided by Paul Alley.  Also,  Lew Coggins was instramental
#in setting up this script.
setwd("~/Documents/CONSULTING/HumpbackChub/HBC_2011_Assessment/R/2012HBCData")
fname= "raw2012data"#paste("Scenario",1:7, sep="")
require(Hmisc)


#if(!exists("HBCdata"))
get.df=function(fn="~/Documents/CONSULTING/HBC ASMR Performance/HBC Data/GCMRC_ASMR_2009_DATASETS/Scenario1.csv")
{
	HBCdata = read.csv(fn, header=T)
	df=HBCdata
	print(paste("Total number of records =",dim(df)[1]))

	#Fix the Date column to proper m/d/Y format
	x=as.character(df$START_DATETIME)
	x=as.POSIXct(strptime(x,"%m/%d/%Y"))
	df$DATES=x
	df$TAGNO = df$TH_ENCOUNTER_RANKING
	df$TL    = df$TOTAL_LENGTH
	
	#Reset global objects that go into the data file
	nyr<<-length(1989:2012)
	mta<<-matrix(0, nrow=nyr, ncol=49)
	rta<<-matrix(0, nrow=nyr, ncol=49)
	rcta<<-array(0, c(nyr, nyr, 49))
	
	return(df)
}



fn<-function(tagno)
{
	#Subset the data frame for tagno
	sdf=cbind(subset(df, df$TAGNO == tagno), 
		RECAP=F, YR=factor(NA, levels=1989:2012), AGE=NA)
	sdf=sdf[order(sdf$TAGNO, sdf$DATE), ]
	
	#Get year for each record
	sdf$YR=getyr(sdf$DATE, sdf)
	
	#Get age at time of tagging
	sdf$AGE[1] = getage(sdf$TL[1])
	iage=sdf$YR-sdf$YR[1]
	sdf$AGE=sdf$AGE[1]+iage
	
	#Recaptures & remove duplicates in same recapture YR
	sdf$RECAP[-1]=T
	sdf=sdf[!duplicated(sdf$YR), ]
	
	#New marks released
	mdf=sdf[sdf$RECAP==F, ]
	mdf$AGE=factor(mdf$AGE, levels=2:50)
	mdf$YR=factor(mdf$YR, levels=1989:2012)
	mta<<-mta+t(as.matrix(table(mdf[, 15:14])))
	
	#Recaptured marks
	rdf=sdf[sdf$RECAP==T, ]
	rdf$AGE=factor(rdf$AGE, levels=2:50)
	rdf$YR=factor(rdf$YR, levels=1989:2012)
	r1 = t(as.matrix(table(rdf[, 15:14])))
	rta<<-rta+r1
	
	#Recaptures by tag year
	iyr=sdf$YR[1]-1989+1
	rcta[iyr, , ]<<-rcta[iyr, , ]+r1
}

fn.m<-function(tagno)
{
	#This function does the same operations as fn, 
	#but restricts the data to Spring or Fall period only.
	## RULES:
	##	IGNORE NO LONGER IMPLEMENTED -1 fish must be tagged in the period
	##	IGNORE NO LONGER IMPLEMENTED -2 ignore recaptures that are not in the period.
	##	-Spring months:	Jan-Jun
	##	-Fall months:	Jul-Dec
	
	#Subset the data frame for tagno
	sdf=cbind(subset(df, df$TAGNO == tagno), 
		RECAP=F, YR=factor(NA, levels=1989:2012), AGE=NA, Month=NA)
	sdf=sdf[order(sdf$TAGNO, sdf$DATE), ]
	
	#Get year and month for each record
	sdf$YR=getyr(sdf$DATE, sdf)
	sdf$Month=getmth(sdf$DATE)
	#print(sdf)
	#if(sdf$Month[1]>6)return()	#Rule #1
	#sdf=sdf[sdf$Month<=6, ]	 	#Rule #2
	#print(sdf)
	
	#Get age at time of tagging
	sdf$AGE[1] = getage(sdf$TL[1])
	iage=sdf$YR-sdf$YR[1]
	sdf$AGE=sdf$AGE[1]+iage
	
	#Recaptures & remove duplicates in same recapture YR
	sdf$RECAP[-1]=T
	sdf=sdf[!duplicated(sdf$YR), ]
	
	#New marks released
	mdf=sdf[sdf$RECAP==F, ]
	mdf$AGE=factor(mdf$AGE, levels=2:50)
	mdf$YR=factor(mdf$YR, levels=1989:2012)
	mta<<-mta+t(as.matrix(table(mdf[, 9:8])))
	
	#Recaptured marks
	rdf=sdf[sdf$RECAP==T, ]
	rdf$AGE=factor(rdf$AGE, levels=2:50)
	rdf$YR=factor(rdf$YR, levels=1989:2012)
	r1 = t(as.matrix(table(rdf[, 9:8])))
	rta<<-rta+r1
	
	#Recaptures by tag year
	iyr=sdf$YR[1]-1989+1
	rcta[iyr, , ]<<-rcta[iyr, , ]+r1
}



getage=function(len)
{
	alast=30
	## this function does new growth curve  (LEW COGGINS pers comm.)
	ifelse(round(-1.05952E-08*len^4+1.25697E-05*len^3-0.005058694*len^2
		+0.878738367*len^1-53.55077549)>alast,alast,
		round(-1.05952E-08*len^4+1.25697E-05*len^3-0.005058694*len^2
			+0.878738367*len^1-53.55077549))
	
	
	ifelse(round(-1.05952E-08*len^4+1.25697E-05*len^3-0.005058694*len^2
		+0.878738367*len^1-53.55077549)==1,2,round(-1.05952E-08*len^4
			+1.25697E-05*len^3-0.005058694*len^2
			+0.878738367*len^1-53.55077549))
}

getyr=function(datee, df){
     ifelse(strptime(datee,"%Y")$year+1900==1992 | 
			strptime(datee,"%Y")$year+1900==1996 |
            strptime(datee,"%Y")$year+1900==2000 | 
			strptime(datee,"%Y")$year+1900==2004 | 
			strptime(datee,"%Y")$year+1900==2008,
            strptime(df$DATE-7862400,"%Y")$year+1900,
            strptime(df$DATE-7776000,"%Y")$year+1900)
}

getmth=function(datee){
	return(as.numeric(format(datee,"%m")))
}

#Write temporary output data file
write.data.file=function(fn="HBCdata.output")
{
	
	write("#mta",file=fn)
	write.table(mta, file=fn, row.names=F, col.names=F, append=T)
	write("#rta", file=fn, append=T)
	write.table(rta, file=fn, row.names=F, col.names=F, append=T)
	yrs=1989:2012; j=0
	for(i in yrs)
	{
		j=j+1
		write(paste("#rcta", i), file=fn, append=T)
		write.table(rcta[j, , ], file=fn, row.names=F, 
			col.names=F, append=T)
	}
	write("#eof", file=fn, append=T)
	write(999, file=fn, append=T)
}
month=c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")


table.1<-function()
{
	#A latex table for the total number of records by month and year
	m=c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec", "Total")
	d=data.frame(date=df$DATES)
	d$yr=as.numeric(format(d$date,"%Y"))
	d$mth=as.numeric(format(d$date, "%m"))
	xx=table(d[, 2:3])
	xx=cbind(xx, rowSums(xx))
	xx=rbind(xx, Total=colSums(xx))
	cap="Total number of records in the database (new marks and recaptures) by year and month.
	Note that the numbers here include multiple recaptures of the same individual in the same
	year."
	latex(xx, title="Table.1", rowlabel="Year", colheads=m, caption=cap, here=T)
	return(xx)
}

table.2<-function()
{
	#Latex table for the unmarked percentage
	d=data.frame(date=df$DATES, id=df$TAGNO)
	d=d[order(d$id, d$date), ]
	
	d$recap=FALSE
	comptags=function(index){d$id[index]==d$id[index-1]}
	d$recap[c(FALSE,sapply(2:dim(d)[1],comptags))]=TRUE
	
	
	d$yr=as.numeric(format(d$date,"%Y"))
	d$mth=as.numeric(format(d$date, "%m"))
	m=d[d$recap==FALSE, ] #new marks
	mm=table(m[, 4:5])
	xx=100*mm/table(d[, 4:5])
	xx=cbind(xx, "mark rate" =100-rowMeans(xx, na.rm=T))
	xx=round(rbind(xx, Average=colMeans(xx, na.rm=T)), 1)
	cap="Percentage of number of records that are new marks released by year and month.
	 The mark rate column is the average mark rate calculated as proportion of the total number
	of fish caught that were previously marked."
	latex(xx, title="Table.2", rowlabel="Year", colheads=c(month, "mark rate"), caption=cap, here=T)
	return(xx)
}

## ******************************
##           MAIN
## ******************************

main<-function()
{
	dfile=paste(fname, ".csv", sep="")
	#for(j in 1	:length(dfile)){
		#df=get.df(dfile[j])
		df=get.df(dfile)
		tag.ids=unique(df$TAGNO)
		t1=Sys.time()
		#for(i in tag.ids[1:1000]) fn(i)

		for(i in tag.ids) fn(i)

		print(Sys.time()-t1)
		ofile=paste(fname, ".dat", sep="")
		write.data.file(ofile)
		rm(df)
	#}
	
}
