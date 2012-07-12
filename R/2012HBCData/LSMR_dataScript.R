# ------------------------------------------------------------------------------- #
# R-script for creating summary data for LSMR 
# based on the data-based extraction from GCMRC
# Author: Steven Martell
# Date: May,  2012
# LOCATION: MAUI
# 
# DESCRIPTION OF ROUTINES:
#   get.DF: reads csv file with raw data,  appends YEAR MONTH TAGNO TL columns
#   
#   getGrowthIncrement: extracts initial tag length and most recent recapture 
#                       length information from DF, at-least 1-year at large
# 
#   getCapturesByLength: extacts capture history by length interval at each time
#                        step.
#   
#   tableCaptures: extracts the monthly captures each year
# 
#   tableGear: extacts the number of HBC capture each year by gear type.
# 
#   tableLf: extracts the catch at length (10 mm bins) each year.
#   
#   tableMarks: marks released and recaptured by GILL and HOOP nets
# ------------------------------------------------------------------------------- #
require(ggplot2)
# require(hacks) #deprecated
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
        if(dt<=365 || dt>730 ) return(NULL)        # at least 1-year at liberty
        l1   <- sDF[1 , 2]                         # length at tagging
        l2   <- sDF[ni, 2]                         # length at recapture
        yr   <- sDF[1 , 3]                         # release year
        tl   <- sDF[1 , 4]                         # release location
        dl   <- (l2-l1)/(dt/365.25)                # annual growth increment
        return (c(yr, tl, l1, l2, dt, round(dl, 2)))
    }
    idNum    <- unique(DF$TAGNO)
    TRdata   <- data.frame()
    for(i in idNum) TRdata <- rbind(TRdata, fnGI(i))
    colnames(TRdata) <- c("YEAR","RIVER","L1","L2","DT","DL")
    
    TRdata <- na.omit(TRdata)      #Remove missing recapture measurements
    ii <- which(TRdata$DL <= -10)  #Remove records with growth increments < -5 mm
    #TRdata$RIVER[TRdata$RIVER==1] <- "COR"
    #TRdata$RIVER[TRdata$RIVER==2] <- "LCR"
    return(as.data.frame(TRdata[-ii, ]))
}

# ------------------------------------------------------------------------------- #
# Annual growth increment data based on captures and recaps in subsequent year.
# ------------------------------------------------------------------------------- #
annualGrowthIncrement <- function(DF)
{
    # Goal: for each year compute the length-at-capture & growth increment into
    #       the following year.
    # PSEUDOCODE:
    # - get vector of unique years.
    # - find unique individuals recaptured year i and year i+1
    # - get TL from each capture event
    
    ATR    <- data.frame()
    iyr    <- unique(DF$YEAR)
    for(i in iyr)
    {
        iDF  <- rbind(subset(DF,DF$YEAR==i), subset(DF, DF$YEAR==i+1))
        iDF  <- iDF[order(iDF$TAGNO, iDF$DATE), ]
        irc  <- unique(iDF$TAGNO[duplicated(iDF$TAGNO)])
        
        for(j in irc)
        {
            ijDF <- subset(iDF, iDF$TAGNO==j)
            nj   <- dim(ijDF)[1]
            dt   <- as.integer(ijDF$DATE[nj] - ijDF$DATE[1])
            l1   <- ijDF$TL[1]
            l2   <- ijDF$TL[nj]
            loc  <- ijDF$RIVER_CODE[1]
            ATR  <- rbind(ATR, c(YEAR=i, RIVER=loc, TAGNO=j, l1=l1, l2=l2, dt=dt))
        }
    }
    colnames(ATR)=c("YEAR","RIVER","TAGNO","l1","l2","dt")
    return(ATR)
}


# ------------------------------------------------------------------------------- #
# Get captures & recaptures by length interval
# ------------------------------------------------------------------------------- #
getCaptureByLength <- function(DF)
{
    xbin    <- seq(50, 600, by=5)
    nbin    <- length(xbin)
    nyr     <- length(unique(DF$YEAR))
    ctx     <- matrix(0, nrow=nyr,ncol=nbin-1)
    fnCH    <- function(idNum)
    {
        sDF <- subset(DF, DF$TAGNO==idNum,select=c(YEAR,MONTH,RIVER_CODE,TL))
        sDF <- cbind(sDF, RECAP=FALSE)
        sDF <- sDF[order(sDF$YEAR,sDF$MONTH), ]
        sDF <- sDF[!duplicated(sDF,MARGIN=c(1,2)),]
        
        sDF$RECAP[-1] <- TRUE
        
        sDF$TL  <- factor(sDF$TL, levels=xbin)
        t_TL    <- table(sDF$TL)
        ctx     <- ctx + t_TL
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
                 , size="footnotesize",cgroup=cgrp, n.cgroup=ncgrp
                 , label="table:Captures")
    
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
    cap       <- "Number of fish captured by gear type listed in the"
    cap       <- paste(cap, " GCMRC database for each year.")
    cgrp      <- c("", "YEAR")
    ncgrp     <- c(1, 10)
    d1        <- latex(tx,file=fn,rowname=NULL, caption=cap
                 , size="footnotesize",cgroup=cgrp, n.cgroup=ncgrp
                 , label="table:Gear")
    
}

# ------------------------------------------------------------------------------- #
# Table of number of TRIP_IDs and CPUE per year
# ------------------------------------------------------------------------------- #
tableEffort <- function(DF)
{
    names(DF) <- toupper(names(DF))
    DFm       <- melt(DF, id=c("YEAR","GROUP"), na.rm=FALSE)
    tx        <- cast(DFm,YEAR~.,function(x) length(unique(x))
                      ,subset=variable=="TRIP_ID",margins=TRUE )
    td        <- cast(DFm, YEAR~., function(x) length(unique(x))
                      ,subset=variable=="DATE",margins=TRUE )
    tn        <- cast(DFm, YEAR~GROUP, function(x) length(unique(x))
                      ,subset=variable=="START_DATETIME"
                      ,margins=TRUE )
    tl        <- cast(DFm, YEAR~GROUP, length, subset=variable=="TL"
                      ,margins=TRUE)

    # Effort table
    Effort    <- cbind(tx, td[,-1],tn[,c(10, 6, 5)]
                       ,tn[,10]-rowSums(tn[,c(6, 5)]) , tl[,c(6,5)]
                       ,round(tl[,c(6,5)]/tn[,c(6,5)], 2))
    colnames(Effort) <- c("Year","Trips","Days","Sets","Hoop"
                       ,"Tramel","Other","Hoop Ct","Tramel Ct"
                       ,"Hoop CPUE","Tramel CPUE")

    fn        <- "../../TABLES/LSMR/tableEffort.tex"
    cap       <- "Number of trips,  days fished,  unique sets, 
                  hoop and tramel net sets,  other gears,  catch of humpback
                  chub in hoop nets and tramel nets,  and the corresponding
                  arithmatic CPUE for each gear. Note that these CPUE trends should
                  not be used as a relative abundance index because zero catch of HBC 
                  have been excluded from the effort data."
    d1        <-latex(Effort[-25,], file=fn, caption=cap
                , label="table:Effort", size="scriptsize", rowname=NULL)
}


# ------------------------------------------------------------------------------- #
# Table of annual captures by gear type
# ------------------------------------------------------------------------------- #
tableLf <- function(DF)
{
    names(DF) <- toupper(names(DF))
    DFm       <- melt(DF, id=c("YEAR","XI","GROUP"), na.rm=TRUE)
	tx        <- cast(DFm,YEAR~XI,length, subset=variable=="TL")
    rtx       <- cast(DFm,YEAR~XI|GROUP,length, subset=variable=="TL")
    gear      <- c("GILL", "HOOP") # Add gear group here to plot results

    plotBubbles(t(tx[,c(-1,-48)]),xval=1989:2011,yval=sort(unique(DF$XI,na.rm=T))
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
    LI        <- t(tx)
    d1        <- latex(LI,file=fn, caption=cap
                 , size="tiny",cgroup=cgrp, n.cgroup=ncgrp
                 , lable="table:Lf")
    # use sed to make sideways table
    sed.exp <- paste("sed -i~ 's/table/sidewaystable/g' ", fn)
    system(sed.exp)

	return(rtx)
}

# ------------------------------------------------------------------------------- #
# Table of marks released and recaptured each year (RECAP=FALSE)
# ------------------------------------------------------------------------------- #
tableMarks <- function(DF)
{
    names(DF) <- toupper(names(DF))
    DFm       <- melt(DF, id=c("YEAR","MONTH","XI","RECAP","GROUP"), na.rm=TRUE)
    tx        <- cast(DFm,YEAR~XI~RECAP|GROUP, length, subset=c(variable=="TL"), fill=0)
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
    cap      <- c("Number of new marks released by year and size interval "
                , "Number of recaptured marks by year and size interval ")
    fdir     <- "../../TABLES/LSMR/"
    cgrp      <- c("YEAR")
    

    fn.write <- function(x)
    {
        for(i in 1:2)
        {
            ncgrp <- c(dim(x)[1])
            pcap  <- paste(cap[i], "(", gear[jj], ").", sep="")
            fi    <- paste(fdir,"table:",event[i],":", gear[jj],".tex", sep="")
            LI    <- t(x[, , i])
            d1    <- latex(LI, file=fi, caption=pcap, size="tiny"
                     , cgroup=cgrp, n.cgroup=ncgrp
                     , label=paste("table:",event[i],"_",gear[jj], sep="")
                     )
            
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
	
	return(tx)
}
# ------------------------------------------------------------------------------- #
# Create LSMR datafile
# ------------------------------------------------------------------------------- #
write.LSMRdatafile <- function()
{
	dfn <- "../../ADMB/srcLSMR/HBC2011.dat"
	gr  <- c("HOOP", "GILL")
	
	C   <- tableLf(DF)        # Total number of captures C[year, length]
	MR  <- tableMarks(DF)     # Total number of markes released and recaptured.
	ic  <- names(C) %in% gr
	im  <- names(MR) %in% gr
	
	write("#Data for HBC 1989:2011", file=dfn)
	write("#syr, nyr, dt", file=dfn, append=TRUE)
	write(c(1989, 2011, 1), file=dfn, append=TRUE)
	write("#Number of gears",file=dfn, append=TRUE)
	write(length(gr), file=dfn, append=TRUE)
	xbin = seq(50, 500, by=10)
	write("#nbin", file=dfn, append=TRUE)
	write(length(xbin), file=dfn, append=TRUE)
	write("#xbin", file=dfn, append=TRUE)
	write(xbin, file=dfn, append=TRUE)
	
	write("#Array dimensions(C)", file=dfn, append=TRUE)
	write.table(matrix(unlist(lapply(C[ic],dim)),nrow=2,byrow=TRUE),file=dfn,row.names=F, col.names=FALSE, append=TRUE)
	write("#Array dimensions(MR)", file=dfn, append=TRUE)
	write.table(matrix(unlist(lapply(MR[im],dim)),nrow=2,byrow=TRUE),file=dfn,row.names=F, col.names=FALSE,  append=TRUE)
	
	write("#Length Intervals for (C)", file=dfn, append=TRUE)
	x = lapply(C[ic],colnames)
	for(i in 1:length(x))
	  write(as.numeric(x[[i]][-1]), file=dfn, append=TRUE)
	write("#Length Intervals for (MR)", file=dfn, append=TRUE)
	write(as.numeric(colnames(MR$GILL[,,1])), file=dfn, append=TRUE)
	write(as.numeric(colnames(MR$HOOP[,,1])), file=dfn, append=TRUE)
	
	write("#Captures by year (row) and length (col)", file=dfn, append=TRUE)
	lapply(C[ic],write.table,file=dfn,row.names=F,col.names=F,append=T)
	
	write("#Marks and Recaptures by year (row) and length (col)", file=dfn, append=TRUE)
	lapply(MR[im],write.table,file=dfn,row.names=T,col.names=F,append=T, quote=FALSE)
	
	write("#End of file", file=dfn, append=TRUE)
	write(999, file=dfn, append=TRUE)
	
}


# ------------------------------------------------------------------------------- #
# Read in data frames and obtain tag-recapture data
# ------------------------------------------------------------------------------- #
if(!exists("DF")) DF <- get.DF(dfile)
if(!exists("TR"))
{
    TR <- getGrowthIncrement(DF)
    print(head(TR))
    TR <- TR[order(TR$YEAR,TR$RIVER),]
    write(dim(TR)[1],file="HBC_Tag_RecaptureII.dat")
    write(c("#YEAR","RIVER","L1","L2","DT","DL"), ncol=6, 
            file="HBC_Tag_RecaptureII.dat", append=TRUE)
    write.table(TR,file="HBC_Tag_RecaptureII.dat",row.names=F, 
            col.names=F, append=TRUE)
    
    fig="../../FIGS/LSMR/fig:GrowthIncrements.pdf"
    plot.gi(TR, file=fig)
}
if(!exists("ATR"))
{
    ATR  <- annualGrowthIncrement(DF)
    O    <- na.omit(ATR[ATR$dt>365, -3])
    O$dl <- round((O$l2-O$l1)/(O$dt/365.25), 2)
    
    print(head(O))
    fn  <- "HBC_Annual_GI.dat"
    write(dim(O)[1], file=fn)
    write(c("#YEAR","RIVER","l1","l2","dt","dl"), ncol=6, file=fn, append=TRUE)
    write.table(O, file=fn, row.names=FALSE, col.names=FALSE, append=TRUE)
    
    fig <-"../../FIGS/LSMR/fig:AnnualGrowthIncrements.pdf"
    plot.atr(O, file=fig)
}

# ------------------------------------------------------------------------------- #
# Plotting routines
# ------------------------------------------------------------------------------- #
plot.gi <- function(TR, file=NULL,  ...)
{
    p<-ggplot(TR,aes(L1,DL))
    p<-p + geom_point(aes(colour=factor(RIVER), levels=2),alpha=0.5)
    p<-p + stat_quantile(alpha=0.7, col="black")+facet_wrap(~YEAR) 
    p<-p + labs(x="Release Length (mm)", y="Annual growth increment (mm)")
    p<-p + labs(colour="River")+theme_grey(base_size =  12, base_family = "")
    if(!is.null(file))
        dev.copy2pdf(file=file)
    
    return(p)
}

plot.atr <- function(ATR, file=NULL, ...)
{
    O    <- ATR[ATR$dt>365, ]
    O$gi <- (O$l2-O$l1)/(O$dt/365.25)
    p    <- ggplot(O, aes(l1, gi))
    p    <- p + geom_point(aes(color=factor(RIVER), levels=2), alpha=0.5)
    p    <- p + stat_quantile(alpha=0.7, col="black") + facet_wrap(~YEAR)
    p    <- p + labs(colour="River")+theme_grey(base_size =  12, base_family = "")
    p
    if(!is.null(file))
        dev.copy2pdf(file=file)
    
    return(p)
    
}