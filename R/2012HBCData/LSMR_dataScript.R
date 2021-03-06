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
source("../read.admb.R")
require(reshape)
require(Hmisc)
require(PBSmodelling)
dfile    <- "raw2012data.csv"
xbin     <- seq(50, 600, by=10)
seas     <- seq(1,12,by=3)

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
    DF$QTR   <- findInterval(DF$MONTH, seas)
    
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
# Create LSMR datafile VPA type
# ------------------------------------------------------------------------------- #
write.LSMRdatafileVPA <- function(DF, ...)
{
	# In this case write annual data of marks-at-large,  recaptures
	# in the following calendar year. It will be necessary to remove
	# duplicates of individuals captured multiple times in each year.
	
	# Remove duplicates by combining YEAR and TAGNO and use duplicated
	# to get boolean index.
	XX  <- as.numeric(paste(DF$YEAR,DF$TAGNO,sep=""))
	ii  <- duplicated(XX)
	iDF <- DF[!ii,]
	
	# Sort data frame by date the tagno,  add recapture field
    iDF        <- iDF[order(iDF$DATE, iDF$TAGNO), ]
    iDF$RECAP  <- duplicated(iDF$TAGNO)
    
	
	# Melt dataframe 
	DFm <- melt(iDF, id=c("YEAR","MONTH","QTR","TRIP_ID","RECAP","XI","GROUP"), na.rm=TRUE)
	
	# Cast DFm to extract catch-at-length,  marks-at-length, recaps-at-length
	C   <- cast(DFm, YEAR ~ XI|GROUP, length, subset=variable=="TL", add.missing=TRUE)
	T   <- cast(DFm, YEAR ~ XI ~ RECAP|GROUP, length, subset=variable=="TL", add.missing=TRUE)
	fn  <- function(x)length(unique(x))
	E   <- cast(DFm[DFm$YEAR!=2012,],YEAR~GROUP,fn, subset=variable=="START_DATETIME", add.missing=TRUE)
	
	# Ensure R matrix and M matrix have the same dimensions at C
	#ii  <- colnames(C$HOOP) %in% colnames(R$HOOP)
	#tmp <- C$HOOP; tmp[, ii] <- R$HOOP; tmp[, !ii] <- 0; R$HOOP <- tmp
	#ii  <- colnames(C$GILL) %in% colnames(R$GILL)
	#tmp <- C$GILL; tmp[, ii] <- R$GILL; tmp[, !ii] <- 0; R$GILL <- tmp
	
	
	# Plot the data to be used
	graphics.off()
	quartz(width=9,height=6.5)
	
	plotSidebars(t(C$HOOP), scale=2, cpro=TRUE); 
	dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Hoop_C_Bar.pdf", width=9, height=6.5)
	plotSidebars(t(T$HOOP[,,1]), scale=2, cpro=TRUE); 
	dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Hoop_M_Bar.pdf", width=9, height=6.5)
	plotSidebars(t(T$HOOP[,,2]), scale=2, cpro=TRUE); 
	dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Hoop_R_Bar.pdf", width=9, height=6.5)
	                                                                                                
	plotSidebars(t(C$GILL), scale=2, cpro=TRUE); 
	dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Gill_C_Bar.pdf", width=9, height=6.5)
	plotSidebars(t(T$GILL[,,1]), scale=2, cpro=TRUE); 
	dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Gill_M_Bar.pdf", width=9, height=6.5)
	plotSidebars(t(T$GILL[,,2]), scale=2, cpro=TRUE); 
	dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Gill_R_Bar.pdf", width=9, height=6.5)
	
	# Write the data file.
	dfn <- "../../ADMB/srcLSMR/HBC2011_ANNUAL.dat"
	gr  <- c("HOOP", "GILL")
	# Get index for gears
	ic  <- names(C) %in% gr
	ie  <- names(E) %in% gr
	write("#Data for HBC 1989:2011", file=dfn)
	write("#syr, nyr, dt", file=dfn, append=TRUE)
	write(c(1989, 2011, 1), file=dfn, append=TRUE)
	write("#Number of gears",file=dfn, append=TRUE)
	write(length(gr), file=dfn, append=TRUE)
	xbin = as.numeric(colnames(C$HOOP)[-1])/10  #units cm
	write("#nbin", file=dfn, append=TRUE)
	write(length(xbin), file=dfn, append=TRUE)
	write("#xbin", file=dfn, append=TRUE)
	write(xbin, file=dfn, append=TRUE)
	
	write("#Array dimensions (rows,  cols) for each gear in (C M R)", file=dfn, append=TRUE)
	write.table(matrix(unlist(lapply(C[ic],dim)),nrow=2,byrow=TRUE),file=dfn,row.names=F, col.names=FALSE, append=TRUE)
	
	# effort data (# of sets)
	write("#Number of sets by gear for each time step", file=dfn, append=TRUE)
	lapply(E[ie], write, file=dfn, append=TRUE)
	
	#order is GILL then HOOP		
	write("#Captures by year (row) and length (col)", file=dfn, append=TRUE)
	lapply(C[ic],write.table,file=dfn,row.names=FALSE,col.names=FALSE, quote=FALSE,append=TRUE)
	
	write("#New Marks by year (row) and length (col)", file=dfn, append=TRUE)
	write.table(T$GILL[, , 1], file=dfn, append=TRUE, quote=FALSE, row.names=TRUE, col.names=FALSE)
	write.table(T$HOOP[, , 1], file=dfn, append=TRUE, quote=FALSE, row.names=TRUE, col.names=FALSE)
	#lapply(T[ic],write.table,file=dfn,row.names=FALSE,col.names=FALSE, quote=FALSE,append=TRUE)
	
	write("#Recaps by year (row) and length (col)", file=dfn, append=TRUE)
	write.table(T$GILL[, , 2], file=dfn, append=TRUE, quote=FALSE, row.names=TRUE, col.names=FALSE)
	write.table(T$HOOP[, , 2], file=dfn, append=TRUE, quote=FALSE, row.names=TRUE, col.names=FALSE)
	#lapply(T[ic][, , 2],write.table,file=dfn,row.names=FALSE,col.names=FALSE, quote=FALSE,append=TRUE)
	
	write("#End of file", file=dfn, append=TRUE)
	write(999, file=dfn, append=TRUE)
	
}


# ------------------------------------------------------------------------------- #
# Create LSMR datafile
# ------------------------------------------------------------------------------- #
write.LSMRdatafile <- function(DF, ...)
{
	dfn <- "../../ADMB/srcLSMR/HBC2011qtr.dat"
	gr  <- c("HOOP", "GILL")
	
	# Melt dataframe 
	DFm <- melt(DF, id=c("YEAR","MONTH","QTR","TRIP_ID","RECAP","XI","GROUP"), na.rm=TRUE)
	
	# Cast DFm to extract Catch-at-length,  Marks-at-length
	C   <- cast(DFm,YEAR + QTR~XI|GROUP,length,subset=variable=="TL", add.missing=TRUE)
	M   <- cast(DFm[DFm$RECAP==FALSE, ], YEAR+QTR~XI|GROUP, length, subset=variable=="TL", add.missing=TRUE)
	R   <- cast(DFm[DFm$RECAP==TRUE, ], YEAR+QTR~XI|GROUP, length, subset=variable=="TL", add.missing=TRUE)
	E   <- cast(DFm, YEAR+QTR~GROUP, function(x) length(unique(x)),subset=variable=="START_DATETIME", add.missing=TRUE)
	E   <- E[1:92, ]
	# Plot the data to be used
	plot.b <- function(X)
	{
		yr <- seq(1989, 2011.75, by=0.25)
		plotBubbles(t(X[,-1:-2]),xval=yr,yval=colnames(X[,-1:-2]), prettyaxis=TRUE
		           , hide0=T,size=0.15,frange=0.01,cpro=FALSE,powr=1
		           , xlab="Year", ylab="Size class (Total Length mm)")
	}
	graphics.off()
	quartz(width=9,height=6.5)
	plot.b(C$HOOP); dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Hoop_C.pdf", width=9, height=6.5)
	plot.b(M$HOOP); dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Hoop_M.pdf", width=9, height=6.5)
	plot.b(R$HOOP); dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Hoop_R.pdf", width=9, height=6.5)
	
	plot.b(C$GILL); dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Gill_C.pdf", width=9, height=6.5)
	plot.b(M$GILL); dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Gill_M.pdf", width=9, height=6.5)
	plot.b(R$GILL); dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Gill_R.pdf", width=9, height=6.5)
	
	# Effort
	barplot(rbind(E$HOOP,E$GILL),beside=T,names.arg=paste(yr)
	       , xlab="Year", ylab="Effort (number of sets)",
	       , legend.text=c("Hoop","Tramel net"), args.legend=list(x="top", bty="n") )
	dev.copy2pdf(file="../../FIGS/LSMR/DATA/fig:Effort.pdf", width=9, height=6.5)
	
	# Get index for gears
	ic  <- names(C) %in% gr
	im  <- names(M) %in% gr
	ir  <- names(R) %in% gr
	ie  <- names(E) %in% gr
	
	# Ensure R matrix and M matrix have the same dimensions at C
	ii  <- colnames(C$HOOP) %in% colnames(R$HOOP)
	tmp <- C$HOOP; tmp[, ii] <- R$HOOP; tmp[, !ii] <- 0; R$HOOP <- tmp
	ii  <- colnames(C$GILL) %in% colnames(R$GILL)
	tmp <- C$GILL; tmp[, ii] <- R$GILL; tmp[, !ii] <- 0; R$GILL <- tmp
	
	
	
	write("#Data for HBC 1989:2011", file=dfn)
	write("#syr, nyr, dt", file=dfn, append=TRUE)
	write(c(1989, 2011, 1/4), file=dfn, append=TRUE)
	write("#Number of gears",file=dfn, append=TRUE)
	write(length(gr), file=dfn, append=TRUE)
	xbin = as.numeric(colnames(C$HOOP)[-1:-2])/10  #units cm
	write("#nbin", file=dfn, append=TRUE)
	write(length(xbin), file=dfn, append=TRUE)
	write("#xbin", file=dfn, append=TRUE)
	write(xbin, file=dfn, append=TRUE)
	
	write("#Array dimensions (rows,  cols) for each gear in (C M R)", file=dfn, append=TRUE)
	write.table(matrix(unlist(lapply(C[ic],dim)),nrow=2,byrow=TRUE),file=dfn,row.names=F, col.names=FALSE, append=TRUE)

	# effort data (# of sets)
	write("#Number of sets by gear for each time step", file=dfn, append=TRUE)
	lapply(E[ie], write, file=dfn, append=TRUE)

	#order is GILL then HOOP		
	write("#Captures by year (row) and length (col)", file=dfn, append=TRUE)
	lapply(C[ic],write.table,file=dfn,row.names=F,col.names=F, quote=F,append=T)

	write("#Marks by year (row) and length (col)", file=dfn, append=TRUE)
	lapply(M[im],write.table,file=dfn,row.names=F,col.names=F, quote=F,append=T)
	
	write("#Recaps by year (row) and length (col)", file=dfn, append=TRUE)
	lapply(R[ir],write.table,file=dfn,row.names=F,col.names=F, quote=F,append=T)
	
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

plotSidebars <- function(z, scale=1, cpro=FALSE, frange=0.05, ...)
{
	plot.new()
	# set up graphical parameters
	par(mar=c(1, 1, 1, 1), oma=c(1, 1, 1, 1)*4, xaxs='r')
	dz  <- dim(z)
	dn  <- dimnames(z)
	if(is.character(names(dn)))
	{
		lbl <- rev(names(dn))
	}
	else
	{
		lbl <- c("YEAR", "LENGTH")
	}
	
	
	ny = nr = dz[1]
	nx = nc = dz[2]
	y  = as.numeric(unlist(dn[1]))
	x  = as.numeric(unlist(dn[2]))
	yl = extendrange(y, f=frange)
	xl = extendrange(x, f=frange)
	
	par(usr=c(xl, yl))
	
	if(cpro)
	{
		zs = apply(z, 2, sum, na.rm=TRUE)
		zz = sweep(z, 2, zs, "/")
	}
	else
	{
		zz = z
	}
	
	box()
	mtext(lbl, side=1:2, outer=TRUE, line=2)
	axis(1); axis(2)
	scale = scale/diff(range(zz,na.rm=TRUE))
	for( i in 1:nx )
	{
		dy = 0.5*diff(y[1:2])
		xx = as.vector(rbind(y-dy, y+dy))
		xx = c(min(xx), xx, max(xx))
		yy = as.vector(rbind(zz[,i], zz[,i]))*scale
		yy = c(0, yy, 0)
		polygon(x[i]-yy, xx,col=colr("grey", 0.5))
	}
}


PBSplotSideBars <- function (z, xval = FALSE, yval = FALSE, dnam = FALSE, rpro = FALSE, 
    cpro = FALSE, rres = FALSE, cres = FALSE, powr = 0.5, size = 0.2, 
    lwd = 1, clrs = c("black", "red", "blue"), hide0 = FALSE, 
    frange = 0.1, prettyaxis = FALSE, ...)
{
    if (is.data.frame(z)) {
        use = !sapply(z, is.factor) & sapply(z, is.numeric)
        z = z[, use, drop = FALSE]
        if (ncol(z) == 0) {
            showAlert("data frame not useable")
            return()
        }
        z = as.matrix(z)
    }
    dz <- dim(z)
    ny = ny1 = dz[1]
    nx = nx1 = dz[2]
    if (length(dz) > 2) {
        showAlert("Input matrix must have only 2 dimensions")
        return()
    }
    xval1 <- 1:nx
    yval1 <- 1:ny
    if (mode(xval) == "logical") {
        if (xval[1]) {
            xval1 <- z[1, ]
            ny1 <- ny - 1
        }
    }
    if (mode(yval) == "logical") {
        if (yval[1]) {
            yval1 <- z[, 1]
            nx1 <- nx - 1
        }
    }
    xind <- (nx - nx1 + 1):nx
    x2 = xlabel = xval1[xind]
    yind <- (ny - ny1 + 1):ny
    y2 = ylabel = yval1[yind]
    if ((mode(xval) != "logical") & (length(xval) == nx1)) {
        if (mode(xval) == "numeric") 
            x2 = xval
        xlabel = xval
    }
    if ((mode(yval) != "logical") & (length(yval) == ny1)) {
        if (mode(yval) == "numeric") 
            y2 = yval
        ylabel = yval
    }
    zz <- array(z[yind, xind], dim = c(length(yind), length(xind)), 
        dimnames = dimnames(z[yind, xind, drop = FALSE]))
    dots = list(...)
    xlab = dots$xlab
    if (is.null(xlab)) 
        xlab = ""
    ylab = dots$ylab
    if (is.null(ylab)) 
        ylab = ""
    if (dnam & !is.null(dimnames(zz))) {
        warn = options()$warn
        options(warn = -1)
        if (!is.null(dimnames(zz)[[2]])) {
            xpos = try(as.numeric(dimnames(zz)[[2]]), silent = TRUE)
            if (all(is.na(xpos))) 
                xlabel = dimnames(zz)[[2]]
            else if (!any(is.na(xpos)) && all(diff(xpos) > 0 | 
                all(diff(xpos) < 0))) {
                xlabel = as.character(xpos)
                x2 = xpos
            }
        }
        if (!is.null(dimnames(zz)[[1]])) {
            ypos = try(as.numeric(dimnames(zz)[[1]]), silent = TRUE)
            if (all(is.na(ypos))) 
                ylabel = dimnames(zz)[[2]]
            else if (!any(is.na(ypos)) && all(diff(ypos) > 0 | 
                all(diff(ypos) < 0))) {
                ylabel = as.character(ypos)
                y2 = ypos
            }
        }
        options(warn = warn)
    }
    xx <- rep(x2, each = length(y2))
    yy <- rep(y2, length(x2))
    minz <- min(zz, na.rm = TRUE)
    maxz <- max(zz, na.rm = TRUE)
    if (rpro | cpro) {
        if (minz < 0) {
            zz <- zz - minz
            minz <- 0
            maxz <- max(zz, na.rm = TRUE)
        }
    }
    if (rpro) {
        zs <- apply(zz, 1, sum, na.rm = TRUE)
        zz <- sweep(zz, 1, zs, "/")
    }
    if (cpro) {
        zs <- apply(zz, 2, sum, na.rm = TRUE)
        zz <- sweep(zz, 2, zs, "/")
    }
    if (rres) {
        zm <- apply(zz, 1, mean, na.rm = TRUE)
        zz <- sweep(zz, 1, zm, "-")
    }
    if (cres) {
        zm <- apply(zz, 2, mean, na.rm = TRUE)
        zz <- sweep(zz, 2, zm, "-")
    }
    zNA <- is.na(zz) | is.nan(zz) | is.infinite(zz)
    zz[zNA] <- 0
    z0 <- sign(zz) * abs(zz)^abs(powr)
    z1 <- z3 <- z0
    z1[z0 <= 0] <- NA
    z3[z0 < 0 | z0 > 0] <- NA
    z2 <- -z0
    z2[z0 >= 0] <- NA
    za <- max(z0, na.rm = TRUE)
    zb <- min(z0, na.rm = TRUE)
    zM <- max(abs(z0))
    sz1 <- max(za * size/zM, 0.001)
    sz2 <- max(-zb * size/zM, 0.001)
    evalCall(plot, argu = list(x = 0, y = 0, xlim = extendrange(x2, 
        f = frange), ylim = extendrange(y2, f = frange), type = "n", 
        axes = FALSE, xlab = xlab, ylab = ylab), ..., checkdef = TRUE, 
        checkpar = TRUE)
    if (prettyaxis) {
        if (length(min(x2):max(x2)) <= 5) 
            xshow = is.element(x2, x2)
        else xshow = is.element(x2, pretty(x2, n = 10))
        yshow = is.element(y2, pretty(y2, n = 10))
    }
    else {
        xshow = rep(TRUE, length(x2))
        yshow = rep(TRUE, length(y2))
    }
    if (!all(xshow)) 
        axis(1, at = x2[!xshow], labels = FALSE, tcl = ifelse(is.null(dots$tcl), 
            par()$tcl, dots$tcl)/3)
    if (!all(yshow)) 
        axis(2, at = y2[!yshow], labels = FALSE, tcl = ifelse(is.null(dots$tcl), 
            par()$tcl, dots$tcl)/3)
    evalCall(axis, argu = list(side = 1, at = x2[xshow], labels = xlabel[xshow]), 
        ..., checkpar = TRUE)
    evalCall(axis, argu = list(side = 2, at = y2[yshow], labels = ylabel[yshow]), 
        ..., checkpar = TRUE)
    if (!hide0 && !all(is.na(z3))) {
        evalCall(symbols, argu = list(x = xx, y = yy, circles = as.vector(z3), 
            inches = 0.001, fg = clrs[3], lwd = lwd, add = TRUE), 
            ..., checkpar = TRUE)
    }
    if (!all(is.na(z2))) {
        evalCall(symbols, argu = list(x = xx, y = yy, circles = as.vector(z2), 
            inches = sz2, fg = clrs[2], lwd = lwd, add = TRUE), 
            ..., checkpar = TRUE)
    }
    if (!all(is.na(z1))) {
        evalCall(symbols, argu = list(x = xx, y = yy, squares = as.vector(z1), 
            inches = sz1, fg = clrs[1], lwd = lwd, add = TRUE), 
            ..., checkpar = TRUE)
    }
    box()
    invisible(z0)
}