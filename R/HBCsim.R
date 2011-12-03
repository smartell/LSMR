## ___________________________________________________________ ##
## Simulation model for humpback chub                          ##
## Author: Steven Martell                                      ##
## DATE: September 26,  2011                                   ##
## Contact: martell.steve@gmail.com                            ##
##                                                             ##
## FUNCTIONS:                                                  ##
##   - .ibm => given N recruits returns list with the capture  ##
##             history of each recruit. It also fills the      ##
##             global variables C M and R for the assessment   ##
##             model.                                          ##
##   - .runIBM Loops over syr:nyr and applies the .ibm function##
##             to each of Rt recruits.                         ##
##                                                             ##
##                                                             ##
## ___________________________________________________________ ##


## __ Required libraries _____________________________________ ##
require( PBSmodelling )
require( Hmisc )
require( ggplot2 )


## __ Simulation controls ____________________________________ ##

#  Model Dimensions
syr     <- 1950     #initial year
nyr     <- 2011     #final year
dt      <- 3/12     #time step
tmax    <- 50/dt    #max longevity index (50 yrs)
dyr     <- seq(syr, nyr, by=dt)
xbin    <- seq(3, 55, by=1)

#  Growth (units in cm)
lmin        <- 4
lmin.sig    <- 0.4
m.linf      <- 0.08
linf        <- 40   #Asymptotic length (cm)
k           <- 0.18 #Growth coefficient

#  Selectivity
lh          <- 10   #length @ 50% selectivity
gamma       <- 1.5  #std in selectivity

sample.yrs  <- seq(1989, 2011, by=dt)
fishing     <- dyr %in% sample.yrs
min.size    <- 15.0 #minimum size for tagging (15 cm)



#  Parameter vector for the IBM model
THETA <- list( lmin=lmin, 
               lmin.sig=lmin.sig, 
               m.linf=m.linf, 
               lh=lh, 
               gamma=gamma )

set.seed(991)
iclr <- rev(topo.colors(length(syr:nyr),0.5))
Rt<-floor(rlnorm(length(syr:nyr), log(300), 0.9))
Et<-runif(length(sample.yrs), 0.05, 0.2)  #Effort

#  Global objects for storing data
fish.id     <- 0    #Unique id for individual fish
tag.id      <- 0    #Unique tag no. for individual
#  Total catch, marked, recap, unmarked by year at length interval
C   <- matrix(0,nrow=length(sample.yrs), ncol=length(xbin))
M   <- matrix(0,nrow=length(sample.yrs), ncol=length(xbin))
R   <- matrix(0,nrow=length(sample.yrs), ncol=length(xbin))

## ___________________________________________________________ ##





## ___________________________________________________________ ##
## PROCEDURES                                                  ##
## ___________________________________________________________ ##


.ibm <- function( Nt, jyr, THETA )
{
    #ARGUMENTS:
        #Nt    = Number of recruits from cohort
        #jyr   = brood year (year the cohort was age-0)
        #THETA = a list of simulation parameters for the ibm.
    
    #RETURNS:
        #df    = a list w records of capture date, len,  tag.no

    with(THETA, {
        #jyr <- syr+(ii-1)
        ## indexs i = individual; j=time

        ## Calculate length & mortality for time j.
        lj <- linf*(1.-exp(-k*dt*1:tmax))
        mj <- m.linf*linf/lj
        
        ## Draw initial recruit size (gamma).
        scale   <- (lmin.sig)^2/lmin
        shape   <- lmin/scale
        li <- rgamma(Nt, shape, 1/scale)
        
        
        ## Growth trajectory
        growth <- function ( li )
        {
            linf.i <- rnorm(1, linf, 0.1*linf)
            li+(linf.i-li)*(1-exp(-k*1:tmax*dt))
        }
        L <- sapply(li, growth)
        
        ## Survival trajectory
        ## S is a list of lengths of individuals at
        ## the time they were alive.
        survival <- function ( lj )
        {
            sj <- exp(-m.linf*linf/lj*dt)
            pj <- rbinom(tmax, 1, sj)
            nj <- which.min(pj)
            return(lj[1:nj])
        }
        S <- apply(L, 2, survival)
        

        ## Capture history
        pcap <- function ( lj )
        {
            #Arg: lj is the length at times j
            #Algorithm:
            #-1. determine index for sampling time periods (fyr)
            #-2. determine capture probability at sampling times
            #-3. if captured fish > 150 mm then assign a tag no.
            if(!is.null(lj))
            {
                # 1. Index for sampling time periods (fyr)
                n   <- length(lj)
                d1  <- which(dyr==jyr)  #min index of dyr for jyr
                d2  <- min(d1+n-1, length(dyr))
                jj  <- d1:d2
                iyr <- dyr[jj]          #years alive
                fyr <- fishing[jj]      #years fished fyr==1
                
                
                ## Did not survive into sampling years
                if(length(iyr[fyr])==0) return(NULL)
                
                # -2. Selectivity at time of sampling
                len <- (lj[1:length(jj)])[fyr]
                iyr <- iyr[fyr]
                # TODO Add annual effort to sj capture probability
                eyr <- findInterval(iyr, sample.yrs)
                #print(cbind(sample.yrs[eyr],Et[eyr]))
                sj  <- Et[eyr]*plogis(len, lh, gamma)
                
                cj  <- rbinom(length(sj), 1, sj)
                rid <- which(cj==1)  #row index for capture
                if(length(rid)==0) return(NULL) #never captured
                
                tmp <- c(iyr[rid], len[rid],
                    tag.no=rep(NA,length=length(rid)))
                tmp <- matrix(tmp, ncol=3)
                

                # -3. Assign tag.id if greater than 150 mm
                if(max(tmp[,2])>=min.size)
                {
                    tag.id <<- tag.id + 1
                }
                tmp[len[rid]>=min.size, 3] <- tag.id
                
                # -4. Populate global variables with catch history
                tt <- findInterval(tmp[,1],sample.yrs)
                xx <- findInterval(tmp[,2],xbin)
                for(i in 1:length(tt))
                    C[tt[i], xx[i]] <<- C[tt[i], xx[i]] + 1  

                #print(cbind(tmp,sample.yrs[tt],xbin[xx],tmp[,3]))
                
                # - Newly marked fish into matrix M
                im   <- !is.na(tmp[,3])
                if(sum(im)>0) #fish was large enough to tag
                {
                    it <- min(which(im))
                    iy <- findInterval(tmp[it,1],sample.yrs)
                    ix <- findInterval(tmp[it,2],xbin)

                    #print(cbind(tmp, sample.yrs[iy],xbin[ix]))
                    M[iy, ix] <<- M[iy, ix] + 1
                    
                    ir <- which(im)[-1]
                    ry <- findInterval(tmp[ir,1],sample.yrs)
                    rx <- findInterval(tmp[ir,2],xbin)
                    for(i in 1:length(ry))
                        R[ry[i], rx[i]] <<- R[ry[i], rx[i]] + 1
                }
                
                return(tmp)
            }
        }
        P <- sapply(S, pcap)
        return(P)
    })
}

.runIBM <- function()
{
    j <- 1
    A <- NULL
    for(i in syr:nyr)
    {
        A <- c(A, .ibm(Rt[j], i, THETA))
        j <- j+1
    }
    return(A)
}

.writeLSMRdata <- function()
{
	# This function writes the data file for LSMR model
	fn <- "simLSMR.dat"
	C  <- cbind(floor(sample.yrs),sample.yrs%%1*4+1,C)
	M  <- cbind(floor(sample.yrs),sample.yrs%%1*4+1,M)
	R  <- cbind(floor(sample.yrs),sample.yrs%%1*4+1,R)
	
	write("#Simulated data from HBCsim.R", fn)
	write("#Model Dimensions", fn, append=T)
	write("#Years, timestep", fn, append=T)
	write(c(range(sample.yrs), dt), fn, append=T)
	write("#Array dimensions(C, M, R)", fn, append=T)
	write(dim(C), fn, append=T)
	write("#Number of length intervals", fn, append=T)
	nbin=length(xbin)+1
	write(nbin, fn, append=T)
	write("#Length intervals (cm)", fn, append=T)
	write(c(xbin,xbin[nbin-1]+diff(xbin[(nbin-2):(nbin-1)])), fn, append=T)
	
	write("#Catch at length by period", fn, append=T)
	write("#Year, Period, Count ...", fn, append=T)
	write.table(C, fn, append=T, col.names=F, row.names=F)
	write("#Number Marked at length by period", fn, append=T)
	write("#Year, Period, Count ...", fn, append=T)
	write.table(M, fn, append=T, col.names=F, row.names=F)
	write("#Number Recaptured at length by period", fn, append=T)
	write("#Year, Period, Count ...", fn, append=T)
	write.table(R, fn, append=T, col.names=F, row.names=F)
	
	write("#End of file", fn, append=T)
	write(999, fn, append=T)
	
}

## ___________________________________________________________ ##
## MAIN                                                        ##
## ___________________________________________________________ ##
#A  <- .ibm(1000, 1980, THETA)
A <- .runIBM()
.writeLSMRdata()

## ___________________________________________________________ ##

par(mfcol=c(6, 4), las=1, mar=c(2, 2, 1, 1), oma=c(3, 3, 1, 1))
for(i in seq(1,length(sample.yrs)-1/dt,by=1/dt))
{
	
	ut <- colSums(C[i:(i+3),]-M[i:(i+3),]-R[i:(i+3),])
	mt <- colSums(M[i:(i+3),])
	rt <- colSums(R[i:(i+3),])
	tmp <- rbind(ut, mt, rt)
	
    barplot(tmp,names.arg=xbin,xlab=""
    ,ylab="", main=sample.yrs[i])
}
mtext(c("Length (cm)","Frequency"), 
    c(1, 2), outer=T, las=0, line=1)