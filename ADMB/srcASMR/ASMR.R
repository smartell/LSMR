#ASMR.R
#R Code for ASMR results.
#SET UP GRAPHICS DEVICE FOR PLOTTING (page-up page-down to scroll through figs)
#graphics.off()
#if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
#windows(record=T); 
#par(las=1,ps=12,cex.main=14/12)
quartz()
require(Hmisc)

## ------------------------------------------------------------------------------
## Documentation for ASMRsturgeon Model:
## Authors:  Steve Martell, Lew Coggins and general whippings by CJW.
## Date:  Many years ago to today
## Instructions:  Source this R code and make sure the ASRMsturgeon1.exe
##			    is in the same directory as this R-code.
##			   Note that you have to modify the flags in the data files to 
##			   control the likelihood and ASMR1-3 parameterizations
##
## Description of functions:
##	-run.model() 	-calls the ASMR program to carry out the estimation.
##			-use the mcmc=T flag to run the mcmc stuff.
##	-get.results()  -returns a list object with the ASMR output from ASMR.
##	-plot.summary() -produces the summary plots
##	-run(df="suwsturgdat.22.dat")	-calls run.model -get.results -plot.summary
## ------------------------------------------------------------------------------

#Source R utilities
setwd("~/Documents/CONSULTING/HBC ASMR Performance/ASRM 2.0")
source("Rutils.R")
fn <- "ASMR.rep"  #Default file name for results of ASMR.
mcon <- T

A = reptoRlist(fn)
B = reptoRlist("~/Documents/CONSULTING/HBC ASMR Performance/Lew's Code & Data/asmr3t.rep")

plot.prt<-function()
{
	#quartz(width=9,height=5)
	par(mfcol=c(1, 2), cex.lab=1.5)
	#This function used the data from the simulation 
	#trials to plot estimated recruitment from 2009:2015
	true.rt=rep(2000, length=7)
	true.rt[3]=4000
	plot(2009:2015, true.rt, type="l", lwd=3, xlab="Year", ylab="Age-2 recruits", ylim=c(0, 5000), xlim=c(2009,2016))
	
	tt=read.table("Trials24.txt", header=F)
	colnames(tt)=c("Scenario", "Seed", "pcap",2009:2015)
	#subset scenario 1
	r1=subset(tt, subset=Scenario==1)
	p1=subset(r1,pcap==0.05, select=c(-Scenario, -Seed, -pcap))
	boxplot(p1,  add=T, at=2009:2015, boxwex=1/10, pch=".")
	px=1:5/5*0.25
	for(i in 2:5){
		
		p2=subset(r1,pcap==px[i], select=c(-Scenario, -Seed, -pcap))
		print(px[i])
		print(head(p2))
		boxplot(p2,  add=T, at=2009:2015+(i-1)/8, boxwex=1/10, col=i, xaxt="n", yaxt="n", pch=".")	
	}
	#legend("topright",paste("Capture probability", px), fill=c(0, 2:5), bty="n")
	#dev.copy2eps(file="../Report/fig.rt1.eps")
	## ------------------------------------------------------------------------------------
	## ------------------------------------------------------------------------------------
	true.rt=rep(2000, length=7)
	true.rt[3]=1000
	plot(2009:2015, true.rt, type="l", lwd=3, xlab="Year", ylab="Age-2 recruits", ylim=c(0, 5000), xlim=c(2009,2016))
	
	tt=read.table("Trials21.txt", header=F)
	colnames(tt)=c("Scenario", "Seed", "pcap",2009:2015)
	#subset scenario 1
	r1=subset(tt, subset=Scenario==1)
	p1=subset(r1,pcap==0.05, select=c(-Scenario, -Seed, -pcap))
	boxplot(p1,  add=T, at=2009:2015, boxwex=1/10, pch=".")
	px=1:5/5*0.25
	for(i in 2:5){
		
		p2=subset(r1,pcap==px[i], select=c(-Scenario, -Seed, -pcap))
		print(px[i])
		print(head(p2))
		boxplot(p2,  add=T, at=2009:2015+(i-1)/8, boxwex=1/10, col=i, xaxt="n", yaxt="n", pch=".")	
	}
	legend("topright",paste("Capture probability", px), fill=c(0, 2:5), bty="n")
	dev.copy2eps(file="../Report/fig.rt1.eps", width=10, height=5)
	
}

plot.rrt<-function()
{
	#This function used the data from the simulation 
	#trials to plot estimated recruitment from 2009:2015
	true.rt=rep(2000, length=7)
	true.rt[3]=1000
	plot(2009:2015, true.rt, type="l", lwd=3, xlab="Year", ylab="Age-2 recruits", ylim=c(0, 5000), xlim=c(2009,2016))
	
	tt=read.table("Trials.txt", header=F)
	colnames(tt)=c("Scenario", "Seed", "pcap",2009:2015)
	#subset scenario 1
	r1=subset(tt, subset=Scenario==1)
	p1=subset(r1,pcap==0.05, select=c(-Scenario, -Seed, -pcap))
	boxplot(p1,  add=T, at=2009:2015, boxwex=1/10, pch=".")
	px=1:5/5*0.25
	for(i in 2:5){
		
		p2=subset(r1,pcap==px[i], select=c(-Scenario, -Seed, -pcap))
		print(px[i])
		print(head(p2))
		boxplot(p2,  add=T, at=2009:2015+(i-1)/8, boxwex=1/10, col=i, xaxt="n", yaxt="n", pch=".")	
	}
	legend("topright",paste("Capture probability", px), fill=c(0, 2:5), bty="n")
	dev.copy2eps(file="../Report/fig.rt3.eps")
}

plot.srt<-function()
{
	#quartz(w=950,h=528)
	par(mfcol=c(1, 2), cex.lab=1.5)
	#This function used the data from the simulation 
	#trials to plot estimated recruitment from 2009:2015
	true.rt=rep(2000, length=7)
	true.rt[3]=4000
	plot(2009:2015, true.rt, type="l", lwd=3, xlab="Year", ylab="Age-2 recruits", ylim=c(0, 5000), xlim=c(2009,2016))
	
	tt=read.table("Trials24.txt", header=F)
	colnames(tt)=c("Scenario", "Seed", "pcap",2009:2015)
	#subset scenario 1
	r1=subset(tt, subset=pcap==0.2)
	p1=subset(r1,Scenario==1, select=c(-Scenario, -Seed, -pcap))
	boxplot(p1,  add=T, at=2009:2015, boxwex=1/10, pch=".")
	px=1:5/5*0.25
	for(i in 2:7){
		
		p2=subset(r1,Scenario==i, select=c(-Scenario, -Seed, -pcap))
		print(px[i])
		print(head(p2))
		boxplot(p2,  add=T, at=2009:2015+(i-1)/8, boxwex=1/10, col=i, xaxt="n", yaxt="n", pch=".")	
	}
	
	#This function used the data from the simulation 
	#trials to plot estimated recruitment from 2009:2015
	true.rt=rep(2000, length=7)
	true.rt[3]=1000
	plot(2009:2015, true.rt, type="l", lwd=3, xlab="Year", ylab="Age-2 recruits", ylim=c(0, 5000), xlim=c(2009,2016))
	
	tt=read.table("Trials21.txt", header=F)
	colnames(tt)=c("Scenario", "Seed", "pcap",2009:2015)
	#subset scenario 1
	r1=subset(tt, subset=pcap==0.2)
	p1=subset(r1,Scenario==1, select=c(-Scenario, -Seed, -pcap))
	boxplot(p1,  add=T, at=2009:2015, boxwex=1/10, pch=".")
	px=1:5/5*0.25
	for(i in 2:7){
		
		p2=subset(r1,Scenario==i, select=c(-Scenario, -Seed, -pcap))
		print(px[i])
		print(head(p2))
		boxplot(p2,  add=T, at=2009:2015+(i-1)/8, boxwex=1/10, col=i, xaxt="n", yaxt="n", pch=".")	
	}
	
	
	legend("topright",paste("Scenario", 1:7), fill=c(0, 2:7), bty="n")
	dev.copy2eps(file="../Report/fig.rt2.eps", width=10, height=5)
}


compare<-function()
{   
	par(mfrow=c(3, 3), las=1)
	plot(A$yr, A$nt4/1000, type="l", xlab="Year", ylab="Adults (age 4+)")
	title(main="Martell code")
	plot(A$yr, B$adults/1000, type="l", xlab="Year", ylab="Adults (age 4+)")
	title(main="Coggins code")
	plot(A$yr, B$adults-A$nt4, type="h", ylim=c(-1, 1), xlab="Year",ylab="Difference")
	title(main="Difference (Coggins-Martell)")
	
	plot(A$byr, A$rt/1000, type="l", xlab="Brood year", ylab="Recruits (age-2)")
	plot(A$byr, B$rt2/1000, type="l", xlab="Brood year", ylab="Recruits (age-2)")
	plot(A$byr, B$rt2-A$rt, type="h", ylim=c(-1, 1),xlab="Year",  ylab="Difference")
	
	matplot(A$yr, A$pta, type="l", xlab="Brood year", ylab="Capture probability", col="grey")
	matplot(A$yr, B$ploc, type="l", xlab="Brood year", ylab="Capture probability", col="grey")
	matplot(A$yr, B$ploc-A$pta, type="h", ylim=c(-1, 1), xlab="Year",ylab="Difference")
	
	dev.copy2eps(file="Comparison.eps")
	
	
	#Now compare retrospective estimates
	#Martell Model
	SM = run.retro("S1.dat", retyr=9)
	LC = retro(retyr=9)
	
	par(mfrow=c(2, 2), las=1)
	matplot(A$yr, SM$nt4/1000, type="l", xlab="Year", ylab="Adults (1000s age-4+)")
	title(main="Martell code") 
	matplot(A$yr, LC$nt4/1000, type="l", xlab="Year", ylab="Adults (1000s age-4+)")
	title(main="Coggins code")
	
	matplot(A$byr, SM$rt/1000, type="l", xlab="Brood year", ylab="Recruits (1000s age-2)") 
	matplot(B$byr, LC$rt/1000, type="l", xlab="Brood year", ylab="Recruits (1000s age-2)")
	dev.copy2eps(file="RetroComparison.eps")
}                     


#some code to plot up the life-history age-schedule information
life.history <- function()
{
	linf = 2177; k=0.13; m =0.08
	age=1:50
	x = seq(1, linf)
	len=linf*(1-exp(-k*age))
	vx = plogis(x, 500, 700)
	ax = -1./k*log(1-x/linf)	#age from length
	t1 = exp(k*(age+1.))-1.;
	t2 = exp(k*(age))-1.;
	sa  = (t1/t2)^(-m/k) #value(pow(elem_div(t1,t2),-m/vonbk));
	par(mfcol=c(2, 2), las=1)
	plot(age, len, xlab="Age (years)", ylab="Length (mm)", ylim=c(0, linf), type="l")
	plot(linf*(1-exp(-k*ax)), as.integer(ax), pch=".", xlab="Length (mm)", ylab="Age (years)")
	plot(ax, vx, xlab="Age (years)", ylab="Vulnerability to gear", type="l")
	plot(age, sa, col=2, xlab="Age (years)", ylab="Survival rate", type="l")
	
	#legend("bottomright", c("Length-at-age", "Vulnerability", "Survival"), 
	#	lty=c(-1, 1, 1), pch=c(1, -1, -1),col=c(1, 1, 2),  bty="n")

}
#life.history()

retro<-function(retyr=0)
{    
	mcon<<-FALSE
	setwd("~/Documents/CONSULTING/HBC ASMR Performance/Lew's Code & Data/") 
	rt=matrix(nrow=length(B$rt2), ncol=retyr+1)
	nt4=matrix(nrow=length(B$adults), ncol=retyr+1)
	#arg <- paste("./ASMR3t -ind chubdat.dat")
	for(i in retyr:0)
	{
		arg <- paste("./ASMR3t -ind chubdat2009.dat -est -nox -retro",  i)
		system(arg)
		trt=get.results("ASMR3t.rep")$rt2
		rt[1:length(trt),i+1]=trt
		trt=get.results("ASMR3t.rep")$adults 
		nt4[1:length(trt), i+1]=trt
	}
	matplot(B$byr, rt, type="l")
	return(list(nt4=nt4, rt=rt))
}

run.retro<-function(dfile="S1.dat", retyr=2)
{	
	#setwd("~/Documents/CONSULTING/HBC ASMR Performance/ASRM 2.0")
	mcon<<-FALSE
	arg <- paste("./ASMR -ind",dfile,"-est -nox -retro", 0)
	system(arg)
	A=get.results(fn)
	nt4 = pt = matrix(nrow=length(A$nt4),  ncol=retyr+1)
	rt = matrix(nrow=length(A$rt), ncol=retyr+1)
	M = rep(0, length=retyr+1)
	for(i in retyr:0)
	{
		arg <- paste("./ASMR -ind",dfile,"-est -nox -retro", i)
		system(arg)
		if(i==0)
		{
			yr=get.results(fn)$yr
			rt.yrs=get.results(fn)$byr#c((yr[1]-10):(yr[1]-1),yr)
		}
		tmp = get.results(fn)$nt4
		nt4[1:(length(tmp)), i+1]=tmp
		tmp = get.results(fn)$rt
		rt[1:(length(tmp)), i+1]=tmp
		tmp = get.results(fn)$pt
		pt[1:(length(tmp)), i+1]=tmp
		M[i+1] = get.results(fn)$M
	}
	
	arg <- paste("./ASMR -ind",dfile,"-est -nox -retro", 0)
	system(arg)
	
	par(mfcol=c(2, 2), las=1)
	matplot(yr, nt4/1000, type="l", xlab="Year", ylab="Age 4+ numbers (1000s)", ylim=c(0, 12))
	matplot(rt.yrs-2, rt/1000, type="l", xlab="Brood year", ylab="Age-2 recruits (1000s)", ylim=c(0, 6))
	abline(v=1999)
	matplot(yr, pt, type="l", xlab="Year", ylab="Capture probability")
	barplot(rev(M), names.arg=seq(max(yr)-retyr, max(yr), length=retyr+1)
		,ylab="Natural mortality rate", xlab="Retrospective year")
	return(list(nt4=nt4, rt=rt))
}

run.mle.scenarios<-function()
{
	#Use this routine to run and compare the mle estimates
	#from the 7 alternative scenarios (data subsets)
	df=paste("S", 1:7, ".dat", sep="") 
	nt4=NULL
	rt=phi=va=NULL
	
	
	for(idf in df)
	{
		arg <- paste("./ASMR -ind", idf)
		print(arg)
		system(arg) 
		
		A=get.results(fn)
		nt4=cbind(nt4, A$nt4) 
		rt=cbind(rt, A$rt)
		tmp=(lowess(colMeans(A$pta[13:21,]),f=1/4))$y
		tmp[tmp<=0]=0
		tmp=tmp/max(tmp)
		va=cbind(va,tmp)
		
		F=read.fit("ASMR")
		est=F$est[13:16]
		std=F$std[13:16]
		phi = rbind(phi,c(est, std))

	} 
	par(las=1, cex=1.5)
	matplot(A$yr, nt4/1000, type="o", col=1, lty=1,xlab="Year",ylab="Adults (1000s age 4+)")
	dev.copy2eps(file="../Report/Fig.1.eps")
	matplot(A$yr, t(t(nt4)/nt4[1,]), type="o", col=1, lty=1,xlab="Year",ylab="Adult depletion")
	dev.copy2eps(file="../Report/Fig.1b.eps")
	
	
	matplot(A$byr,rt/1000,type="o", col=1, lty=1,xlab="Year",ylab="Age-2 recruits (1000s)") 
	dev.copy2eps(file="../Report/Fig.2.eps")
	write.table(phi, "Spars.txt")
	matplot(2:50,va, type="o",lty=1, xlab="Age", ylab="Relative age-specific selectivity (2000-2009)")
	dev.copy2eps(file="../Report/Fig.3.eps")
	write.table(va, file="Selectivity.dat", row.names=FALSE, col.names=FALSE)
	return(phi)
}

simulation<-function(seed=999)
{
	arg <- paste("./ASMR -est -sim", seed, "10")
	system(arg)
	A=get.results(fn)
	S=get.results("SimValues.dat")
	
	par(mfcol=c(2, 2), las=1)
	n=length(S$yr)
	ix=50:n
	plot(S$yr[ix],S$nt4[ix]/1000, type="l", xlab="Year", ylab="Abundance (1000 age 4+)", lwd=2, col="grey")
	lines(A$yr, A$nt4/1000)
	
	
	ix=40:n
	plot(S$yr[ix],S$rt[ix]/1000, type="l", xlab="Year", ylab="Abundance (1000 age-2 recruits)", lwd=2, col="grey")
	lines(S$yr[ix], A$rt/1000)
	
	#Now do retrospective runs.
	retyr=1:8
	for(i in retyr)
	{
		arg <- paste("./ASMR -est -ind simASMR.dat -retro", i)
		system(arg)
	}
}


run.model <- function(data.file="asmr.dat",mcmc=0)
{
	#This function calls the ASMR exe
	#for mcmc I would recommend 100,000 iterations
	arg <- paste("./ASMR -ind",data.file,"-nox")
	if(!mcmc){
		mcon <<- FALSE
	}
	if(mcmc){
		arg <- paste(arg,"-mcmc",format(mcmc,scientific=F),"-mcsave 50","-nosdmcmc")
		mcon <<- TRUE
	}
	system(arg)
	if(mcon) system(paste(arg,"-mceval"))
}




get.results <- function(fn)
{
	#This function reads the report file,
	#and if mcon=T reads the posterior samples
	A <- reptoRlist(fn)
	if(mcon)
	{
		A$mc = mc=read.table("asmr.mcmc",header=T)
		A$mc.nt = mc.nt=read.table("nt4.mcmc",header=F)
		A$mc.rt = mc.rt=read.table("rt.mcmc",header=F)
		A$mc.pt = mc.pt = read.table("pt.mcmc", header=F)
	}
	return(A)
}

plot.adults<-function()
{
	#Plot median age-4+ numbers and 95% CI
	nt=t(apply(A$mc.nt, 2, quantile, probs=c(0.025, 0.5, 0.95)))
	matplot(A$yr, nt/1000,type="l",  lty=c(2, 1, 2), col=1, xlab="Year", ylab="Adult abundance (1000s)")
}

plot.recruits<-function()
{
	#Plot median age-2 recruits and 95% CI 
	rt.yrs=c((A$yr[1]-10):(A$yr[1]-1),A$yr)
	rt=t(apply(A$mc.rt, 2, quantile, probs=c(0.025, 0.5, 0.95)))
	matplot(rt.yrs, rt/1000,type="l",  lty=c(2, 1, 2), col=1, xlab="Year", ylab="Age-2 recruits (1000s)")
	boxplot( (A$mc.rt)/1000, names.arg=rt.yrs, ylim=c(0, 7), pch="")
}



plot.summary=function(A,saveplots=T,figfile="ASMRsummary.pdf")
{
	
	pl1 <- function(){
		par(las=1,ps=12,cex.main=14/12)
		rt.yrs=c((A$yr[1]-10):(A$yr[1]-1),A$yr)
		plot(A$yr,A$nt4/1000,type="o",ylim=c(0,max(A$nt2)/1000),xlab="Year",ylab="Abundance (4+ numbers)",main="Abundance")
		#lines(A$yr,A$nt2,col=2)
		#lines(A$yr,A$nt3,col=3)
		
		plot(A$yr,A$pt,type="l",ylim=c(0,1),xlab="Year",ylab="Capture probability")
		plot(1987:2004,A$rt[9:26]/1000,type="o", pch=18,xlab="Year",ylab="Age-2 recruits (1000)", ylim=c(0, max(A$rt)/1000))

		#raw data
		bubble.plot(A$yr,A$age[1:30],A$mta[, 1:30],0.051,log.scale=T,main="Logarithm of New Marks Released")
		bubble.plot(A$yr,A$age[1:31],A$rta[, 1:31],0.051,log.scale=T,main="Logarithm of Marks Recaptured")

		#residuals for new marks
		bubble.plot(A$yr,A$age,A$epsilon,0.2,log.scale=F,ylim=c(0,30),main="Residuals for Marks released", leg=T)

		#residuals for recaps
		bubble.plot(A$yr,A$age,A$delta,0.2,log.scale=F,ylim=c(0,30),main="Residuals for Recaptured Marks", leg=T)

		## MCMC stuff
	
		if(with(A,exists("mc"))){
			pairs(A$mc,pch=".",main="Posterior samples",gap=0,upper.panel=NULL,diag.panel=panel.hist)

			#mc.nt=read.table("nt2.mcmc",header=F)
			ci = apply(A$mc.nt,2,quantile,probs=c(0.025,0.5,0.975))
			matplot(A$yr,t(ci),type="l",xlab="Year",ylab="4+ numbers",lty=c(2,1,2),col=1,ylim=c(0,max(ci)))

			ci = apply(A$mc.rt,2,quantile,probs=c(0.025,0.5,0.975))
			matplot(rt.yrs,t(ci),type="l",xlab="Year",ylab="Age-1 recruits",lty=c(2,1,2),col=1)
			
			#posterior and prior for M
			d.m = density(A$mc$m, adjust=1.5)
			plot(d.m, xlab="M", ylab="Marginal posterior density", xlim=c(0, 0.25), main="")
			x=seq(0, 0.25, by=d.m$bw)
			lines(x, dlnorm(x, log(0.13), 0.2), col="grey", lwd=2)
		}
	}
	pl1();
	if(saveplots)
	{
		pdf(file=figfile)
		pl1()
		dev.off()
	}

}


## __MAIN SECTION__

#A=get.results(fn)
#plot.summary(A, F) 
#plot.adults() 
#plot.recruits()


## _________________