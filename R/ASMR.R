# --------------------------------------------------------------------------- #
# R-script for plotting output from ASMR.rep
# Author: Steven Martell
# DATE: March 27, 2012
# LOCATION: Seattle Washington
# DEPENDENCIES: read.admb.R






# --------------------------------------------------------------------------- #

source("read.admb.R", echo=FALSE)


.FILENAME	<- "../ADMB/srcASMR/asmr"

obj			<- read.admb(.FILENAME)
class(obj)	<- c(class(obj), "asmr")



# --------------------------------------------------------------------------- #
# S3 method for class 'asmr'                                                  #
# --------------------------------------------------------------------------- #
plot.asmr <- function(obj, ask=dev.interactive(), ...)
{
	devAskNewPage(ask)
	with(obj,{
		plot(yr, nt4, type="l", xlab="Year", ylab="Abundance (age-4+)"
			,ylim=c(0, 1.2*max(nt4)),  ...)
			
		plot(byr, rt, type="l", xlab="Brood year", ylab="Age-2 recruits", 
			,ylim=c(0, 1.2*max(rt)), ...)
	})
}

print.asmr <- function(obj, ...)
{
	print(attributes(obj))
}
