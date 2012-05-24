# R-script for creating summary data for LSMR based on the data-based extraction from GCMRC

dfile    <- "raw2012data.csv"
# ------------------------------------------------------------------------------- #
# Read data file and create data.frame object for use. Returns (DF)
# ------------------------------------------------------------------------------- #
get.DF   <- function(dfile)
{
	DF       <- read.csv(dfile, header=T)
	t1       <- as.POSIXct(strptime(as.character(DF$START_DATETIME), "%m/%d/%Y"))
	DF$TAGNO <- DF$TH_ENCOUNTER_RANKING
	DF$TL    <- DF$TOTAL_LENGTH
	
	return (DF)
}

# ------------------------------------------------------------------------------- #
# Construct growth increment data.frame for estimating growth using Laslett model.
# ------------------------------------------------------------------------------- #


if(!exists("DF")) DF <- get.DF(dfile)