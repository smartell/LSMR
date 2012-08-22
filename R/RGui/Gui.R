.VIEWTRCK <- "./lsmrTracker.txt"
.WINTEXT  <- "./lsmrWin.txt"
.OVERLAY  <- FALSE



guiView <-
function( winName )
{
	# Requires libraries
	require(PBSmodelling)
	
	# Close any open graphics devices or X11 devices
	graphics.off()
	closeWin()
	
	# Check for global .VIEWTRCK and read header files
	trckExists <- file.exists( .VIEWTRCK )
	if (trckExists)
	{
		tmp <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
		cat( "MSG (.lsmrView): Viewer tracking file ",.VIEWTRCK, " found.\n" )
		ifiles <- tmp
		gomenu <- TRUE
	}
	else
	{
		cat( "Error MSG (.iscamView): missing", .VIEWTRCK, "file.\n")
		gomenu <- FALSE
	}
	
	if(gomenu)
	createWin(.WINTEXT)
	
}


.viewPlot <-
function()
{
	print(".mpdView")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	# List of gui changes
	guiChanges <- list()
	
	# Determine which files have been selected
	# and read files and store into model object (M)
	hdr	<- ifiles[ ifiles$Select, ]
	fn	<- hdr$Control.File
	M	<- lapply(fn, read.admb)
	names(M) <- hdr$Model
	
	# Overlay
	if(overLay)
	{
		.OVERLAY <<- TRUE
	}
	else
	{
		.OVERLAY <<- FALSE
	}
	
	# use plotType to determine which function to call.
	switch(plotType, 
		abun50={
			ylbl <- "Abundance (> 50 mm)"
			.lsmrPlotNt(M, what="Nt", ylbl=ylbl)
		}, 
		abun100={
			ylbl <- "Abundance (> 100 mm)"
			.lsmrPlotNt(M, what="N100", ylbl=ylbl)
		}, 
		abun150={
			ylbl <- "Abundance (> 150 mm)"
			.lsmrPlotNt(M, what="N150", ylbl=ylbl)
		}, 
		agecomp={
			print("agecomp")
			.plotAgeComps(M)
		},
		agehist={
			print("agehist")
			.plotAgeHist(M)
		},  
		meanwt={
			print("meanwt")
			.plotMeanWt(M)
		}, 
		vbio={
			print("vbio")
			.plot_bt(M, "bt")
		}, 
		sbio={
			print("sbio")
			.plot_bt(M, "sbt")
		}, 
		urate={
			print("urate")
			.plot_bt(M, "ut")
		}, 
		ct_res={
			print("ct_res")
			.plot_resid(M, "eta")
		}, 
		it_res={
			print("it_res")
			.plot_resid(M, "epsilon")
		}
		)
	
	# Return a global model object that are in play
	iMod <<- M
}

